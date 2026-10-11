"""Run the unmodified SimilarityRegression (SR) scoring code under Python 3, and assign tiers.

SR is fetched as a tarball of a pinned commit and extracted under reference/. Two
things stop the 2019 package from importing as-is:
- PairwiseAlignment.py imports `Bio.SubsMat.MatrixInfo.blosum62`, removed in
  Biopython 1.80. A stand-in module exposing the same BLOSUM62 values from
  `Bio.Align.substitution_matrices` is registered before importing it.
- PredictSimilarity.py mixes tabs and spaces in a helper that is never called, a
  TabError in Python 3. The two functions needed, ReadSRModel and
  ScoreAlignmentResult, come before that helper, so the file's text up to it is
  executed verbatim.
"""

import json
import sys
import tarfile
import types
from pathlib import Path

import pandas as pd
from Bio.Align import substitution_matrices

PREDICT_SIMILARITY_CUT = "def SeqDictIterator_XoverY"

TIER_PASS = "pass"
TIER_BELOW = "below-threshold"
TIER_NONE = "none"


def extract(tarball, dest):
    """Unpack the SR tarball into `dest` once; return the repository root inside it."""
    dest = Path(dest)
    if not dest.exists():
        partial = dest.with_suffix(".part")
        with tarfile.open(tarball) as archive:
            archive.extractall(partial, filter="data")
        partial.rename(dest)

    (root,) = dest.iterdir()
    return root


def install_blosum62_shim():
    """Register a module that looks like the old Bio.SubsMat.MatrixInfo, holding BLOSUM62."""
    if "Bio.SubsMat.MatrixInfo" in sys.modules:
        return

    matrix = substitution_matrices.load("BLOSUM62")
    matrix_info = types.ModuleType("Bio.SubsMat.MatrixInfo")
    matrix_info.blosum62 = {(a, b): int(matrix[a][b]) for a in matrix.alphabet for b in matrix.alphabet}

    package = types.ModuleType("Bio.SubsMat")
    package.MatrixInfo = matrix_info
    sys.modules["Bio.SubsMat"] = package
    sys.modules["Bio.SubsMat.MatrixInfo"] = matrix_info


def load(sr_root):
    """(PairwiseAlignment module, namespace holding ReadSRModel and ScoreAlignmentResult)."""
    install_blosum62_shim()
    if str(sr_root) not in sys.path:
        sys.path.insert(0, str(sr_root))
    from similarityregression import PairwiseAlignment

    # running PredictSimilarity.py up to the helper that does not compile under Python 3
    source_path = Path(sr_root) / "similarityregression" / "PredictSimilarity.py"
    source = source_path.read_text()
    if source.count(PREDICT_SIMILARITY_CUT) != 1:
        raise RuntimeError(f"unexpected layout of {source_path}")

    namespace = {"__name__": "similarityregression.PredictSimilarity_prefix"}
    prefix = source.split(PREDICT_SIMILARITY_CUT)[0]
    exec(compile(prefix, str(source_path), "exec"), namespace)
    return PairwiseAlignment, namespace


def family_models(sr_root, families):
    """{Family_ID: SR model} for a families table (columns Family_Name, sr_model, cutoff).

    `sr_model` names a file in SR's SRModels/. Each model's HSim threshold must equal
    the family cutoff CisBP records.
    """
    _, predict = load(sr_root)

    models = {}
    for family_id, row in families.iterrows():
        path = Path(sr_root) / "SRModels" / f"{row['sr_model']}.json"

        # a family-specific model must really be for this family
        if row["sr_model"] != "NO_THRESHOLD":
            if json.loads(path.read_text())["Family_Name"] != row["Family_Name"]:
                raise ValueError(f"{path.name} is not family {row['Family_Name']}")

        model = predict["ReadSRModel"](str(path))
        if abs(model["Threshold.HSim"] - row["cutoff"]) > 1e-9:
            raise ValueError(f"{family_id}: SR threshold {model['Threshold.HSim']} "
                             f"!= CisBP cutoff {row['cutoff']}")
        models[family_id] = model

    return models


class Scorer:
    def __init__(self, sr_root, models):
        self.alignment, predict = load(sr_root)
        self._score = predict["ScoreAlignmentResult"]
        self.models = models

    def score_all(self, query_segments, family_id, references):
        """One row per (reference TF id, segments): SR score, SR class (HSim/Amb/Dis), alignment %ID."""
        rows = []
        for reference_id, reference_segments in references:
            result = self.alignment.AlignDBDArrays(("query", list(query_segments)),
                                                   ("reference", list(reference_segments)))
            score, sr_class = self._score(resultDict=result, scoreDict=self.models[family_id])
            rows.append((reference_id, float(score), sr_class, float(result["PctID_L"])))

        return pd.DataFrame(rows, columns=["reference_tf", "sr_score", "sr_class", "pct_id"])


def decide(scores, motifs):
    """(tier, reason, best row, runner-up row or None) for one target TF's scores against its family.

    `motifs` maps reference TF -> the PWM copied from it. pass needs a reference TF that
    SR calls highly similar (HSim) and that has a PWM; the best is the highest SR score
    among those (then %ID, then id), and its PWM is copied. Otherwise the best is the
    most similar reference TF overall, and no PWM is copied.
    """
    if scores.empty:
        return TIER_NONE, "no reference TF in family", None, None

    ranked = scores.sort_values(["sr_score", "pct_id", "reference_tf"], ascending=[False, False, True])
    hsim = ranked["sr_class"] == "HSim"
    with_pwm = ranked[hsim & ranked["reference_tf"].isin(motifs)]

    if len(with_pwm):
        tier, reason, best = TIER_PASS, "highly similar reference TF with a PWM", with_pwm.iloc[0]
    elif hsim.any():
        tier, reason, best = TIER_BELOW, "highly similar reference TFs have no PWM matrix", ranked.iloc[0]
    else:
        tier, reason, best = TIER_BELOW, "no reference TF is highly similar", ranked.iloc[0]

    others = ranked[ranked["reference_tf"] != best["reference_tf"]]
    runner_up = others.iloc[0] if len(others) else None
    return tier, reason, best, runner_up
