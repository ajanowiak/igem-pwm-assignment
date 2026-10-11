"""Check that this pipeline reproduces SR and CisBP's own SR inferences.

Usage: python scripts/check_sr_vs_cisbp.py [--log_path logs/check_sr_vs_cisbp.log]

Needs the reference built by build_reference.py. Optional: the pipeline does not use its output.

1. SR's shipped example: the SR code, run here, gives the 55 pairwise scores saved in
   SR's example notebook, and hmmalign (with the semiglobal end placement in
   srmotifs/hmmer.py) and the reference's Pfam Homeobox model gives the same
   match-state strings as SR's own aphid alignments.
2. CisBP's own PWM copies: CisBP marks each TF without a lab-measured PWM as I (CisBP
   copied a PWM to it by SR) or N (no reference TF cleared the family cutoff), and names
   the copied PWM for I. For the species in config [check], those TFs are resolved from
   CisBP's own protein table exactly like reference TFs, scored, and compared.
Agreement is evidence that the pipeline implements CisBP's method. It is not evidence
that any inferred motif is biologically correct.

Writes results/checks/sr_vs_cisbp/{summary.txt, comparison.tsv, agreement_by_margin.tsv}.
"""

import argparse
import json
import logging
import tempfile
from pathlib import Path

import pandas as pd

from srmotifs import cisbp, files, hmmer, reference, sr
from srmotifs.constants import (CHECKS_DIR, CISBP_TF_TABLE, DBD_MODELS_HMM, LOGS_DIR, SR_DIR, SR_HSIM,
                                TIER_PASS)
from srmotifs.log import configure_logging

MARGIN_BINS = [0.0, 0.025, 0.05, 0.1, 0.2, float("inf")]

logger = logging.getLogger(__name__)


def shipped_matchpos(sr_root):
    """{DBD sequence: match-state string} as SR's own aphid alignment produced them."""
    path = sr_root / "Example" / "DBDMatchPos_aphid" / "HDSeqs.matchpos_semiglobal.fa"
    lines = path.read_text().split()
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines), 2)}


def notebook_scores(sr_root):
    """{frozenset(motif pair): (PctID_L, SR score, class)} from the saved output of SR's example notebook."""
    notebook = json.loads((sr_root / "Example" / "Example Analysis Notebook.ipynb").read_text())
    cell = next(c for c in notebook["cells"] if "ReadSRModel" in "".join(c["source"]))
    table = "".join(cell["outputs"][0]["data"]["text/plain"]).splitlines()[1:]

    # pandas printed the first motif of each pair only on the first row of its group
    expected, first = {}, None
    for row in table:
        parts = row.split()
        if len(parts) == 7:
            first, parts = parts[0], parts[1:]
        expected[frozenset((first, parts[0]))] = (float(parts[1]), float(parts[3]), parts[4])
    return expected


def example_segments(sr_root):
    """{motif id: SR segments} for SR's homeodomain example, built as its notebook does."""
    data = sr_root / "Example" / "HomeodomainData"
    matchpos = shipped_matchpos(sr_root)
    proteins = pd.read_csv(data / "prot_seq.txt", sep="\t", index_col=0)
    first_dbds = proteins.groupby(level=0)["DBD_seqs"].first()

    segments = {}
    for line in (data / "TF_Information.txt").read_text().splitlines()[1:]:
        fields = line.split("\t")
        for motif in (m.strip().rstrip(",") for m in fields[6:]):
            if motif:
                segments[motif] = [matchpos[dbd] for dbd in first_dbds[fields[0]].split(",")]
    return segments


def check_sr_example(sr_root):
    """(largest difference from the notebook's scores and %IDs, class mismatches, pairs compared)."""
    alignment, predict = sr.load(sr_root)
    model = predict["ReadSRModel"](str(sr_root / "SRModels" / "F223_1.97d.json"))
    segments = example_segments(sr_root)
    expected = notebook_scores(sr_root)

    worst, mismatches = 0.0, 0
    for pair, (pct_id, score, sr_class) in expected.items():
        first, second = sorted(pair)
        result = alignment.AlignDBDArrays((first, segments[first]), (second, segments[second]))
        ours, our_class = predict["ScoreAlignmentResult"](resultDict=result, scoreDict=model)
        worst = max(worst, abs(ours - score), abs(result["PctID_L"] - pct_id))
        mismatches += our_class != sr_class

    return worst, mismatches, len(expected)


def check_match_states(hmm_path, model_length, sr_root):
    """(SR example DBDs whose hmmalign match states equal SR's aphid strings, total)."""
    shipped = shipped_matchpos(sr_root)
    sequences = {f"seq{i}": seq for i, seq in enumerate(sorted(shipped))}
    ours = hmmer.align_to_model(hmm_path, sequences, model_length)

    identical = sum(ours[name] == shipped[seq] for name, seq in sequences.items())
    return identical, len(sequences)


def score_targets(targets, resolved, by_family, motif_of, scorer):
    """One row per CisBP I/N TF with segments: CisBP's call next to ours."""
    status = targets.set_index("TF_ID")

    rows = []
    for row in resolved[resolved["segments"] != ""].itertuples():
        references = by_family.get(row.Family_ID, [])
        scores = scorer.score_all(row.segments.split(","), row.Family_ID, references)
        tier, _, best, _ = sr.decide(scores, motif_of)

        source = status.loc[row.TF_ID, "cisbp_source_tf"]
        hsim = set(scores.loc[scores["sr_class"] == SR_HSIM, "reference_tf"])
        rows.append({"TF_ID": row.TF_ID, "TF_Species": row.TF_Species, "Family_Name": row.Family_Name,
                     "cisbp_status": status.loc[row.TF_ID, "TF_Status"], "cisbp_source_tf": source,
                     "cisbp_source_usable": source in {tf_id for tf_id, _ in references},
                     "our_tier": tier,
                     "our_best_tf": "" if best is None else best["reference_tf"],
                     "our_best_score": None if best is None else best["sr_score"],
                     "cisbp_source_in_our_hsim": source in hsim})

    return pd.DataFrame(rows)


def summarise(table):
    i_rows = table[table["cisbp_status"] == "I"]
    n_rows = table[table["cisbp_status"] == "N"]
    lines = [pd.crosstab(table["cisbp_status"], table["our_tier"]).to_string()]

    if len(i_rows):
        share = (i_rows["our_tier"] == TIER_PASS).mean()
        lines.append(f"CisBP I -> our pass: {share:.1%} of {len(i_rows)}")
    if len(n_rows):
        share = (n_rows["our_tier"] != TIER_PASS).mean()
        lines.append(f"CisBP N -> our no pass: {share:.1%} of {len(n_rows)}")

    # where both infer a motif, do we pick CisBP's source TF?
    both = i_rows[(i_rows["our_tier"] == TIER_PASS) & i_rows["cisbp_source_usable"]]
    if len(both):
        same_source = (both["our_best_tf"] == both["cisbp_source_tf"]).mean()
        lines.append(f"both pass, CisBP's source usable ({len(both)}): our best hit is CisBP's source "
                     f"{same_source:.1%}; CisBP's source is HSim for us {both['cisbp_source_in_our_hsim'].mean():.1%}")
    return "\n".join(lines)


def agreement_by_margin(table, families):
    """Share of our pass calls that CisBP also infers (I), by SR score margin above the family cutoff."""
    passed = table[table["our_tier"] == TIER_PASS].copy()
    cutoff = families.set_index("Family_Name")["cutoff"]
    margin = passed["our_best_score"] - passed["Family_Name"].map(cutoff)
    passed["bin"] = pd.cut(margin, MARGIN_BINS, right=False)

    return passed.groupby("bin", observed=False).agg(
        n_pass=("TF_ID", "size"),
        cisbp_agrees=("cisbp_status", lambda status: round((status == "I").mean(), 4)))


def sr_example_lines(sr_root, ref_dir, models, config):
    """Part 1: report lines, and whether the example is reproduced exactly."""
    worst, mismatches, n_pairs = check_sr_example(sr_root)

    with tempfile.TemporaryDirectory() as tmp:
        homeobox = Path(tmp) / "Homeobox.hmm"
        homeobox.write_text(files.run(["hmmfetch", str(ref_dir / DBD_MODELS_HMM), "Homeobox"]))
        same, total = check_match_states(homeobox, models["Homeobox"][0], sr_root)

    lines = ["== 1. SR's shipped example",
             f"SR scores vs SR's notebook: {n_pairs} pairs, max |difference| {worst:.1e}, "
             f"class mismatches {mismatches}",
             f"hmmalign match states vs SR's aphid strings (Pfam {config['pfam']['release']} Homeobox): "
             f"{same}/{total} identical"]
    passed = worst <= 1e-6 and not mismatches and same == total
    return lines, passed


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--log_path", default=LOGS_DIR / "check_sr_vs_cisbp.log")
    args = parser.parse_args()
    configure_logging(args.log_path)

    config = files.reference_config()
    ref_dir = files.reference_dir(config)
    CHECKS_DIR.mkdir(parents=True, exist_ok=True)
    _, families, models, by_family, motif_of = reference.load(ref_dir)
    (sr_root,) = (ref_dir / SR_DIR).iterdir()

    # part 1: SR's own example
    logger.info("reproducing SR's shipped example")
    lines, example_passed = sr_example_lines(sr_root, ref_dir, models, config)

    # part 2: resolving CisBP's I/N TFs like references, and scoring them
    logger.info("scoring CisBP's own I/N TFs of the check species")
    tf_table = cisbp.load_tf_table(cisbp.raw_dir(config) / CISBP_TF_TABLE)
    targets = cisbp.check_targets(tf_table, config["check"]["species"])
    resolved, _ = reference.resolve_tfs(targets, cisbp.load_proteins(config), families,
                                        ref_dir / DBD_MODELS_HMM, models,
                                        config["rules"]["segment_max_i_evalue"])
    scorer = sr.Scorer(sr_root, sr.family_models(sr_root, families))
    table = score_targets(targets, resolved, by_family, motif_of, scorer)
    agreement = agreement_by_margin(table, families.reset_index())

    table.to_csv(CHECKS_DIR / "comparison.tsv", sep="\t", index=False)
    agreement.to_csv(CHECKS_DIR / "agreement_by_margin.tsv", sep="\t")

    # writing the summary, overall and per species
    lines += ["", "== 2. CisBP's own inferences",
              f"species: {', '.join(config['check']['species'])}",
              f"CisBP I/N TFs: {len(targets)}; with SR segments: {(resolved['segments'] != '').sum()}",
              summarise(table), "",
              "agreement of our pass calls with CisBP, by margin above the cutoff:", agreement.to_string()]
    for species, group in table.groupby("TF_Species"):
        lines += ["", f"-- {species}", summarise(group)]

    summary = "\n".join(lines)
    (CHECKS_DIR / "summary.txt").write_text(summary + "\n")
    logger.info(f"summary, also written to {CHECKS_DIR / 'summary.txt'}:\n{summary}")
    if not example_passed:
        logger.error("SR example check failed")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
