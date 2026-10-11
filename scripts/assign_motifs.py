"""Copy a CisBP PWM to each target TF gene by Similarity Regression, and write one MEME file.

Usage: python scripts/assign_motifs.py --target tatroviride_p1.toml [--log_path logs/assign_motifs/<label>.log]

Needs the reference (build_reference.py) and the TF calls of the target (call_tfs.py).

Each target TF's SR segments are scored against every usable reference TF (a CisBP TF
with a lab-measured PWM) of the same CisBP family, with CisBP's SR model for that
family. Tiers:
  pass             some reference TF is highly similar (SR's HSim: score at or above the
                   family cutoff, or identical DBDs) and has a PWM; the PWM of the best
                   such reference TF is copied
  below-threshold  reference TFs were scored but none qualifies; no PWM
  none             nothing could be scored (no family, no segment, or no reference TF); no PWM
A copied PWM is a prediction of binding from DBD similarity, not evidence of it.

Writes results/<label>/assignments.tsv (one row per target TF gene) and motifs.meme
(each copied PWM once, taken from the reference's motifs.meme).
"""

import argparse
import logging

import pandas as pd

from srmotifs import files, meme, reference, sr
from srmotifs.constants import (ASSIGNMENTS_TSV, FAMILIES_TSV, LOGS_DIR, MOTIF_SOURCE_COLUMNS, MOTIFS_MEME,
                                REFERENCES_TSV, SR_DIR, SR_HSIM, TF_CALLS_TSV, TIER_NONE, TIER_PASS)
from srmotifs.log import configure_logging

COLUMNS = ["gene_id", "protein_id", "family_id", "family_name", "tier", "reason", "motif_id",
           "best_reference_tf", "best_reference_name", "best_reference_species", "sr_score", "cutoff", "margin",
           "sr_class", "pct_id", "n_segments", "n_references", "n_hsim", "runner_up_tf", "runner_up_score",
           "motif_experiment", "motif_study", "motif_year", "motif_pmid"]

logger = logging.getLogger(__name__)


def unscorable_reason(gene, references):
    """Why a gene cannot be scored, or None if it can."""
    if not gene.family_id:
        return "no CisBP family for its DBD types"
    if not gene.n_segments:
        return "no DBD segment at the segment threshold"
    if not references:
        return "no reference TF in family"
    return None


def assign_one(gene, references, cutoff, motif_of, scorer):
    """Score one target TF gene against its family's reference TFs and decide its tier."""
    scores = scorer.score_all(gene.segments.split(","), gene.family_id, references)
    tier, reason, best, runner_up = sr.decide(scores, motif_of)

    return {"tier": tier, "reason": reason,
            "motif_id": motif_of[best["reference_tf"]] if tier == TIER_PASS else "",
            "best_reference_tf": best["reference_tf"], "sr_score": best["sr_score"],
            "cutoff": cutoff, "margin": best["sr_score"] - cutoff,
            "sr_class": best["sr_class"], "pct_id": best["pct_id"],
            "n_references": len(references), "n_hsim": int((scores["sr_class"] == SR_HSIM).sum()),
            "runner_up_tf": "" if runner_up is None else runner_up["reference_tf"],
            "runner_up_score": None if runner_up is None else runner_up["sr_score"]}


def assign(genes, by_family, families, motif_of, scorer):
    rows = []
    for gene in genes.itertuples():
        row = {"gene_id": gene.gene_id, "protein_id": gene.protein_id, "family_id": gene.family_id,
               "family_name": gene.family_name, "n_segments": gene.n_segments}
        references = by_family.get(gene.family_id, [])

        reason = unscorable_reason(gene, references)
        if reason:
            row.update({"tier": TIER_NONE, "reason": reason})
        else:
            cutoff = families.loc[gene.family_id, "cutoff"]
            row.update(assign_one(gene, references, cutoff, motif_of, scorer))
        rows.append(row)

    return pd.DataFrame(rows)


def add_reference_details(assignments, references):
    """Name and species of the best reference TF; source experiment of the copied PWM."""
    info = references.set_index("TF_ID")
    assignments["best_reference_name"] = assignments["best_reference_tf"].map(info["TF_Name"])
    assignments["best_reference_species"] = assignments["best_reference_tf"].map(info["TF_Species"])

    passed = assignments["tier"] == TIER_PASS
    for column in MOTIF_SOURCE_COLUMNS:
        assignments[column] = assignments["best_reference_tf"].map(info[column]).where(passed, "")

    return assignments.astype({"n_references": "Int64", "n_hsim": "Int64"})[COLUMNS]


def check_partition(assignments, genes):
    """Every TF gene appears once, every scored gene has a score, and only pass genes have a motif."""
    if sorted(assignments["gene_id"]) != sorted(genes["gene_id"]) or assignments["gene_id"].duplicated().any():
        raise AssertionError("tiers do not cover every TF gene exactly once")

    scored = assignments["tier"] != TIER_NONE
    if assignments.loc[scored, "sr_score"].isna().any():
        raise AssertionError("a scored TF lacks a numeric SR score")

    has_motif = assignments["motif_id"].fillna("") != ""
    if not (has_motif == (assignments["tier"] == TIER_PASS)).all():
        raise AssertionError("motifs must be given to exactly the pass tier")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--target", required=True, help="target config, a path or a file name in config/")
    parser.add_argument("--log_path", default=None, help="default: logs/assign_motifs/<label>.log")
    args = parser.parse_args()

    config = files.reference_config()
    target = files.load_toml(args.target)
    label = target["organism"]["label"]
    configure_logging(args.log_path or LOGS_DIR / "assign_motifs" / f"{label}.log")
    logger.info(f"assigning PWMs to the TFs of {target['organism']['name']} ({label})")

    ref_dir = files.reference_dir(config)
    out = files.results_dir(label)

    # loading the reference, the TF calls and the SR models
    references, families, _, by_family, motif_of = reference.load(ref_dir)
    genes = pd.read_csv(out / TF_CALLS_TSV, sep="\t", dtype=str, keep_default_na=False)
    genes["n_segments"] = genes["n_segments"].astype(int)
    (sr_root,) = (ref_dir / SR_DIR).iterdir()
    scorer = sr.Scorer(sr_root, sr.family_models(sr_root, families))

    # scoring every TF against its family's references
    logger.info(f"scoring {len(genes)} TF genes against their families' reference TFs")
    assignments = assign(genes, by_family, families, motif_of, scorer)
    assignments = add_reference_details(assignments, references)
    check_partition(assignments, genes)
    assignments.to_csv(out / ASSIGNMENTS_TSV, sep="\t", index=False)

    # one MEME file with each copied PWM once
    blocks = meme.read_blocks(ref_dir / MOTIFS_MEME)
    copied = sorted(set(assignments["motif_id"].dropna()) - {""})
    meme.write(out / MOTIFS_MEME, {motif_id: blocks[motif_id] for motif_id in copied})

    inputs = [out / TF_CALLS_TSV, ref_dir / REFERENCES_TSV, ref_dir / MOTIFS_MEME, ref_dir / FAMILIES_TSV]
    files.write_provenance(out, inputs, name="provenance_assign_motifs")

    passed = assignments[assignments["tier"] == TIER_PASS]
    logger.info(f"TF genes: {len(assignments)}; tiers: {assignments['tier'].value_counts().to_dict()}")
    logger.info(f"{MOTIFS_MEME}: {len(copied)} distinct PWMs for {len(passed)} TFs")
    logger.info(f"pass by family: {passed['family_name'].value_counts().to_dict()}")
    logger.info(f"source species: {passed['best_reference_species'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
