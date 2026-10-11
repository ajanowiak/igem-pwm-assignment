"""Find the transcription factors (TFs) in a proteome.

Usage: python scripts/call_tfs.py --target tatroviride_p1.toml [--log_path logs/call_tfs/<label>.log]

Needs the reference built by build_reference.py.

What it does, in order:
1. reads the proteins: translated from genome + GFF3 with gffread, or read from a
   protein FASTA (srmotifs/proteome.py)
2. scans every protein with hmmscan against the Pfam models of the DNA-binding domain
   (DBD) types CisBP recognises; hmmscan only reports which Pfam domain matched where
3. calls a protein a TF if one of its DBD hits passes Pfam's own score cutoff for that
   domain (the "gathering threshold") and is a DBD type CisBP uses to define TFs;
   where two hits overlap, the better one is kept. Optionally (`calling_i_evalue` in
   config/reference.toml), genes with no such hit may be called from weaker hits;
   the default 0 keeps the gathering threshold only
4. keeps one protein per gene: among the gene's isoforms with DBDs, the one with the
   most DBDs (ties: highest domain score, then longest, then smallest id)
5. gives each TF gene its CisBP family: the family defined by exactly the same set of
   DBD types (one homeobox -> Homeodomain); no matching family -> no family
6. turns that protein's DBDs into SR segments, exactly as for the reference TFs

Writes results/<label>/tf_calls.tsv (one row per TF gene) and dbd_hits.tsv (every hit
hmmscan reported, for inspection).
"""

import argparse
import logging

import pandas as pd

from srmotifs import files, hmmer, proteome, reference
from srmotifs.constants import (DBD_HITS_TSV, DBD_MODELS_HMM, DBD_MODELS_TSV, FAMILIES_TSV, LOGS_DIR,
                                TF_CALLS_TSV)
from srmotifs.log import configure_logging

COLUMNS = ["gene_id", "protein_id", "n_dbd_isoforms", "protein_length", "n_dbds", "architecture",
           "dbd_types", "family_id", "family_name", "n_segments", "segments"]

logger = logging.getLogger(__name__)


def call_domains(hits, models, calling_types, proteins, calling_i_evalue):
    """DBD hits that make a protein a TF, of DBD types that call TFs, overlaps resolved.

    A hit counts if it passes Pfam's gathering threshold (GA). If `calling_i_evalue` is
    above 0, genes with no GA hit at all may also be called from hits at or below that
    i-Evalue. Genes with a GA hit keep only their GA hits, so a looser value adds genes
    but never changes the domains, and so the family, of a gene GA already calls.
    """
    hits = hits[hits["model"].isin(calling_types)]
    passing = hmmer.gathering_pass(hits, models)
    if not calling_i_evalue:
        return hmmer.resolve_overlaps(passing)

    # genes GA already calls are left exactly as they are
    gene_of = proteins.set_index("protein_id")["gene_id"]
    called_genes = set(passing["protein"].map(gene_of))
    uncalled = ~hits["protein"].map(gene_of).isin(called_genes)

    # other genes may be called from hits at or below the i-Evalue cutoff
    loose = hits[uncalled & (hits["i_evalue"] <= calling_i_evalue)]
    return hmmer.resolve_overlaps(pd.concat([passing, loose]))


def representatives(called, proteins):
    """One row per TF gene: its representative protein and that protein's DBDs."""
    # summarising each TF protein's domains
    summary = called.groupby("protein").agg(
        n_dbds=("model", "size"),
        architecture=("model", ",".join),
        dbd_types=("model", lambda models: ",".join(sorted(set(models)))),
        dbd_score_sum=("dom_score", "sum"),
    ).reset_index()

    info = proteins.set_index("protein_id")
    summary["gene_id"] = summary["protein"].map(info["gene_id"])
    summary["protein_length"] = summary["protein"].map(info["sequence"].str.len())
    summary["n_dbd_isoforms"] = summary.groupby("gene_id")["protein"].transform("size")

    # picking one protein per gene
    summary = summary.sort_values(["gene_id", "n_dbds", "dbd_score_sum", "protein_length", "protein"],
                                  ascending=[True, False, False, False, True])
    genes = summary.drop_duplicates("gene_id").rename(columns={"protein": "protein_id"})
    return genes.reset_index(drop=True)


def add_families(genes, families):
    """The CisBP family whose DBD types are exactly the gene's; empty if no family has them."""
    family_of_types = {frozenset(row.DBDs.split(",")): family_id for family_id, row in families.iterrows()}
    genes["family_id"] = [family_of_types.get(frozenset(types.split(",")), "") for types in genes["dbd_types"]]
    genes["family_name"] = genes["family_id"].map(families["Family_Name"]).fillna("")
    return genes


def add_segments(genes, hits, sequences, families, ref_dir, models, max_i_evalue):
    """SR segments of each representative protein, by the same rule as for references."""
    family_dbds = {gene.protein_id: families.loc[gene.family_id, "DBDs"].split(",")
                   for gene in genes.itertuples() if gene.family_id}

    domains = hmmer.segment_hits(hits[hits["protein"].isin(family_dbds)], family_dbds, max_i_evalue)
    segment_map = hmmer.segments(sequences, domains, family_dbds, ref_dir / DBD_MODELS_HMM, models)

    genes["n_segments"] = [len(segment_map.get(protein, [])) for protein in genes["protein_id"]]
    genes["segments"] = [",".join(segment_map.get(protein, [])) for protein in genes["protein_id"]]
    return genes


def report(genes, notes):
    logger.info(f"proteins: {notes['proteins']} ({notes['genes']} genes; dropped for an internal stop: "
                f"{notes['dropped_internal_stop']})")
    if notes["proteins_without_gene_id"]:
        logger.warning(f"{notes['proteins_without_gene_id']} proteins have no gene id and are each "
                       "counted as their own gene; isoforms of one gene may be counted twice")

    no_family = genes["family_id"] == ""
    no_segment = ~no_family & (genes["n_segments"] == 0)
    logger.info(f"TF genes: {len(genes)}; with more than one DBD: {int((genes['n_dbds'] > 1).sum())}; "
                f"no CisBP family for their DBD types: {int(no_family.sum())}; "
                f"family but no segment at the segment threshold: {int(no_segment.sum())}")

    top_families = genes["family_name"].replace("", "<no family>").value_counts().head(15)
    logger.info(f"TF genes per family (top 15):\n{top_families.to_string()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--target", required=True, help="target config, a path or a file name in config/")
    parser.add_argument("--log_path", default=None, help="default: logs/call_tfs/<label>.log")
    args = parser.parse_args()

    config = files.reference_config()
    target = files.load_toml(args.target)
    label = target["organism"]["label"]
    configure_logging(args.log_path or LOGS_DIR / "call_tfs" / f"{label}.log")
    logger.info(f"calling TFs of {target['organism']['name']} ({label})")

    ref_dir = files.reference_dir(config)
    out = files.results_dir(label)
    out.mkdir(parents=True, exist_ok=True)

    # loading the reference and the target proteome
    _, families, models, _, _ = reference.load(ref_dir)
    dbd_table = pd.read_csv(ref_dir / DBD_MODELS_TSV, sep="\t")
    calling_types = set(dbd_table.loc[dbd_table["calls_tfs"], "model"])
    proteins, inputs, notes = proteome.load(target)
    sequences = proteins.set_index("protein_id")["sequence"].to_dict()

    # calling DNA-binding domain hits
    hits = hmmer.scan(ref_dir / DBD_MODELS_HMM, sequences)
    hits.to_csv(out / DBD_HITS_TSV, sep="\t", index=False)
    calling_i_evalue = config["rules"]["calling_i_evalue"]
    called = call_domains(hits, models, calling_types, proteins, calling_i_evalue)

    # one representative protein per TF gene, its family and its SR segments
    genes = representatives(called, proteins)
    genes = add_families(genes, families)
    genes = add_segments(genes, hits, sequences, families, ref_dir, models,
                         config["rules"]["segment_max_i_evalue"])

    genes[COLUMNS].to_csv(out / TF_CALLS_TSV, sep="\t", index=False)
    files.write_provenance(out, inputs + [ref_dir / DBD_MODELS_HMM, ref_dir / FAMILIES_TSV],
                           name="provenance_call_tfs", proteome=notes, calling_i_evalue=calling_i_evalue)
    report(genes, notes)


if __name__ == "__main__":
    main()
