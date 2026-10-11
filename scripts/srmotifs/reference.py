"""Resolve CisBP TFs to one isoform with SR segments, and load a built reference."""

import pandas as pd

from srmotifs import hmmer
from srmotifs.constants import DBD_MODELS_HMM, FAMILIES_TSV, REFERENCES_TSV


def choose_isoform(keys, expected, found, counts, sequences):
    """(chosen key or None, outcome text) for one TF's candidate isoforms.

    A candidate is valid if its DBD hits reproduce exactly the family's DBD types.
    Among valid candidates: most segments, then longest, then smallest id.
    """
    if not keys:
        return None, "CisBP lists no protein"

    valid = [key for key in keys if found.get(key) == expected]
    if not valid:
        seen = sorted({",".join(sorted(found[key])) for key in keys if key in found})
        return None, f"no isoform has DBD types {','.join(sorted(expected))} (found {seen or 'none'})"

    valid.sort(key=lambda key: (-counts[key], -len(sequences[key]), key))
    return valid[0], f"resolved ({len(valid)} valid of {len(keys)} isoforms)"


def resolve_tfs(tfs, proteins, families, library, models, max_i_evalue):
    """One row per TF: the chosen CisBP isoform and its SR segments, or why there are none.

    CisBP lists every isoform of a gene without marking the one it scored, so each
    isoform is scanned and one is chosen by `choose_isoform`.
    """
    # every listed isoform, keyed '<TF_ID>|<Protein_ID>'
    candidates = proteins[proteins["TF_ID"].isin(tfs["TF_ID"])]
    sequences = {f"{row.TF_ID}|{row.Protein_ID}": row.Protein_seq for row in candidates.itertuples()}

    keys_by_tf = {}
    for key in sequences:
        keys_by_tf.setdefault(key.split("|")[0], []).append(key)

    # each isoform's family DBD types, in CisBP's order
    family_of = tfs.set_index("TF_ID")["Family_ID"]
    family_dbds = {key: families.loc[family_of[key.split("|")[0]], "DBDs"].split(",") for key in sequences}

    # scanning all isoforms and keeping the hits that count as SR segments
    hits = hmmer.scan(library, sequences)
    domains = hmmer.segment_hits(hits, family_dbds, max_i_evalue)
    found = domains.groupby("protein")["model"].agg(frozenset)
    counts = domains.groupby("protein").size()

    # picking one isoform per TF
    chosen, outcome = {}, {}
    for tf_id in tfs["TF_ID"]:
        keys = keys_by_tf.get(tf_id, [])
        expected = frozenset(family_dbds[keys[0]]) if keys else None
        key, outcome[tf_id] = choose_isoform(keys, expected, found, counts, sequences)
        if key:
            chosen[tf_id] = key

    # aligning the chosen isoforms' domains into SR segments
    chosen_domains = domains[domains["protein"].isin(chosen.values())]
    segment_map = hmmer.segments(sequences, chosen_domains, family_dbds, library, models)

    table = tfs.copy()
    table["protein_id"] = [chosen[t].split("|", 1)[1] if t in chosen else "" for t in table["TF_ID"]]
    table["n_isoforms"] = [len(keys_by_tf.get(t, [])) for t in table["TF_ID"]]
    table["n_segments"] = [len(segment_map[chosen[t]]) if t in chosen else 0 for t in table["TF_ID"]]
    table["segments"] = [",".join(segment_map[chosen[t]]) if t in chosen else "" for t in table["TF_ID"]]
    table["outcome"] = table["TF_ID"].map(outcome)

    chosen_sequences = {tf_id: sequences[key] for tf_id, key in chosen.items()}
    return table, chosen_sequences


def load(ref_dir):
    """Load a built reference.

    Returns (references table, families table, {model: (length, seq GA, dom GA)},
    {Family_ID: [(TF_ID, segments)]}, {TF_ID: Motif_ID}).
    """
    references = pd.read_csv(ref_dir / REFERENCES_TSV, sep="\t", dtype=str, keep_default_na=False)
    families = pd.read_csv(ref_dir / FAMILIES_TSV, sep="\t", dtype={"cutoff": float}).set_index("Family_ID")
    models = hmmer.read_library(ref_dir / DBD_MODELS_HMM)

    # only references with segments can be scored
    usable = references[references["segments"] != ""]
    by_family = {}
    for row in usable.itertuples():
        by_family.setdefault(row.Family_ID, []).append((row.TF_ID, row.segments.split(",")))

    motif_of = {row.TF_ID: row.motif_id for row in usable.itertuples() if row.motif_id}
    return references, families, models, by_family, motif_of
