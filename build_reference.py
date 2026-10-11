"""Step A: download CisBP, Pfam and SR, and build the reference SR searches against.

Usage: python build_reference.py

Downloads (write-once, SHA-256 recorded in data/raw/SOURCES.tsv): CisBP's TF table,
PWMs and protein sequences; Pfam-A HMMs of the pinned release; SR at a pinned commit.

Writes reference/cisbp_<build>_pfam_<release>/:
  references.tsv   one row per reference TF (a CisBP TF with at least one lab-measured
                   PWM): its family, the PWM that gets copied from it and that PWM's
                   source experiment, the CisBP isoform used, and its SR DBD segments
                   (blank segments = unusable, with the reason in `outcome`)
  motifs.meme      the PWMs named in references.tsv
  families.tsv     per family: DBD types, SR model, CisBP's similarity cutoff
  dbd_models.hmm   the Pfam models of CisBP's DBD types (pressed), and dbd_models.tsv
                   saying which of them call TFs
  sr/              the extracted SR code and family models
  provenance.json

Isoform choice for each reference: see srmotifs/reference.py.
"""

import gzip
import json

import pandas as pd

from srmotifs import cisbp, files, hmmer, meme, reference, sr

# where each PWM comes from (experiment, study), copied from CisBP's TF table into references.tsv
MOTIF_COLUMNS = {"motif_experiment": "MSource_Type", "motif_study": "MSource_Author",
                 "motif_year": "MSource_Year", "motif_pmid": "PMID"}


def fetch_all(config):
    """Download everything the reference needs; return the paths and the TF table."""
    tf_zip, pwm_zip = cisbp.fetch_tables(config)

    pfam_dest = files.DATA_RAW / f"pfam_{config['pfam']['release']}" / "Pfam-A.hmm.gz"
    pfam = files.fetch(config["pfam"]["hmm_url"], pfam_dest)

    commit = config["similarity_regression"]["commit"]
    sr_url = config["similarity_regression"]["tarball_url"].format(commit=commit)
    sr_tarball = files.fetch(sr_url, files.DATA_RAW / f"SimilarityRegression_{commit[:7]}.tar.gz")

    # protein sequences for the references, and for the TFs the SR check scores
    tf_table = cisbp.load_tf_table(tf_zip)
    reference_ids = cisbp.reference_tfs(tf_table)["TF_ID"]
    check_ids = cisbp.check_targets(tf_table, config["check"]["species"])["TF_ID"]
    missing = {"references": cisbp.fetch_proteins(config, "references", reference_ids),
               "check": cisbp.fetch_proteins(config, "check", check_ids)}

    return tf_table, pwm_zip, pfam, sr_tarball, missing


def pfam_names(pfam_hmm_gz):
    with gzip.open(pfam_hmm_gz, "rt") as handle:
        return {line.split()[1] for line in handle if line.startswith("NAME")}


def dbd_model_table(config, tf_table, families, pfam, out):
    """Build the DBD HMM library; return (models, table of DBD types saying which call TFs)."""
    names = {dbd for dbds in families["DBDs"] for dbd in dbds.split(",")}
    models = hmmer.build_library(pfam, names & pfam_names(pfam), out / "dbd_models.hmm")
    species, calling = cisbp.calling_dbd_types(tf_table, config["rules"]["min_cisbp_species"])

    rows = []
    for name in sorted(names):
        rows.append({"model": name, "in_pfam": name in models,
                     "length": models[name][0] if name in models else None,
                     "cisbp_species": species.get(name, 0),
                     "calls_tfs": name in models and name in calling})
    return models, pd.DataFrame(rows).astype({"length": "Int64"})


def family_table(families, models, sr_root):
    """families.tsv: each family's DBD types, SR model file and cutoff.

    Follows CisBP's SR_Model / SR_NoThreshold columns, and asserts that each positional
    SR model's weight vector spans exactly the match states of the family's DBD models.
    """
    rows = []
    for family_id, row in families.iterrows():
        # which SR model file CisBP uses for this family
        if row["SR_NoThreshold"] == "True":
            model_name = "NO_THRESHOLD"
        else:
            model_name = f"{family_id.split('_')[0]}_1.97d"
        path = sr_root / "SRModels" / f"{model_name}.json"
        if not path.exists():
            raise ValueError(f"{family_id} ({row['Family_Name']}) needs SR model {path.name}, "
                             "which SR does not ship")

        # the model's weights must line up with the Pfam match states, position by position
        weights = json.loads(path.read_text()).get("SR.Weights")
        dbds = row["DBDs"].split(",")
        span = sum(models[d][0] for d in dbds) if all(d in models for d in dbds) else None
        if isinstance(weights, list) and len(weights) != span:
            raise ValueError(f"{family_id}: SR weights {len(weights)} != DBD match states {span}")

        rows.append({"Family_ID": family_id, "Family_Name": row["Family_Name"], "DBDs": row["DBDs"],
                     "sr_model": model_name, "cutoff": row["Cutoff"]})

    table = pd.DataFrame(rows).set_index("Family_ID")
    sr.family_models(sr_root, table)
    return table


def add_motifs(references, tf_table, motif_of):
    """Add the PWM copied from each reference TF, and that PWM's source experiment."""
    references.insert(5, "motif_id", references["TF_ID"].map(motif_of).fillna(""))

    measured = tf_table[tf_table["TF_Status"] == "D"].drop_duplicates("Motif_ID").set_index("Motif_ID")
    for column, source in MOTIF_COLUMNS.items():
        references[column] = references["motif_id"].map(measured[source]).fillna("")
    return references


def report(references, dbd_table, missing, config, out):
    usable = references[references["segments"] != ""]
    print(f"reference TFs: {len(references)}; with SR segments: {len(usable)}; "
          f"with segments and a motif: {(usable['motif_id'] != '').sum()}")

    reasons = references.loc[references["segments"] == "", "outcome"]
    reasons = reasons.str.replace(r" \(.*", "", regex=True).str.replace(r"types .*", "types of its family", regex=True)
    print("without segments:", reasons.value_counts().to_dict())

    not_calling = dbd_table.loc[dbd_table["in_pfam"] & ~dbd_table["calls_tfs"], "model"].tolist()
    not_in_pfam = dbd_table.loc[~dbd_table["in_pfam"], "model"].tolist()
    print(f"motifs: {meme.validate(out / 'motifs.meme')}; DBD models: {int(dbd_table['in_pfam'].sum())}, "
          f"calling TFs: {int(dbd_table['calls_tfs'].sum())}; not calling: {not_calling}; "
          f"CisBP DBD types not in Pfam {config['pfam']['release']}: {not_in_pfam}")
    print("CisBP lists no protein for:", {name: len(ids) for name, ids in missing.items()})


def main():
    config = files.reference_config()
    out = files.reference_dir(config)
    out.mkdir(parents=True, exist_ok=True)

    # downloading CisBP, Pfam and SR
    tf_table, pwm_zip, pfam, sr_tarball, missing = fetch_all(config)
    sr_root = sr.extract(sr_tarball, out / "sr")

    # DBD models, and which SR model each family uses
    families = cisbp.families(tf_table)
    models, dbd_table = dbd_model_table(config, tf_table, families, pfam, out)
    dbd_table.to_csv(out / "dbd_models.tsv", sep="\t", index=False)
    family_table(families, models, sr_root).to_csv(out / "families.tsv", sep="\t")

    # one isoform and its SR segments per reference TF
    references, sequences = reference.resolve_tfs(
        cisbp.reference_tfs(tf_table), cisbp.load_proteins(config), families,
        out / "dbd_models.hmm", models, config["rules"]["segment_max_i_evalue"])

    # the one PWM copied from each reference TF
    matrices = cisbp.load_matrices(pwm_zip)
    references = add_motifs(references, tf_table, cisbp.representative_motifs(tf_table, set(matrices)))

    # writing the reference
    references.to_csv(out / "references.tsv", sep="\t", index=False)
    (out / "proteins.fasta").write_text("".join(f">{tf}\n{seq}\n" for tf, seq in sorted(sequences.items())))

    names = references.drop_duplicates("motif_id").set_index("motif_id")["TF_Name"]
    blocks = {motif_id: meme.block(motif_id, names[motif_id], matrices[motif_id])
              for motif_id in names.index if motif_id}
    meme.write(out / "motifs.meme", blocks)

    # the reference depends on every downloaded file except target organisms' inputs
    raw_files = [files.DATA_RAW / path for path in files.read_ledger() if not path.startswith("targets/")]
    files.write_provenance(out, raw_files,
                           cisbp_tfs_without_protein={name: sorted(ids) for name, ids in missing.items()})

    report(references, dbd_table, missing, config, out)


if __name__ == "__main__":
    main()
