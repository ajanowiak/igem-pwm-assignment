"""Download CisBP and read its TF table, PWMs and protein sequences."""

import io
import logging
import re
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from srmotifs import files
from srmotifs.constants import CISBP_PROTEIN_TABLE, CISBP_PWM_ARCHIVE, CISBP_TF_TABLE, DATA_RAW

# TFs without a recognised DNA-binding domain cannot be reached by DBD homology, so they
# are never references or check targets
UNREACHABLE_FAMILIES = {"Unknown"}

logger = logging.getLogger(__name__)


def raw_dir(config):
    return DATA_RAW / f"cisbp_{config['cisbp']['build']}"


def fetch_tables(config):
    """The TF table and the PWM archive, from CisBP's bulk downloads."""
    url = config["cisbp"]["bulk_url"]
    tf_table = files.fetch(url + CISBP_TF_TABLE, raw_dir(config) / CISBP_TF_TABLE)
    pwm_archive = files.fetch(url + CISBP_PWM_ARCHIVE, raw_dir(config) / CISBP_PWM_ARCHIVE)
    return tf_table, pwm_archive


def fetch_proteins(config, name, tf_ids):
    """Fetch CisBP's protein table for `tf_ids` in fixed-size batches; return the ids it has no protein for.

    The cart archive form names its zip by the minute only, so a batch is accepted only
    if every TF it lists was requested in that batch.
    """
    # splitting the sorted ids into batches, so reruns request exactly the same batches
    ids = sorted(tf_ids)
    size = config["cisbp"]["protein_batch_size"]
    batches = [ids[i:i + size] for i in range(0, len(ids), size)]

    missing = set()
    for number, batch in enumerate(batches, start=1):
        dest = raw_dir(config) / "proteins" / f"{name}_{number:02d}_of_{len(batches):02d}.zip"
        if not files.is_recorded(dest):
            logger.info(f"requesting CisBP proteins, {name} batch {number} of {len(batches)}")
            download_protein_batch(config["cisbp"]["site_url"], batch, dest)

        listed = set(read_protein_table(dest)["TF_ID"])
        if not listed or listed - set(batch):
            raise files.FetchError(f"{dest.name} lists TFs that were not requested in this batch")
        missing |= set(batch) - listed

    return missing


def download_protein_batch(site_url, tf_ids, dest):
    """POST one batch of TF ids to CisBP's cart archive form and keep the returned zip."""
    form = "studyID=0&dataSourceTF=0&prot_seq=on&submit=Download"
    form += "".join(f"&tfID[]={tf_id}" for tf_id in tf_ids)

    with tempfile.TemporaryDirectory() as tmp:
        form_file = Path(tmp) / "form.txt"
        page = Path(tmp) / "page.html"
        archive = Path(tmp) / "archive.zip"
        form_file.write_text(form)

        # the form answers with a page that links to a freshly built zip
        files.curl(site_url + "archive.php", page, form_file)
        link = re.search(r'href="(tmp/[^"]+\.zip)"', page.read_text(errors="replace"))
        if not link:
            raise files.FetchError("CisBP's archive form returned no archive link")
        files.curl(site_url + link.group(1), archive)

        # someone else's download in the same minute would overwrite ours
        listed = set(read_protein_table(archive)["TF_ID"])
        if not listed or listed - set(tf_ids):
            raise files.FetchError(f"archive {link.group(1)} is not the one requested; retry in a minute")

        dest.parent.mkdir(parents=True, exist_ok=True)
        archive.rename(dest)

    source = (f"POST {site_url}archive.php prot_seq for {len(tf_ids)} TF ids "
              f"{tf_ids[0]}..{tf_ids[-1]} (archive {link.group(1)})")
    files.record(dest, source)


def load_tf_table(path):
    with zipfile.ZipFile(path) as archive:
        (name,) = [n for n in archive.namelist() if n.endswith(".txt")]
        with archive.open(name) as handle:
            return pd.read_csv(handle, sep="\t", dtype=str, keep_default_na=False, low_memory=False)


def families(tf_table):
    """One row per family: DBD types, SR model type, and CisBP's similarity cutoff."""
    columns = ["Family_ID", "Family_Name", "DBDs", "SR_Model", "SR_NoThreshold", "Cutoff"]
    table = tf_table[columns].drop_duplicates()
    if table["Family_ID"].duplicated().any():
        raise ValueError("a family has inconsistent rows in the TF table")

    table["Cutoff"] = table["Cutoff"].astype(float)
    return table.set_index("Family_ID")


def reference_tfs(tf_table):
    """One row per TF with at least one directly measured motif (TF_Status D)."""
    measured = tf_table["TF_Status"] == "D"
    reachable = ~tf_table["Family_Name"].isin(UNREACHABLE_FAMILIES)
    columns = ["TF_ID", "TF_Name", "TF_Species", "Family_ID", "Family_Name"]

    references = tf_table.loc[measured & reachable, columns].drop_duplicates("TF_ID")
    return references.reset_index(drop=True)


def check_targets(tf_table, species):
    """CisBP's I and N TFs of `species`, with the reference TF each I motif came from."""
    selected = (tf_table["TF_Species"].isin(species)
                & tf_table["TF_Status"].isin(["I", "N"])
                & ~tf_table["Family_Name"].isin(UNREACHABLE_FAMILIES))
    targets = tf_table[selected].drop_duplicates("TF_ID").copy()

    # an inferred motif belongs to exactly one measured TF: that is the source
    owner = tf_table[tf_table["TF_Status"] == "D"].set_index("Motif_ID")["TF_ID"]
    targets["cisbp_source_tf"] = targets["Motif_ID"].map(owner).where(targets["TF_Status"] == "I", "")

    columns = ["TF_ID", "TF_Name", "TF_Species", "Family_ID", "Family_Name",
               "TF_Status", "Motif_ID", "cisbp_source_tf"]
    return targets[columns].reset_index(drop=True)


def calling_dbd_types(tf_table, min_species):
    """({DBD type: number of species CisBP calls TFs with it}, DBD types used in >= min_species)."""
    species = {}
    for dbds, group in tf_table.drop_duplicates("TF_ID").groupby("DBDs"):
        for dbd in dbds.split(","):
            species.setdefault(dbd, set()).update(group["TF_Species"])

    counts = {dbd: len(seen) for dbd, seen in species.items()}
    calling = {dbd for dbd, count in counts.items() if count >= min_species}
    return counts, calling


def read_protein_table(zip_path):
    with zipfile.ZipFile(zip_path) as archive:
        with archive.open(CISBP_PROTEIN_TABLE) as handle:
            return pd.read_csv(handle, sep="\t", dtype=str, keep_default_na=False)


def load_proteins(config):
    """Every fetched protein row, one per (TF_ID, Protein_ID).

    Surrounding whitespace and a terminal '*' (stop) are dropped; any other
    non-residue character is an error.
    """
    paths = sorted((raw_dir(config) / "proteins").glob("*.zip"))
    if not paths:
        raise FileNotFoundError("no CisBP protein tables fetched")
    table = pd.concat([read_protein_table(p) for p in paths], ignore_index=True)

    # cleaning the sequences
    table["Protein_seq"] = table["Protein_seq"].str.strip().str.removesuffix("*")
    invalid = table[~table["Protein_seq"].str.fullmatch(r"[A-Z]+")]
    if len(invalid):
        raise ValueError(f"unexpected characters in {len(invalid)} CisBP proteins, "
                         f"e.g. {invalid['TF_ID'].iloc[0]}")

    # the same protein can be listed in two batches; it must be the same sequence
    if (table.groupby(["TF_ID", "Protein_ID"])["Protein_seq"].nunique() > 1).any():
        raise ValueError("CisBP lists one protein id with two different sequences")
    return table.drop_duplicates(["TF_ID", "Protein_ID"]).reset_index(drop=True)


def load_matrices(pwm_zip):
    """{Motif_ID: matrix (width x 4, columns A C G T)}, rows renormalized to sum to 1.

    CisBP lists some motifs with a header but no matrix rows; those are skipped.
    """
    matrices = {}
    with zipfile.ZipFile(pwm_zip) as archive:
        for name in archive.namelist():
            if not name.endswith(".txt"):
                continue

            table = pd.read_csv(io.BytesIO(archive.read(name)), sep="\t")
            if table.empty:
                continue

            matrix = table[["A", "C", "G", "T"]].to_numpy(dtype=float)
            if (matrix < 0).any() or np.any(matrix.sum(axis=1) <= 0):
                raise ValueError(f"{name}: invalid probabilities")
            matrices[Path(name).stem] = matrix / matrix.sum(axis=1, keepdims=True)

    return matrices


def representative_motifs(tf_table, motifs_with_matrix):
    """{reference TF_ID: Motif_ID} — the one PWM that gets copied from each reference TF.

    A reference TF can have several lab-measured PWMs (different studies, different
    experiments). CisBP's TF table lists every PWM it copied to an unmeasured TF (rows
    with TF_Status I), and each PWM id belongs to exactly one reference TF. So for each
    reference TF we count how many times CisBP copied each of its PWMs, and take the
    most-copied one. If CisBP never copied any of them, we take the first one it lists.
    Only PWMs with matrix rows count; a reference TF with none gives nothing.
    """
    has_matrix = tf_table["Motif_ID"].isin(motifs_with_matrix)
    measured = tf_table[(tf_table["TF_Status"] == "D") & has_matrix]
    copied = tf_table[(tf_table["TF_Status"] == "I") & has_matrix]

    # counting how often CisBP copied each PWM of each reference TF
    owner = measured.set_index("Motif_ID")["TF_ID"]
    copied = copied.assign(source=copied["Motif_ID"].map(owner)).dropna(subset=["source"])
    counts = copied.groupby(["source", "Motif_ID"]).size().reset_index(name="n")
    counts = counts.sort_values(["source", "n", "Motif_ID"], ascending=[True, False, True])
    chosen = counts.drop_duplicates("source").set_index("source")["Motif_ID"].to_dict()

    # reference TFs CisBP never copied from keep their first-listed PWM
    for tf_id, motif_id in measured.groupby("TF_ID")["Motif_ID"].first().items():
        chosen.setdefault(tf_id, motif_id)
    return chosen
