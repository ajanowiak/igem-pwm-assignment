"""Constants shared across the pipeline: project layout, file names, tier names.

Values that only one module uses stay in that module. Settings a user may change
(builds, URLs, thresholds) live in config/*.toml, not here.
"""

from pathlib import Path

# project layout; this file sits at scripts/srmotifs/constants.py
ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "config"
REFERENCE_CONFIG = CONFIG_DIR / "reference.toml"
DATA_RAW = ROOT / "data" / "raw"
LEDGER = DATA_RAW / "SOURCES.tsv"
REFERENCE_ROOT = ROOT / "reference"
RESULTS_ROOT = ROOT / "results"
CHECKS_DIR = RESULTS_ROOT / "checks" / "sr_vs_cisbp"
LOGS_DIR = ROOT / "logs"

# files of the built reference, written by build_reference.py and read by the other scripts
REFERENCES_TSV = "references.tsv"
FAMILIES_TSV = "families.tsv"
DBD_MODELS_HMM = "dbd_models.hmm"
DBD_MODELS_TSV = "dbd_models.tsv"
MOTIFS_MEME = "motifs.meme"
SR_DIR = "sr"

# files of one target organism, written by call_tfs.py and assign_motifs.py
TF_CALLS_TSV = "tf_calls.tsv"
DBD_HITS_TSV = "dbd_hits.tsv"
ASSIGNMENTS_TSV = "assignments.tsv"

# CisBP file names: the bulk downloads, and the table inside each cart archive
CISBP_TF_TABLE = "TF_Information.txt.zip"
CISBP_PWM_ARCHIVE = "PWMs.zip"
CISBP_PROTEIN_TABLE = "prot_seq.txt"

# where each PWM comes from: references.tsv column -> CisBP TF table column
MOTIF_SOURCE_COLUMNS = {"motif_experiment": "MSource_Type", "motif_study": "MSource_Author",
                        "motif_year": "MSource_Year", "motif_pmid": "PMID"}

# SR's class for a pair of highly similar DBDs, and the model of families without a cutoff
SR_HSIM = "HSim"
SR_NO_THRESHOLD_MODEL = "NO_THRESHOLD"

# assignment tiers of a target TF gene
TIER_PASS = "pass"
TIER_BELOW = "below-threshold"
TIER_NONE = "none"
