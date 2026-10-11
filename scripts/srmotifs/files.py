"""Configs, output folders, write-once downloads into data/raw, and provenance records."""

import datetime
import hashlib
import json
import logging
import subprocess
import tomllib
from pathlib import Path

from srmotifs.constants import (CONFIG_DIR, DATA_RAW, LEDGER, REFERENCE_CONFIG, REFERENCE_ROOT,
                                RESULTS_ROOT, ROOT)

LEDGER_FIELDS = ["path", "source", "retrieved_utc", "bytes", "sha256"]

logger = logging.getLogger(__name__)


class FetchError(RuntimeError):
    pass


def load_toml(path):
    """Read a TOML config; a bare file name is looked up in config/."""
    path = Path(path)
    if not path.is_absolute() and not path.exists():
        path = CONFIG_DIR / path

    with open(path, "rb") as handle:
        return tomllib.load(handle)


def reference_config():
    return load_toml(REFERENCE_CONFIG)


def reference_dir(config):
    name = f"cisbp_{config['cisbp']['build']}_pfam_{config['pfam']['release']}"
    return REFERENCE_ROOT / name


def results_dir(label):
    return RESULTS_ROOT / label


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def run(command):
    """Run an external tool and return its output; raise with its error message if it fails."""
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{command[0]} failed ({result.returncode}): {result.stderr.strip()}")
    return result.stdout


def read_ledger():
    """{path relative to data/raw: ledger row} for every fetched file."""
    if not LEDGER.exists():
        return {}

    rows = {}
    for line in LEDGER.read_text().splitlines()[1:]:
        row = dict(zip(LEDGER_FIELDS, line.split("\t")))
        rows[row["path"]] = row
    return rows


def record(path, source):
    """Add one fetched file to the ledger. Each file may be recorded only once."""
    relative = Path(path).relative_to(DATA_RAW).as_posix()
    if relative in read_ledger():
        raise RuntimeError(f"{relative} is already recorded; data/raw is write-once")

    if not LEDGER.exists():
        LEDGER.write_text("\t".join(LEDGER_FIELDS) + "\n")

    row = [relative, source, utc_now(), str(Path(path).stat().st_size), sha256(path)]
    with open(LEDGER, "a") as handle:
        handle.write("\t".join(row) + "\n")


def is_recorded(path):
    """True if `path` exists and still matches its ledger hash.

    Raises FetchError if the file exists but was never recorded, or has changed since.
    """
    path = Path(path)
    if not path.exists():
        return False

    row = read_ledger().get(path.relative_to(DATA_RAW).as_posix())
    if row is None:
        raise FetchError(f"{path} exists but is not in {LEDGER.name}; refusing to adopt it")
    if sha256(path) != row["sha256"]:
        raise FetchError(f"{path} no longer matches its recorded SHA-256")
    return True


def curl(url, dest, form_file=None):
    """Download `url` to `dest`, POSTing the contents of `form_file` if one is given."""
    command = ["curl", "--fail", "--location", "--silent", "--show-error",
               "--retry", "5", "--retry-delay", "10", "--connect-timeout", "60",
               "--output", str(dest), url]
    if form_file:
        command += ["--data-binary", f"@{form_file}"]

    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise FetchError(f"curl failed ({result.returncode}) for {url}: {result.stderr.strip()}")
    if Path(dest).stat().st_size == 0:
        raise FetchError(f"empty download from {url}")


def fetch(url, dest):
    """Fetch `url` into data/raw once; later calls only check the recorded hash."""
    dest = Path(dest)
    if is_recorded(dest):
        logger.info(f"already in data/raw, hash checked: {dest.relative_to(DATA_RAW)}")
        return dest

    # download under a temporary name, so an interrupted download never looks finished
    logger.info(f"downloading {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name + ".part")
    curl(url, partial)
    partial.rename(dest)

    record(dest, url)
    logger.info(f"recorded {dest.relative_to(DATA_RAW)} in {LEDGER.name}")
    return dest


def tool_versions():
    """Version strings of the external tools, as they report them."""
    hmmer_help = run(["hmmscan", "-h"])
    hmmer_line = next(line for line in hmmer_help.splitlines() if "HMMER" in line)
    return {"hmmer": hmmer_line.lstrip("# ").strip(),
            "gffread": run(["gffread", "--version"]).strip()}


def code_version():
    """Git commit of this repository, with '-dirty' if tracked files have uncommitted changes."""
    def git(*args):
        command = ["git", "-C", str(ROOT), *args]
        return subprocess.run(command, capture_output=True, text=True).stdout.strip()

    commit = git("rev-parse", "--verify", "--quiet", "HEAD") or "uncommitted"
    if git("status", "--porcelain", "--untracked-files=no"):
        commit += "-dirty"
    return commit


def write_provenance(out_dir, inputs, name="provenance", **extra):
    """<name>.json: code version, tool versions, and every input with its SHA-256.

    Steps writing into the same folder each pass their own `name`, so none overwrites
    another's record.
    """
    hashes = {}
    for path in map(Path, inputs):
        shown = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
        hashes[str(shown)] = sha256(path)

    provenance = {"written_utc": utc_now(), "code": code_version(), "tools": tool_versions(),
                  "inputs": hashes, **extra}
    (Path(out_dir) / f"{name}.json").write_text(json.dumps(provenance, indent=2) + "\n")
