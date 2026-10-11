"""Write, subset and validate MEME minimal-format (version 4) motif files."""

import re
from pathlib import Path

import numpy as np

DECIMALS = 6
HEADER = ("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
          "Background letter frequencies\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n")

# CisBP ships probabilities, not site counts; nsites is a fixed pseudocount weight for
# tools that use one, not a number of observed sites
NSITES = 20


def format_rows(matrix):
    """Rows printed with DECIMALS places, adjusted so each printed row sums to exactly 1."""
    scale = 10 ** DECIMALS
    units = np.floor(matrix * scale + 0.5).astype(np.int64)

    # rounding can leave a row a unit off; the largest entry absorbs the difference
    for row in units:
        row[np.argmax(row)] += scale - row.sum()
    if (units < 0).any():
        raise ValueError("rounding produced a negative probability")

    lines = [" ".join(f"{value // scale}.{value % scale:0{DECIMALS}d}" for value in row) for row in units]
    return "".join(line + "\n" for line in lines)


def block(motif_id, alt_name, matrix):
    """One MOTIF entry."""
    alt_name = re.sub(r"[^A-Za-z0-9_.-]", "_", alt_name) or "unnamed"
    return (f"MOTIF {motif_id} {alt_name}\n"
            f"letter-probability matrix: alength= 4 w= {len(matrix)} nsites= {NSITES} E= 0\n"
            f"{format_rows(matrix)}\n")


def write(path, blocks):
    """Write {motif id: block} in sorted id order, then validate the written file."""
    Path(path).write_text(HEADER + "".join(blocks[motif_id] for motif_id in sorted(blocks)))
    validate(path)


def read_blocks(path):
    """{motif id: block text} from a file written by `write`."""
    chunks = Path(path).read_text().split("\nMOTIF ")[1:]
    return {chunk.split()[0]: "MOTIF " + chunk.rstrip("\n") + "\n\n" for chunk in chunks}


def validate(path):
    """Check unique filesystem-safe ids, declared widths, 4 columns and rows summing to 1."""
    blocks = read_blocks(path)
    if len(blocks) != Path(path).read_text().count("\nMOTIF "):
        raise ValueError(f"{path}: duplicate motif ids")

    for motif_id, text in blocks.items():
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", motif_id):
            raise ValueError(f"{path}: unsafe motif id {motif_id!r}")

        lines = text.strip().splitlines()
        width = int(re.search(r"w= (\d+)", lines[1]).group(1))
        rows = np.array([[float(value) for value in line.split()] for line in lines[2:]])
        if rows.shape != (width, 4) or np.abs(rows.sum(axis=1) - 1).max() > 1e-6:
            raise ValueError(f"{path}: malformed matrix for {motif_id}")

    return len(blocks)
