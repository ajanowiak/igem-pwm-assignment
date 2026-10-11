"""DBD detection with hmmscan, and SR's per-domain match-state segments with hmmalign.

The same functions encode reference TFs and target TFs, so both sides of every SR
comparison are built identically.
"""

import logging
import os
import tempfile
from pathlib import Path

import pandas as pd

from srmotifs.files import run

DOMTBL_COLUMNS = ["model", "model_acc", "model_len", "protein", "protein_acc", "protein_len",
                  "seq_evalue", "seq_score", "seq_bias", "dom_index", "dom_count", "c_evalue",
                  "i_evalue", "dom_score", "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to",
                  "env_from", "env_to", "acc"]
HIT_COLUMNS = ["protein", "model", "seq_score", "dom_score", "i_evalue", "ali_from", "ali_to"]

# one permissive hmmscan pass reports every hit; the calling rule (Pfam GA) and the
# segment rule (i-Evalue) are filters applied to it afterwards
REPORTING_EVALUE = "1"

logger = logging.getLogger(__name__)


def write_fasta(sequences, path):
    Path(path).write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))


def build_library(pfam_hmm_gz, names, dest):
    """Extract the named models from Pfam into `dest` and press it for hmmscan."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile("w", suffix=".txt") as name_file:
        name_file.write("".join(f"{name}\n" for name in sorted(names)))
        name_file.flush()
        dest.write_text(run(["hmmfetch", "-f", str(pfam_hmm_gz), name_file.name]))

    run(["hmmpress", "-f", str(dest)])
    return read_library(dest)


def read_library(path):
    """{model name: (length, sequence GA, domain GA)} from an HMM file."""
    models, name, length = {}, None, None
    for line in Path(path).read_text().splitlines():
        if line.startswith("NAME"):
            name = line.split()[1]
        elif line.startswith("LENG"):
            length = int(line.split()[1])
        elif line.startswith("GA"):
            seq_ga, dom_ga = (float(value.rstrip(";")) for value in line.split()[1:3])
            models[name] = (length, seq_ga, dom_ga)
    return models


def scan(library, sequences):
    """hmmscan {protein id: sequence} against the library; one row per reported domain hit."""
    logger.info(f"hmmscan: {len(sequences)} proteins against {Path(library).name}")
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "proteins.fa"
        domtbl = Path(tmp) / "hits.domtbl"
        write_fasta(sequences, fasta)

        run(["hmmscan", "-E", REPORTING_EVALUE, "--domE", REPORTING_EVALUE, "--noali",
             "--cpu", str(os.cpu_count() or 1), "--domtblout", str(domtbl), "-o", os.devnull,
             str(library), str(fasta)])

        # the last column (description) may contain spaces, so split a fixed number of times
        lines = [line for line in domtbl.read_text().splitlines() if not line.startswith("#")]
        rows = [line.split(maxsplit=len(DOMTBL_COLUMNS))[:len(DOMTBL_COLUMNS)] for line in lines]

    hits = pd.DataFrame(rows, columns=DOMTBL_COLUMNS)[HIT_COLUMNS]
    hits = hits.astype({"seq_score": float, "dom_score": float, "i_evalue": float,
                        "ali_from": int, "ali_to": int})
    return hits.sort_values(["protein", "ali_from", "model"]).reset_index(drop=True)


def gathering_pass(hits, models):
    """Hits passing Pfam GA the way hmmscan --cut_ga applies it: both scores at or above GA."""
    seq_ga = hits["model"].map(lambda model: models[model][1])
    dom_ga = hits["model"].map(lambda model: models[model][2])
    return hits[(hits["seq_score"] >= seq_ga) & (hits["dom_score"] >= dom_ga)]


def resolve_overlaps(hits):
    """Within each protein, keep the best-scoring hit wherever hits overlap."""
    kept = []
    for _, group in hits.groupby("protein"):
        accepted = []
        for hit in group.sort_values(["dom_score", "model"], ascending=[False, True]).to_dict("records"):
            overlaps = any(hit["ali_from"] <= other["ali_to"] and other["ali_from"] <= hit["ali_to"]
                           for other in accepted)
            if not overlaps:
                accepted.append(hit)
        kept.extend(accepted)

    kept = pd.DataFrame(kept, columns=hits.columns)
    return kept.sort_values(["protein", "ali_from"]).reset_index(drop=True)


def segment_hits(hits, family_dbds, max_i_evalue):
    """Hits of each protein's own family DBD types at or below `max_i_evalue`, overlaps resolved.

    `family_dbds` maps protein -> its family's DBD types.
    """
    own_type = [model in family_dbds.get(protein, ())
                for protein, model in zip(hits["protein"], hits["model"])]
    own_type = pd.Series(own_type, index=hits.index, dtype=bool)

    selected = hits[own_type & (hits["i_evalue"] <= max_i_evalue)]
    return resolve_overlaps(selected)


def segments(sequences, domains, family_dbds, library, models):
    """{protein: [segment]} — SR's input: one match-state string per DBD, in protein order.

    Each segment spans every DBD type of the family, in CisBP's order (`family_dbds`
    maps protein -> ordered DBD types): the domain's own block holds its match states
    and the other types' blocks are gaps, as in SR's training code.
    """
    # aligning every domain to its own HMM, one hmmalign run per model
    match_states = {}
    with tempfile.TemporaryDirectory() as tmp:
        for model, group in domains.groupby("model"):
            model_path = Path(tmp) / f"{model}.hmm"
            model_path.write_text(run(["hmmfetch", str(library), model]))

            pieces = {f"{hit.protein}#{hit.ali_from}": sequences[hit.protein][hit.ali_from - 1:hit.ali_to]
                      for hit in group.itertuples()}
            match_states.update(align_to_model(model_path, pieces, models[model][0]))

    # joining each domain's match states with gap blocks for the family's other DBD types
    result = {}
    for protein, group in domains.groupby("protein"):
        protein_segments = []
        for hit in group.sort_values("ali_from").itertuples():
            blocks = [match_states[f"{protein}#{hit.ali_from}"] if dbd == hit.model else "-" * models[dbd][0]
                      for dbd in family_dbds[protein]]
            protein_segments.append("".join(blocks))
        result[protein] = protein_segments
    return result


def align_to_model(hmm_path, sequences, model_length):
    """hmmalign {name: sequence} to one HMM; return {name: match-state string}."""
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "pieces.fa"
        write_fasta(sequences, fasta)
        stockholm = run(["hmmalign", "--outformat", "Stockholm", str(hmm_path), str(fasta)])

    aligned, reference = read_stockholm(stockholm)
    if reference.count("x") != model_length:
        raise ValueError(f"RF line has {reference.count('x')} match states, model has {model_length}")
    if set(sequences) - set(aligned):
        raise ValueError("hmmalign dropped sequences")

    return {name: match_state_string(aligned[name], reference) for name in sequences}


def read_stockholm(text):
    """({sequence name: aligned string}, RF annotation) from one Stockholm block."""
    aligned, reference = {}, []
    for line in text.splitlines():
        if not line.strip() or line.startswith("# STOCKHOLM") or line.startswith("//"):
            continue
        if line.startswith("#=GC RF"):
            reference.append(line.split()[2])
        elif not line.startswith("#"):
            name, chunk = line.split()
            aligned[name] = aligned.get(name, "") + chunk
    return aligned, "".join(reference)


def match_state_string(aligned_row, reference):
    """Residue or '-' at each HMM match state, in the semiglobal convention SR was trained on.

    SR's alignments (R package aphid) place every residue of a DBD segment: residues
    hmmalign leaves unaligned before the first (after the last) occupied match state
    are moved into the empty match states immediately preceding (following) it.
    Flanking residues beyond the available empty states stay unaligned.
    """
    if len(aligned_row) != len(reference):
        raise ValueError("aligned row and RF line differ in length")

    # reading the residue (or gap) at each match state; 'x' in RF marks a match column
    columns = list(zip(aligned_row, reference))
    states = [char.upper() if char not in ".-" else "-" for char, rf in columns if rf == "x"]
    occupied_columns = [i for i, (char, rf) in enumerate(columns) if rf == "x" and char not in ".-"]
    if not occupied_columns:
        return "".join(states)

    # residues hanging off either end of the aligned core
    n_flank = [char.upper() for char, _ in columns[:occupied_columns[0]] if char not in ".-"]
    c_flank = [char.upper() for char, _ in columns[occupied_columns[-1] + 1:] if char not in ".-"]

    occupied_states = [i for i, state in enumerate(states) if state != "-"]
    first_state, last_state = occupied_states[0], occupied_states[-1]

    # moving the N-terminal flank into the empty states just before the core
    n_fill = n_flank[len(n_flank) - min(len(n_flank), first_state):]
    for offset, residue in enumerate(n_fill):
        states[first_state - len(n_fill) + offset] = residue

    # moving the C-terminal flank into the empty states just after the core
    c_fill = c_flank[:len(states) - 1 - last_state]
    for offset, residue in enumerate(c_fill):
        states[last_state + 1 + offset] = residue

    return "".join(states)
