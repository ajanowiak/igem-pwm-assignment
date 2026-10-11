"""Get a target organism's proteins with their gene ids, from genome + GFF3 or from a protein FASTA."""

import gzip
import shutil
import tempfile
from pathlib import Path

import pandas as pd

from srmotifs import files

# what gffread reports as the gene id of a transcript that has no parent gene
NO_GENE = {"", "."}


def input_path(value, label):
    """A URL is fetched once into data/raw/targets/<label>/; a local path is used as given."""
    if value.startswith(("http://", "https://", "ftp://")):
        return files.fetch(value, files.DATA_RAW / "targets" / label / Path(value).name)

    path = Path(value).expanduser()
    if not path.is_absolute():
        path = files.ROOT / path
    if not path.is_file():
        raise FileNotFoundError(f"input not found: {path}")
    return path


def read_fasta(path):
    """{record id: sequence}; the id is the header up to the first space."""
    opener = gzip.open if str(path).endswith(".gz") else open

    records, name, chunks = {}, None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())

    if name is not None:
        records[name] = "".join(chunks)
    return records


def plain_copy(path, dest):
    """Copy a possibly gzipped file uncompressed (gffread needs an indexable genome)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as source, open(dest, "wb") as target:
        shutil.copyfileobj(source, target)
    return dest


def from_genome(genome, annotation):
    """Translate every coding transcript of the annotation with gffread.

    The gene id is the transcript's parent gene. A transcript without one (some
    annotations have no gene features) becomes its own gene, and is counted.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        genome_fa = plain_copy(genome, tmp / "genome.fa")
        gff = plain_copy(annotation, tmp / "annotation.gff")

        # translating coding transcripts, and listing each transcript's parent gene
        files.run(["gffread", "-g", str(genome_fa), "-y", str(tmp / "proteins.faa"), str(gff)])
        files.run(["gffread", "--table", "@id,@geneid", "-o", str(tmp / "genes.tsv"), str(gff)])

        sequences = read_fasta(tmp / "proteins.faa")
        genes = pd.read_csv(tmp / "genes.tsv", sep="\t", header=None, names=["protein_id", "gene_id"],
                            dtype=str, keep_default_na=False)

    # keeping coding transcripts only (gffread lists non-coding ones too)
    table = genes[genes["protein_id"].isin(sequences)].drop_duplicates("protein_id").copy()
    table["sequence"] = table["protein_id"].map(sequences)

    # transcripts without a parent gene stand for themselves
    orphan = table["gene_id"].isin(NO_GENE)
    table.loc[orphan, "gene_id"] = table.loc[orphan, "protein_id"]
    return table, int(orphan.sum())


def from_protein_fasta(path):
    """Each FASTA record is treated as one gene."""
    sequences = read_fasta(path)
    table = pd.DataFrame({"protein_id": list(sequences), "gene_id": list(sequences),
                          "sequence": list(sequences.values())})
    return table, len(table)


def load(target):
    """(proteins table with protein_id, gene_id, sequence; input files used; notes for provenance).

    A terminal stop ('.' or '*') is removed. Proteins with an internal stop are dropped
    and counted in the notes.
    """
    label, source = target["organism"]["label"], target["input"]

    if "genome" in source and "annotation" in source:
        inputs = [input_path(source["genome"], label), input_path(source["annotation"], label)]
        table, without_gene = from_genome(*inputs)
    elif "proteins" in source:
        inputs = [input_path(source["proteins"], label)]
        table, without_gene = from_protein_fasta(inputs[0])
    else:
        raise ValueError("[input] needs genome + annotation, or proteins")

    # removing terminal stops; an internal stop means the protein is not usable
    table["sequence"] = table["sequence"].str.upper().str.rstrip(".*")
    internal_stop = table["sequence"].str.contains(r"[.*]")
    table = table[~internal_stop].reset_index(drop=True)

    notes = {"proteins": len(table) + int(internal_stop.sum()),
             "dropped_internal_stop": int(internal_stop.sum()),
             "proteins_without_gene_id": without_gene,
             "genes": int(table["gene_id"].nunique())}
    return table, inputs, notes
