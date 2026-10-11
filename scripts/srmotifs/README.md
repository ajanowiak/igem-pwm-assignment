# srmotifs

Python package shared by the [pipeline scripts](../README.md). It contains no
command-line entry points. The scripts in `scripts/` import it directly, since Python
places a script's own directory on the import path; no installation is required. Other
code, such as notebooks, can import it after adding `scripts/` to `sys.path`.

| module | purpose |
|---|---|
| `constants.py` | project paths, file names, CisBP file names, SR class and tier names |
| `log.py` | `configure_logging()`: logging to standard output and to an appended log file |
| `files.py` | configuration files, output directories, write-once downloads, provenance records, execution of external tools |
| `cisbp.py` | retrieval of CisBP data and parsing of its TF table, motifs and protein sequences |
| `proteome.py` | proteome of a target organism with gene identifiers, from genome and GFF3 or from a protein FASTA file |
| `hmmer.py` | DBD detection with hmmscan and construction of SR segments with hmmalign |
| `reference.py` | isoform selection for CisBP TFs, and loading of a built reference |
| `sr.py` | execution of the SimilarityRegression scoring code, and tier assignment |
| `meme.py` | writing, subsetting and validation of MEME motif files |

## constants.py

Values used by more than one module or script: the project layout (`ROOT`, `DATA_RAW`,
`REFERENCE_ROOT`, `RESULTS_ROOT`, `LOGS_DIR` and others), the file names of the
reference and of the per-organism results, the CisBP file names, the motif provenance
columns taken from CisBP, the SR class name `HSim` and the tier names `pass`,
`below-threshold` and `none`. Values used by a single module are defined in that module;
user-adjustable settings belong in `config/*.toml`.

## log.py

`configure_logging(log_path)` is called once at the start of each script. It directs the
records of all loggers to standard output and appends them to `log_path`. Modules log
through `logging.getLogger(__name__)`, so that each record identifies its module of
origin.

## files.py

- `load_toml(path)`: reads a configuration file; a bare file name is resolved in
  `config/`.
- `reference_config()`, `reference_dir(config)`, `results_dir(label)`: the shared
  configuration and the output directories.
- `fetch(url, dest)`: downloads a file into `data/raw/` once. The download is written
  under a temporary `.part` name, so that an interrupted transfer is never mistaken for a
  complete file, and its source, size and SHA-256 checksum are recorded in
  `data/raw/SOURCES.tsv`. Subsequent calls only verify the checksum; a modified file, or
  an existing file without a record, raises `FetchError`.
- `run(command)`: executes an external tool and raises an exception containing its
  error output if it fails.
- `write_provenance(out_dir, inputs, name, **extra)`: writes `<name>.json` containing the
  git commit of the code, the HMMER and gffread versions, the SHA-256 checksum of every
  input and any additional settings. Steps writing to the same directory use distinct
  names.

## cisbp.py

- `fetch_tables(config)`: retrieves the TF table and the motif archive from the CisBP
  bulk downloads.
- `fetch_proteins(config, name, tf_ids)`: protein sequences are not part of the CisBP
  bulk downloads and are therefore requested through the download form of the CisBP
  website (`archive.php`), in fixed batches of sorted identifiers. Because the form
  names its archive by the minute only, a returned archive is accepted only if every TF
  it contains was requested in the same batch.
- `families(tf_table)`: one row per family, with its DBD types, SR model and cutoff.
- `reference_tfs(tf_table)`: TFs with an experimentally determined motif (`TF_Status`
  `D`), excluding the family `Unknown`, which has no DBD for comparison.
- `check_targets(tf_table, species)`: the `I` and `N` TFs of the verification species,
  with the reference TF from which each `I` motif was transferred.
- `calling_dbd_types(tf_table, min_species)`: the DBD types CisBP uses to define TFs in
  at least `min_species` species.
- `load_matrices(pwm_zip)`: all motifs as width × 4 matrices (A, C, G, T), with rows
  renormalised to sum to 1; motifs listed without matrix rows are skipped.
- `representative_motifs(tf_table, ...)`: the motif designated for transfer from each
  reference TF, namely the motif CisBP most often transfers to other TFs, or otherwise
  the first motif listed.

## proteome.py

`load(target)` returns the proteins (`protein_id`, `gene_id`, `sequence`), the input
files used and summary counts for the provenance record.

- Genome and GFF3: gffread translates every coding transcript (`-y`) and reports the
  parent gene of each transcript (`--table @id,@geneid`). A transcript without a parent
  gene is treated as a gene of its own and counted.
- Protein FASTA: each record is treated as one gene.
- A terminal stop codon is removed; proteins with an internal stop codon are excluded
  and counted.

## hmmer.py

A single permissive hmmscan run (E-value ≤ 1) reports all hits. The TF identification
and segment rules are applied afterwards as filters, so the scan is never repeated.

- `build_library(...)`: extracts the named models from Pfam and indexes them with
  hmmpress.
- `scan(library, sequences)`: runs hmmscan and returns one row per domain hit.
- `gathering_pass(hits, models)`: retains hits that pass the Pfam gathering threshold as
  `hmmscan --cut_ga` applies it (sequence and domain score both at or above the
  threshold).
- `resolve_overlaps(hits)`: within each protein, retains the highest-scoring hit among
  overlapping hits.
- `segment_hits(...)` and `segments(...)`: construct the SR input. Each DBD hit of the
  family's types is aligned to its Pfam model with hmmalign, and the residues at the
  model's match states are extracted. Each segment spans all DBD types of the family in
  CisBP order: the block of the domain's own type contains its match states, and the
  blocks of the other types are gaps, as in the SR training code.
- `match_state_string(...)`: the SR alignments were produced with the R package aphid,
  which places every residue of a DBD segment. To reproduce them, residues that hmmalign
  leaves unaligned before the first (or after the last) occupied match state are moved
  into the adjacent empty match states. On the SR example, this yields strings identical
  to those of SR for 10 of 10 DBDs.

Reference TFs and target TFs are processed by the same functions, so that both sides of
every SR comparison are constructed identically.

## reference.py

- `resolve_tfs(...)`: CisBP lists all isoforms of a TF without indicating which one it
  scored. All isoforms are scanned; an isoform is valid if it has exactly the family's
  DBD types, and among valid isoforms the one with the most segments, then the greatest
  length, then the lowest identifier is selected. The result has one row per TF with its
  segments, or the reason none were obtained.
- `load(ref_dir)`: returns the built reference as the reference table, the family table,
  the DBD models, the usable reference TFs per family and the designated motif per
  reference TF.

## sr.py

SR is downloaded at a pinned commit and extracted to `reference/.../sr/`. Its scoring
code is executed without modification; two measures make it run under Python 3:

- `PairwiseAlignment.py` imports `Bio.SubsMat.MatrixInfo.blosum62`, which was removed in
  Biopython 1.80. A substitute module providing the same BLOSUM62 values (from
  `Bio.Align.substitution_matrices`) is registered before the import.
- `PredictSimilarity.py` mixes tabs and spaces in an unused helper function, which raises
  a `TabError` under Python 3. The two required functions, `ReadSRModel` and
  `ScoreAlignmentResult`, precede that helper, so the source text up to the helper is
  executed verbatim.

`family_models(...)` loads the SR model of each family and asserts that the model
belongs to that family and that its threshold equals the CisBP cutoff.
`Scorer.score_all(...)` scores the segments of one TF against a list of reference TFs.
`decide(scores, motifs)` assigns the tier: `pass` if an HSim reference TF has a motif
(the best such TF by SR score, then percent identity, then identifier), otherwise
`below-threshold`, or `none` if no reference TF could be scored.

## meme.py

Writes MEME minimal format, version 4, with a uniform background. Probabilities are
written with six decimal places and rounded such that each written row sums to exactly
1. `nsites` is set to a fixed weight of 20 rather than a site count, because CisBP
provides probabilities only. `validate` re-reads a written file and checks for unique,
filesystem-safe identifiers, the declared widths, four columns and rows summing to 1.
