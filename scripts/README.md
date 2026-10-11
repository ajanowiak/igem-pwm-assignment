# scripts

Command-line scripts of the pipeline. Each script performs one step; `run_pipeline.sh`
in the repository root executes them in order. The purpose of the pipeline and its
installation are described in the [root README](../README.md); the shared code is
documented in [srmotifs/](srmotifs/README.md).

All scripts are run from the repository root within the conda environment:

```bash
conda activate igem-pwm-assignment
```

| order | script | frequency | input | output |
|---|---|---|---|---|
| 1 | `build_reference.py` | once | `config/reference.toml`, downloaded data | `reference/cisbp_3.10_pfam_31.0/` |
| 2 | `call_tfs.py` | per organism | target configuration, reference | `results/<label>/tf_calls.tsv`, `dbd_hits.tsv` |
| 3 | `assign_motifs.py` | per organism | `tf_calls.tsv`, reference | `results/<label>/assignments.tsv`, `motifs.meme` |
| – | `check_sr_vs_cisbp.py` | optional | reference, CisBP's own inferences | `results/checks/sr_vs_cisbp/` |

## Contents

- [Script descriptions](#script-descriptions)
- [Configuration](#configuration)
- [Outputs](#outputs)
- [Method rules](#method-rules)
- [Logging, downloads and provenance](#logging-downloads-and-provenance)

## Script descriptions

### build_reference.py

```bash
python scripts/build_reference.py
```

Downloads CisBP 3.10 (TF table, motifs and the protein sequences of the required TFs),
Pfam 31.0 and SimilarityRegression at a pinned commit, and builds the reference against
which all organisms are scored:

1. extracts the Pfam models of all DBD types used by CisBP and determines which of them
   are used for TF identification (see [Method rules](#method-rules));
2. pairs each CisBP family with its SR model, asserting that the model's weight vector
   spans exactly the match states of the family's Pfam models and that the SR threshold
   equals the CisBP cutoff;
3. for each reference TF (a CisBP TF with an experimentally determined motif), scans all
   isoforms listed by CisBP, selects one and aligns its DBDs to the Pfam models (SR
   segments);
4. designates one motif per reference TF for transfer and writes all designated motifs
   to a single MEME file.

The log reports the number of usable reference TFs and the reasons the remainder are
unusable. The first run takes several minutes, mainly for the Pfam download and the
batched retrieval of CisBP protein sequences; subsequent runs reuse the downloaded files.

### call_tfs.py

```bash
python scripts/call_tfs.py --target tatroviride_p1.toml
```

Identifies the TF genes of one organism:

1. obtains the proteome, either by translating genome and GFF3 with gffread or from a
   protein FASTA file;
2. scans all proteins with `hmmscan` against the DBD models of the reference;
3. classifies a protein as a TF if it has a DBD hit of a TF-defining type that passes the
   Pfam gathering threshold; of overlapping hits, the highest-scoring one is retained;
4. selects one representative protein per gene;
5. assigns the CisBP family whose set of DBD types equals that of the protein;
6. converts the protein's DBDs into SR segments, using the same code as for reference
   TFs.

The log reports the numbers of proteins, genes and TF genes and the most frequent
families, and issues a warning if proteins lack a gene identifier.

### assign_motifs.py

```bash
python scripts/assign_motifs.py --target tatroviride_p1.toml
```

Scores the SR segments of each TF gene against all usable reference TFs of its family
with the family's SR model, and assigns each gene to one tier (`pass`,
`below-threshold` or `none`; see the [root README](../README.md#outputs)). Before the
results are written, the script asserts that every TF gene is assigned exactly one tier,
that every scored gene has a numeric SR score and that only `pass` genes carry a motif.
The written MEME file is re-read and validated.

### check_sr_vs_cisbp.py

```bash
python scripts/check_sr_vs_cisbp.py
```

Verifies the implementation of the method; it does not process any target organism.
Part 1 runs the homeodomain example distributed with SR and requires the 55 published
pairwise scores to be reproduced (to within 1e-6, with identical classifications) and
the DBD alignments to be identical to those of SR; otherwise the script exits with an
error. Part 2 selects the TFs of the species listed under `[check]` in
`config/reference.toml` for which CisBP inferred a motif (status `I`) or declined to
infer one (status `N`), processes them in the same way as reference TFs, scores them and
compares the outcome with CisBP's. Agreement indicates that the method is implemented as
CisBP applies it; it is not evidence that any individual motif is biologically correct.

### Common options

| option | scripts | description |
|---|---|---|
| `--target` | `call_tfs.py`, `assign_motifs.py` | target configuration, as a path or as a file name in `config/` |
| `--log_path` | all | log file to append to; defaults listed below |

Default log files: `logs/build_reference.log`, `logs/call_tfs/<label>.log`,
`logs/assign_motifs/<label>.log`, `logs/check_sr_vs_cisbp.log`.

## Configuration

### config/reference.toml

Pins all settings shared by the organisms. Any change to this file affects all results.

| key | description |
|---|---|
| `cisbp.build`, `cisbp.bulk_url` | CisBP build and its bulk-download directory |
| `cisbp.site_url`, `cisbp.protein_batch_size` | CisBP download form used for protein sequences (which are not part of the bulk downloads), and the number of TFs requested per batch |
| `pfam.release`, `pfam.hmm_url` | Pfam release and the location of its `Pfam-A.hmm.gz` |
| `similarity_regression.commit`, `.tarball_url` | SR commit, retrieved as a GitHub tarball |
| `rules.min_cisbp_species` | minimum number of species in which CisBP uses a DBD type to define TFs for that type to be used for TF identification (10) |
| `rules.calling_i_evalue` | 0 (default): a TF requires a DBD hit passing the Pfam gathering threshold. Values above 0: genes without such a hit may additionally be identified from hits with an independent E-value at or below this value. Genes identified by the gathering threshold are never altered |
| `rules.segment_max_i_evalue` | once the family of a protein is known, its SR segments are the hits of the family's DBD types with an independent E-value at or below this value (1e-3); applied identically to reference and target TFs |
| `check.species` | CisBP species used by `check_sr_vs_cisbp.py` |

`calling_i_evalue = 0` is the conservative setting, chosen because an incorrectly
assigned motif is more costly in downstream analyses than an unidentified TF. For
*N. crassa*, a value of 1e-5 adds 27 TF genes (17 of which have an experimentally
determined motif and recover it) without altering any existing assignment.

### Target configurations

One file per organism, for example `config/tatroviride_p1.toml`:

| key | description |
|---|---|
| `organism.name`, `organism.strain` | species and strain, recorded in the logs |
| `organism.label` | name of the output directory `results/<label>/` |
| `input.genome`, `input.annotation` | genome FASTA and GFF3 annotation (URL or path, optionally gzipped) |
| `input.proteins` | alternative to the above: protein FASTA with one record per gene |

URLs are downloaded once into `data/raw/targets/<label>/`. Relative local paths are
resolved from the repository root.

## Outputs

### reference/cisbp_3.10_pfam_31.0/

| file | content |
|---|---|
| `references.tsv` | one row per reference TF (columns below) |
| `motifs.meme` | all motifs listed in `references.tsv` |
| `families.tsv` | per CisBP family: `Family_ID`, `Family_Name`, `DBDs` (DBD types in CisBP order), `sr_model` (SR model file), `cutoff` |
| `dbd_models.hmm` (with index files), `dbd_models.tsv` | Pfam models of the CisBP DBD types; per model: `in_pfam`, `length` (number of match states), `cisbp_species`, `calls_tfs` |
| `proteins.fasta` | selected isoform of each reference TF |
| `sr/` | extracted SR code and models |
| `provenance.json` | code and tool versions, checksums of all downloaded inputs, CisBP TFs without a listed protein sequence |

`references.tsv`:

| column | description |
|---|---|
| `TF_ID`, `TF_Name`, `TF_Species`, `Family_ID`, `Family_Name` | as in CisBP |
| `motif_id` | motif designated for transfer from this TF (empty if CisBP provides no matrix for any of its motifs) |
| `protein_id` | CisBP isoform used; empty if no isoform is usable |
| `n_isoforms` | number of isoforms listed by CisBP |
| `n_segments`, `segments` | SR segments: one match-state string per DBD, comma-separated |
| `outcome` | `resolved (...)`, or the reason no segments were obtained |
| `motif_experiment`, `motif_study`, `motif_year`, `motif_pmid` | experiment type and publication of the motif, as in CisBP |

### results/\<label\>/

`tf_calls.tsv`, one row per TF gene:

| column | description |
|---|---|
| `gene_id`, `protein_id` | gene and its representative protein |
| `n_dbd_isoforms` | number of isoforms of the gene with a TF-defining DBD |
| `protein_length` | length of the representative protein |
| `n_dbds`, `architecture` | number of DBDs and their types in N- to C-terminal order |
| `dbd_types` | distinct DBD types, sorted; these determine the family |
| `family_id`, `family_name` | CisBP family; empty if no family has exactly these DBD types |
| `n_segments`, `segments` | SR segments, as in `references.tsv` |

`assignments.tsv`, one row per TF gene:

| column | description |
|---|---|
| `gene_id`, `protein_id`, `family_id`, `family_name`, `n_segments` | as in `tf_calls.tsv` |
| `tier`, `reason` | assignment tier and the reason for it |
| `motif_id` | assigned CisBP motif (`pass` only) |
| `best_reference_tf`, `best_reference_name`, `best_reference_species` | for `pass`, the reference TF from which the motif was transferred; otherwise the highest-scoring reference TF |
| `sr_score`, `sr_class` | SR score and SR classification of that reference TF (`HSim` highly similar, `Amb` ambiguous, `Dis` dissimilar) |
| `pct_id` | percent identity of the DBD alignment as reported by SR (`PctID_L`) |
| `cutoff`, `margin` | family cutoff, and `sr_score - cutoff` |
| `n_references`, `n_hsim` | number of reference TFs scored, and number classified as HSim |
| `runner_up_tf`, `runner_up_score` | highest-ranked reference TF other than `best_reference_tf`, and its SR score |
| `motif_experiment`, `motif_study`, `motif_year`, `motif_pmid` | experiment type and publication of the assigned motif (`pass` only) |

`dbd_hits.tsv` lists all domain hits reported by hmmscan at E-value ≤ 1 (`protein`,
`model`, sequence and domain scores, independent E-value, alignment start and end)
before the identification and segment rules are applied. `provenance_call_tfs.json`
additionally records proteome statistics (proteins, genes, proteins excluded for
internal stop codons, proteins without a gene identifier) and `calling_i_evalue`.

### results/checks/sr_vs_cisbp/

`summary.txt` (the report also written to the log), `comparison.tsv` (one row per CisBP
`I`/`N` TF with SR segments, listing CisBP's decision alongside the pipeline's) and
`agreement_by_margin.tsv` (fraction of the pipeline's `pass` calls with which CisBP
agrees, binned by margin above the cutoff).

## Method rules

All decision rules applied by the pipeline, with their justification.

| rule | justification |
|---|---|
| CisBP 3.10 is the sole source of motifs, reference sequences, families and cutoffs | a single curated, versioned source; the cutoffs are those CisBP calibrated with SR |
| Pfam 31.0 | the only release whose model lengths equal the SR weight vectors of all 26 positional SR models; asserted by `build_reference.py` |
| SR at commit `b4b27b3`, unmodified | only two Python 3 compatibility measures (`srmotifs/sr.py`); reproduces the published SR example scores to within 5e-7 |
| A protein is a TF if a DBD hit passes the Pfam gathering threshold | the curated per-model threshold provided by Pfam; no threshold is tuned |
| ... restricted to DBD types that CisBP uses to define TFs in at least 10 species | excludes exactly RRM_1 and zf-CCCH (predominantly RNA-binding; used in 1 and 2 species, respectively); all other types are used in 22 or more species |
| One protein per gene: most DBDs, then highest summed domain score, then greatest length, then lowest identifier | a fixed, deterministic order favouring the isoform with the most DBDs |
| Family: the CisBP family whose DBD types equal those of the protein | CisBP's own family definitions |
| SR segments: hits of the family's DBD types at independent E-value ≤ 1e-3, aligned with hmmalign | identical rule for reference and target TFs; combined with SR's terminal residue placement, it yields the same match states as SR's own alignments (10 of 10 on the SR example) |
| Reference isoform: exactly the family's DBD types, then most DBDs, then greatest length, then lowest identifier | CisBP lists isoforms without indicating which one it scored |
| One motif per reference TF: the motif CisBP most often transfers to other TFs (otherwise the first listed) | a reference TF may have several experimentally determined motifs; this follows CisBP's own choice |
| Threshold: CisBP's family cutoff, unchanged | no retuning; a `below-threshold` TF has an unknown, not a demonstrably different, specificity |

## Logging, downloads and provenance

- **Logging.** Each script writes its log to standard output and appends it to a file in
  `logs/` (format: time, module, level, message).
- **Downloads.** Each downloaded file is stored once in `data/raw/` and recorded in
  `data/raw/SOURCES.tsv` with its source, retrieval time, size and SHA-256 checksum.
  Subsequent runs verify the checksum; a modified or unrecorded file aborts the run.
- **Provenance.** Each output directory contains a `provenance*.json` file with the git
  commit of the code (suffixed `-dirty` if tracked files were modified), the HMMER and
  gffread versions, the SHA-256 checksum of every input and the settings used.
- **Reproducibility.** All steps are deterministic. A rerun from `data/raw/` reproduces
  `reference/` and `results/` byte for byte, except for the timestamps in the provenance
  files.
