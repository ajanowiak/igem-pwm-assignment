# cisbp-sr-motifs

Give the transcription factors (TFs) of any proteome a DNA-binding PWM, using CisBP's own
method: **Similarity Regression** (SR; Lambert et al. 2019) against **CisBP build 3.10**.
The output is one MEME file of PWMs, ready for FIMO, and a table saying for every TF
whether it got a PWM, from which lab-measured TF, and with what score.

**A PWM given to a TF here is a prediction of binding from protein similarity, not
evidence of binding.** A TF whose DNA-binding domain is not similar enough to any
lab-measured TF gets no PWM; nothing is forced.

## Words used here

| word | meaning |
|---|---|
| PWM | a position weight matrix: the DNA-binding preference of a TF. CisBP calls these "motifs" and gives each an id like `M01613_3.10`; file and column names here (`motifs.meme`, `motif_id`) follow CisBP |
| DBD | DNA-binding domain, the part of a TF that contacts DNA (homeobox, zinc finger, ...) |
| family | CisBP's TF family, defined by which DBD types a protein has (one homeobox: Homeodomain) |
| reference TF | a TF in CisBP (any species) with at least one PWM measured in a lab experiment; the only TFs PWMs are copied from |
| target TF | a TF of the organism being analysed, which gets its PWM copied from a reference TF |
| SR score | how similar two DBDs of one family are, as SR's family model scores them |
| cutoff | CisBP's SR score threshold for a family; at or above it, two DBDs count as "highly similar" (SR's HSim) |
| margin | SR score minus cutoff; the further above 0, the more trustworthy the copied PWM |

## How it works

```
A  build_reference.py   CisBP 3.10 (TF table, PWMs, protein sequences) + Pfam 31.0 + SR
                        -> reference/: reference TFs with SR-ready DBD segments and their PWMs
B  call_tfs.py          genome + GFF3 (gffread) or protein FASTA -> proteins -> hmmscan
                        -> target TFs with CisBP family and SR-ready DBD segments
   assign_motifs.py     SR score of each target TF vs every reference TF of its family
                        -> assignments.tsv + motifs.meme
   check_sr_vs_cisbp.py optional: reproduce SR's own example and CisBP's own PWM copies
```

**What SR is.** For one DBD family (homeodomain, Zinc cluster, ...), CisBP aligned the
DBDs of many pairs of TFs whose PWMs had been measured, position by position along the
family's Pfam HMM. It then fitted a linear model (ridge regression with non-negative
weights) that predicts from per-position identity (or BLOSUM62 similarity) whether two
TFs bind the same DNA. The result is, in effect, *percent identity with learned weights
per DBD position*, plus a family cutoff chosen for at least 75% precision on held-out
data. CisBP uses exactly this, from build 2.00 on, to copy PWMs from reference TFs to TFs
nobody has measured. SR's code and fitted models are public
(github.com/smlmbrt/SimilarityRegression); this pipeline runs them unmodified.

**The data flow.** SR compares two DBD sequences, so it needs the sequence of every
reference TF. CisBP serves those (every listed isoform per TF) through its cart download,
not its bulk downloads; `build_reference.py` fetches them in batches. Each reference TF
gets one isoform, its DBDs cut out and aligned to Pfam 31.0 HMM match states (SR's input
format), and one PWM. A target TF is scored against every reference TF of its CisBP
family; if the best one clears the family's cutoff, its PWM is copied to the target TF.

## Quick start

```bash
mamba env create -f environment.yml
```

```bash
conda activate cisbp-sr
```

```bash
python build_reference.py
```

```bash
python call_tfs.py --target tatroviride_p1.toml
```

```bash
python assign_motifs.py --target tatroviride_p1.toml
```

```bash
python check_sr_vs_cisbp.py
```

`build_reference.py` runs once for all organisms (about 8 minutes, mostly the Pfam
download). For another organism, copy `config/tatroviride_p1.toml` (RefSeq) or
`config/ncrassa.toml` (Ensembl), set the label and either `genome` + `annotation` (GFF3)
or `proteins`, as URLs or local paths.

- With a GFF3, proteins are grouped into genes by their parent gene feature. If the
  annotation has no gene features, each transcript counts as its own gene and
  `call_tfs.py` prints a warning (isoforms of one gene may then be counted twice).
- With a protein FASTA, each record is one gene; give one protein per gene if you can.
- gffread translates with the standard genetic code, so mitochondrial proteins usually
  show internal stops and are dropped (counted in the output). They are not TFs.

## Outputs

`reference/cisbp_3.10_pfam_31.0/` (shared by every organism):

| file | content |
|---|---|
| `references.tsv` | one row per reference TF: family, the PWM that gets copied from it and that PWM's experiment and study, the CisBP isoform used, its SR segments; `outcome` says why a TF has no segments |
| `motifs.meme` | every PWM named in `references.tsv` |
| `families.tsv` | per CisBP family: DBD types, SR model file, cutoff |
| `dbd_models.hmm`, `dbd_models.tsv` | the Pfam 31.0 models of CisBP's DBD types; which of them call TFs |
| `proteins.fasta`, `sr/`, `provenance.json` | chosen reference isoforms; SR code; input hashes and tool versions |

`results/<label>/`:

| file | content |
|---|---|
| `motifs.meme` | each copied PWM once (id = CisBP motif id, alt name = the reference TF it came from). Uniform background; `nsites= 20` is a fixed weight, not a site count; `E= 0` is a placeholder |
| `assignments.tsv` | one row per target TF gene: tier, reason, PWM, best reference TF and species, SR score, family cutoff, margin, runner-up, the PWM's source experiment |
| `tf_calls.tsv`, `dbd_hits.tsv` | target TF genes with family and SR segments; every hmmscan domain hit |
| `provenance_call_tfs.json`, `provenance_assign_motifs.json` | per step: code version, input SHA-256s, tool versions, settings used |

Tiers: **pass** (a reference TF scores at or above the cutoff and has a PWM; that PWM is
copied), **below-threshold** (scored, no reference TF qualifies; no PWM), **none**
(nothing to score).

## Every rule, and why

| rule | why |
|---|---|
| CisBP 3.10 is the only source of PWMs, reference sequences, families and cutoffs | one curated, versioned source; the cutoffs are the ones CisBP calibrated with SR |
| Pfam 31.0 | the only release whose model lengths equal the SR weight vectors of all 26 positional SR models; `build_reference.py` asserts it |
| SR at commit `b4b27b3`, unmodified | two Python 3 shims only (`srmotifs/sr.py`); reproduces SR's saved example scores to 5e-7 |
| A protein is a TF if a DBD hit passes Pfam's gathering threshold (GA) | Pfam's own curated per-model threshold, no tuned cutoff |
| ... for DBD types CisBP itself uses to call TFs in ≥ 10 species | excludes exactly RRM_1 and zf-CCCH (mostly RNA-binding; 1 and 2 species in CisBP); all other types are in ≥ 22 |
| Family = the CisBP family whose DBD types equal the protein's | CisBP's own family definitions |
| SR segments: the family's DBD types at domain i-Evalue ≤ 1e-3, aligned with hmmalign | same rule for reference and target TFs; hmmalign plus SR's end placement gives the same match states as SR's own (aphid) alignments, 10/10 on SR's example |
| Reference isoform: matches the family's DBD types, then most DBDs, longest, smallest id | CisBP lists isoforms without marking the one it scored |
| One PWM per reference TF: the one CisBP itself copies to other TFs most often | a reference TF can have several lab-measured PWMs (different studies); this follows CisBP's own choice |
| Threshold: CisBP's family cutoff, unchanged | no retuning; a below-threshold TF is *unknown*, not *different* |

## Verified (2026-10-10)

*T. atroviride* P1 (NCBI RefSeq GCF_020647795.1):
- gffread translations equal NCBI's own protein FASTA for 13,567 of 13,568 proteins.
- 431 TF genes; 109 pass (101 distinct PWMs), 321 below-threshold, 1 none.
- With the gene features removed from the GFF, the run warns and counts 457 TF
  "genes" (one per transcript), the same as with NCBI's protein FASTA as input.

*N. crassa* OR74A (Ensembl Fungi 63), a species CisBP covers itself:
- 231 TF genes, 228 of them in CisBP's own *N. crassa* TF list (290); the family agrees
  for 224 of the 228. 22 proteins dropped for internal stops are all mitochondrial.
- 125 of CisBP's 156 lab-measured *N. crassa* TFs are called; 121 of those get their own
  measured PWM back (their best reference TF is themselves). Of the other 4, two are
  reference TFs CisBP lists no protein for, and two get a different DBD combination
  than CisBP records.
- Against CisBP's own calls for these genes (I = CisBP copied a PWM to it, N = CisBP did
  not): I → pass 10 of 12, N → no pass 87 of 91.
- 61 of the 62 CisBP TFs not called have a DBD hit below Pfam's gathering threshold
  (see Known limitations).

Both organisms:
- 4,745 of 4,817 CisBP reference TFs are usable (60 have no isoform with the family's DBD
  types under the segment rule; for 12 CisBP lists no protein).
- `check_sr_vs_cisbp.py`, on CisBP's own I/N calls for 6 fungi (1,318 TFs): CisBP I →
  our pass 94.6%, CisBP N → our no pass 91.7%. Disagreements sit near the cutoff: our pass
  calls agree with CisBP 47% within 0.025 of the cutoff and 100% beyond 0.2.

## Assumed, not verified

- That the isoform chosen for each reference TF is the one CisBP scored.
- That CisBP's I/N labels reflect SR on today's reference set (disagreements cluster at
  the cutoff, consistent with isoform or segment-boundary differences).
- That a PWM copied by SR describes the target TF's binding. SR's own held-out precision
  for HSim is about 0.8; it is weakest for Zinc-cluster TFs.

## Known limitations

- TF calling with Pfam's gathering threshold is stricter than CisBP's own TF calling.
  On *N. crassa* it misses 62 of CisBP's 290 TFs (31 lab-measured); calling at the
  segment rule's i-Evalue ≤ 1e-3 would find 282 of 290, with 36 calls CisBP lacks.
- CisBP's protein table lacks 12 reference TFs whose PWMs CisBP itself copies to other TFs.
- For C2H2 zinc-finger proteins, the segment rule often keeps fewer fingers than CisBP lists.
- Fungal reference TFs are sparse: most copied PWMs come from *N. crassa*,
  *A. nidulans*, *S. cerevisiae* and *F. graminearum*; most Zinc-cluster TFs get none.

## Citation

Lambert SA et al. (2019) Similarity regression predicts evolution of transcription factor
sequence specificity. *Nature Genetics* 51:981–989.
Weirauch MT et al. (2014) Determination and inference of eukaryotic transcription factor
sequence specificity. *Cell* 158:1431–1443.
