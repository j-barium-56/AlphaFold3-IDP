# AlphaFold 3 on IDR–IDR pairs

Given two intrinsically disordered regions and nothing else, does AlphaFold 3 co-folding
distinguish pairs that really interact from pairs that do not?

Not "is the predicted structure right" — for most of these pairs there is no structure to
be right about — but pure **discrimination**, scored on the confidence outputs AF3 already
produces, against a trivial sequence-composition baseline rather than against chance.

> **Stage-1 answer, on 60 of 62 jobs: no.** ipTM does not beat sequence composition on any
> pre-registered test, and in the homotypic arm it is significantly worse. The full result,
> with every grouping and both confounds, is in
> [`af3_idr_results_report.pdf`](af3_idr_results_report.pdf) (26 pages).

## Why the question is open

The closest paper — Mehdiabadi et al., *AlphaFold3 and Intrinsically Disordered Proteins*
(bioRxiv Dec 2025 / Protein Science 2026) — benchmarks AF3 on IDR–**receptor** dimers and
fold-upon-binding complexes. It explicitly does not cover pairs where **both** chains are
disordered, disorder-to-disorder ("fuzzy") complexes, homotypic self-association, or phase
separation, and its headline metric (DockQ) needs a reference structure that does not exist
for most pairs here.

## The dataset

**62 jobs: 21 positives, 31 real-sequence negatives, 10 composition-matched scrambles**,
drawn from 30 IDR regions. 26 of the 30 were verified against UniProt's reported length
*and* average molecular mass when the set was built; FMR1 and the three late additions
(p53 TAD, p21, 4E-BP1) were later checked character-for-character against UniProt. 30/30.

| Positives | n | | Negatives | n |
|---|---|---|---|---|
| Heterotypic, fuzzy | 5 | | `noncognate` — two real IDRs, no reported interaction | 10 |
| Heterotypic, mutual folding | 2 | | `cross_kingdom` — human × plant/yeast | 5 |
| Homotypic, fuzzy, LLPS | 12 | | `charge_repulsive` — two strongly acidic IDRs | 2 |
| Homotypic, mutual folding | 2 | | `hard_coiledcoil` — zipper-competent, do not pair | 3 |
| | | | `homo_negative` — IDRs that do not self-associate | 11 |
| | | | `scramble` — shuffled chain B, **separate arm** | 10 |

Negatives carry a per-row `evidence_strength` grade rather than being assumed uniform:
8 `demonstrated` (3 of them `matched_internal`, where the same chain also appears in a
positive), 11 `inferred`, 10 `absence_only`, 2 `uncertain`. The weakest joint in a
benchmark like this is how you know a negative is negative, so it is graded, not hidden.

## Stage-1 result

Primary metric is ipTM of the top-ranked model. The comparator is a leave-one-out logistic
model on cheap composition features — **not 0.5** — and a test passes only if the paired
bootstrap CI on ΔAUC excludes zero.

| Test | AF3 ipTM AUC | composition | ΔAUC (95% CI) | |
|---|---|---|---|---|
| 1. positives vs 31 real negatives | 0.611 | 0.767 | −0.156 [−0.32, +0.01] | fail |
| 2a. mutual folding vs hard zipper | 0.417 | 0.250 | +0.167 [−0.50, +0.88] | fail |
| 2b. fuzzy positives — *the result* | 0.676 | 0.784 | −0.107 [−0.26, +0.05] | fail |
| 3. vs gold-standard negatives | 0.118 | 0.640 | −0.522 [−0.74, −0.29] | fail |
| Homotypic arm, arm baseline | 0.564 | 0.902 | −0.337 [−0.58, −0.12] | fail |
| Heterotypic arm, arm baseline | 0.485 | 0.500 | −0.015 [−0.32, +0.29] | fail |

Four things that matter for reading those numbers:

- **The pipeline is sound.** The GCN4-p1 canary reaches ipTM 0.73, all four memorised
  controls sit at 0.70–0.80, and all 60 jobs ran the sequences they were supposed to.
  Test 2a fails because the three hard-zipper negatives score 0.70/0.78/0.80 — AUC 0.000
  against the positives. AF3 dimerises any zipper-competent pair, cognate or not.
- **Test 3 cannot adjudicate.** Its 8 negatives are 94–236 residues and its 17 positives
  205–882, so total chain length alone separates that split at AUC 0.985. Length-adjusting
  moves ipTM from 0.118 to 0.449.
- **Nothing in the 16-metric panel rescues it.** Interchain contact-probability sum is the
  one metric that separates on its own (0.787, q = 0.008) and it still does not beat
  composition (Δ = +0.020). Of 464 cells with a baseline, 15 clear zero and 14 are the
  memorised controls.
- **Templates were on for all 95 chain entries** — the server default, never a choice made
  here — and they land on both the memorised positives and the matched-internal negatives.
  The rerun that tests this is the first item in the handoff.

What did work: the three matched-internal contrasts, which rest on no literature absence at
all. AF3 ranks ACTR×NCBD above ACTR×ACTR and NCBD×NCBD, and c-Fos×c-Jun above c-Fos×c-Fos —
3/3 on ipTM. That is the design worth expanding.

## Reproducing it

```bash
python3 baseline_leakage.py af3_idr_pairs_1.csv        # run BEFORE looking at AF3 output
# upload af3_batch_01..04.json to alphafoldserver.com, download results into data/
python3 parse_af3_results.py data af3_idr_pairs_1.csv  # -> af3_scored.csv + the tests
python3 af3_report_stats.py                            # -> af3_report_data.json (~1 min)
python3 af3_report_build.py                            # -> the 26-page PDF + full CSV
python3 af3_workbook_build.py                          # -> af3_idr_workbook.xlsx, every raw + calculated value
python3 viewer/serve.py                                # -> http://127.0.0.1:8765, PAE / contact maps per sample
```

Needs `numpy` (parser, baseline), `reportlab` (report) and `openpyxl` (workbook). The report registers Arial from
`/System/Library/Fonts/Supplemental`, so as written it builds on macOS; both report scripts
currently hard-code the repository path. `parse_af3_results.py` accepts either the result
zips or unzipped job folders, and matches them to rows by job name — check the list of
unmatched jobs it prints before trusting any AUC.

## Layout

| | |
|---|---|
| `af3_idr_pairs_1.csv` | The dataset: one row per job, both sequences, all labels, `evidence_strength`, stage |
| `af3_idr_regions.csv` | 30 IDR regions: accession, boundaries, sequence, boundary confidence |
| `af3_batch_01..04.json` | Upload-ready AlphaFold Server batch files (20/20/20/2) |
| `baseline_leakage.py` | Composition-only baseline and single-feature leakage check |
| `parse_af3_results.py` | Result folders → scored table, pre-registered tests, metric panel |
| `af3_report_stats.py` | The same functions over every grouping × metric → JSON |
| `af3_report_build.py` | JSON → the PDF (vector, no hand-typed numbers) |
| `af3_idr_results_report.pdf` | Stage-1 results |
| `af3_all_statistics.csv` | All 496 cells: AUC, CI, permutation p, BH q, ΔAUC and its CI |
| `af3_scored.csv` | Per-job scored table |
| `af3_workbook_build.py` | Raw server output + every calculated table → one workbook; refuses to write unless it reproduces `af3_scored.csv` |
| `af3_idr_workbook.xlsx` | That workbook. `Jobs` is the pair table with raw and calculated columns appended, `Models` one row per server sample |
| `viewer/` | `serve.py` + `index.html`: local job browser with per-sample PAE and contact-probability heatmaps, read from `data/` |
| `af3_idr_plan.md` | Scientific design, pre-registered criteria, the stage-1 result, limitations |
| `af3_idr_handoff.md` | Operational doc: what to run next, and what not to undo |

## What is not in the repository

`data/` — the raw AlphaFold Server result folders, about 493 MB across 60 jobs (5 models,
5 summary files and 5 full-data files each, plus MSAs and template hits). It is gitignored;
re-download from the server if it is lost. Everything the analysis needs is derived from it
into `af3_scored.csv`, `af3_report_data.json` and `af3_all_statistics.csv`, all of which are
tracked.

## Status and caveats

Stage 1 is scored; batch 4 (2 scramble jobs) is outstanding. The honest limitations live in
[`af3_idr_plan.md`](af3_idr_plan.md) — 15 of them, including that 21 of 31 negatives rest on
inferred or absence-only evidence, that length is not matched by design, that the ΔAUC
intervals are uncorrected for multiplicity, and that five samples per job are not five
seeds. Next steps, in priority order, are in [`af3_idr_handoff.md`](af3_idr_handoff.md).

Research in progress; nothing here is published.
