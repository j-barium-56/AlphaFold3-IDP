# AF3 IDR–IDR benchmark — handoff

**Status:** stage 1 scored. 60 of 62 jobs run, parsed and reported; batch 4 (2 jobs) pending.
**Written:** 14 September 2026, updated the same evening with the stage-1 result.
**Companion doc:** `af3_idr_plan.md` (the scientific design and the numbers in full).
This document is the operational one — enough to pick the project up cold.

---

## 1. What the project is, in three sentences

Given two disordered regions and nothing else, does AlphaFold 3 co-folding distinguish
pairs that really interact from pairs that do not? Not "is the structure right" — for most
of these pairs no structure exists — but pure **discrimination**, scored on AF3's own
confidence outputs. The gap is real: the closest paper (Mehdiabadi et al., bioRxiv Dec 2025 /
Protein Science 2026) covers IDR–**receptor** pairs and fold-upon-binding, and explicitly
not IDR–IDR, not fuzzy complexes, not homotypic, not phase separation.

## 2. Current state

| Thing | State |
|---|---|
| 62 curated jobs, labelled | **done** |
| All sequences verified against UniProt | **done**, 30/30 — 26 by length + mass, then FMR1 and the three late additions (p53 TAD, p21, 4E-BP1) character-exact |
| AlphaFold Server batch JSONs | **done**, 4 files, names validated |
| Composition-baseline leakage check | **done**, run 14 Sep: 0.700 all / 0.529 heterotypic / 0.851 homotypic |
| AlphaFold Server submissions | **60 of 62** — batches 1–3 done, batch 4 (S20, S21) pending |
| Result parser + statistics | **done**, 496 cells over 31 groupings × 16 metrics |
| Results report | **done**, `af3_idr_results_report.pdf` |
| Stage 2 branch decision | **open** — see §4 |

## 3. Files

| File | Role |
|---|---|
| `af3_idr_pairs_1.csv` | The dataset. One row per job: both sequences, all labels, `evidence_strength`, stage |
| `af3_idr_regions.csv` | 30 IDR regions: accession, boundaries, sequence, `boundary_confidence`. p53 TAD, p21 and 4E-BP1 were back-filled after stage 1; the histone id is `H1.0_full` here and `H10_full` in the pair table |
| `af3_batch_01..04.json` | Upload-ready AlphaFold Server batch files (20/20/20/2) |
| `baseline_leakage.py` | Composition-only baseline. **Run before AlphaFold** |
| `parse_af3_results.py` | Result zips → `af3_scored.csv`, ΔAUC tests, corrected metric panel |
| `af3_report_stats.py` | Every grouping × metric → `af3_report_data.json` (uses the parser's own functions) |
| `af3_report_build.py` | `af3_report_data.json` → the 26-page PDF |
| `af3_idr_results_report.pdf` | Stage-1 results |
| `af3_all_statistics.csv` | All 496 cells, machine-readable |
| `af3_scored.csv` | Per-job scored table |
| `af3_idr_plan.md` | Scientific design, pre-registered criteria, stage-1 result, limitations |
| `af3_idr_handoff.md` | This file |
| `data/` | Raw result folders, ~493 MB, **gitignored** — re-download from the server if lost |

Build scripts (`build_dataset.py`, `verify.py`, `raw_fetch_1.fasta`) live in the session
workspace and are only needed if the set is rebuilt. Ask for them if so.

## 4. Do this next — 20 jobs, one day

1. **Submit batch 4** (S20, S21). Finishes the scramble arm at 10 pairs. The paired test is
   currently 6/8 with a sign-test p of 0.29 — two more pairs will not fix the power, but it
   completes the arm as designed.
2. **Templates off.** Re-run P06, P07, P20, P21, N15, N16, N17, N26, N27, N28 with
   `"useStructureTemplate": false` on every chain — 10 jobs. **This is the highest-value
   experiment left**, because templates were on for all 95 chain entries by server default
   and they land squarely on the memorised positives *and* the matched-internal negatives.
   If the controls collapse without templates, test 2a's failure is about templates rather
   than about AF3's priors, and the memorisation caveat changes shape.
3. **Seed replicates.** `"modelSeeds": ["2","3"]` on P06, P07, P20, P21 only — 8 jobs.
   Seed 1 already ran on 59 of 60 jobs, so seeds 2 and 3 are the only new information.
   This is the first actual measurement of run-to-run variance.
4. **Confirm or drop** the two `uncertain` negatives: H1.0 alone (N22) and ERD10 alone
   (N24). They currently sit in every headline number on unconfirmed labels.
5. **Decide how test 3 handles length, and write it down before looking again.** Either a
   length-matched subset or the partialled version the parser already asks for. Deciding
   after seeing the numbers is the thing the pre-registration exists to prevent.

## 5. What stage 1 found

Full numbers in the PDF; the short version:

- **ipTM loses to sequence composition on every pre-registered test.** Test 1 Δ = −0.156
  [−0.32, +0.01]; test 2b, the pre-registered "this is the result", Δ = −0.107 [−0.26, +0.05].
- **The homotypic arm is significantly worse than the trivial baseline**: Δ = −0.337
  [−0.58, −0.12], interval clear of zero on the wrong side. In the heterotypic arm
  composition is at chance (0.500) and AF3 is no better (0.485), on 5 positives.
- **The pipeline is sound.** GCN4-p1 canary at ipTM 0.73, controls 0.70–0.80, all 60 jobs
  ran the right sequences with 5 samples each. Test 2a fails because the three hard-zipper
  negatives score 0.70/0.78/0.80 — AUC 0.000 against the positives. AF3 pairs any two
  zippers.
- **Test 3's inversion is a length artefact.** Length alone separates that split at 0.985;
  length-adjusting moves ipTM from 0.118 to 0.449.
- **One non-control cell in 464 beats composition**: fuzzy homotypic positives on
  interchain contact-probability sum, Δ = +0.091 [+0.011, +0.204]. Uncorrected across 464
  comparisons, and the pre-registered arm does not pass. Hypothesis, not result.
- **The matched-internal contrasts worked**: 3/3 correct orderings on ipTM. That is the
  design to expand — with the template caveat above.

## 6. Decisions already made, and why — do not silently undo these

**Sequences are mass-verified.** The only network route available in the build environment
relays text through a language model, which silently garbled the FMR1 sequence three
separate times in three different places. Every sequence was therefore checked against
UniProt's reported length *and* average molecular mass — 26/26 pass. FMR1 445–632 has since
been confirmed character-for-character against Q06787, and so have the three homotypic
negatives added late (p53 TAD = P04637 1–61, p21 = P38936 full-length, 4E-BP1 = Q13541
full-length), which had never been recorded in `af3_idr_regions.csv` and so were outside the
original tally. 30/30 now. **If you add any sequence, verify it the same way — and add it to
the regions table, which is where the last three went missing.**

**Scrambles are not the primary negative.** A shuffled sequence has no homologs, so AF3
builds it no MSA — measured: 4 to 34 unpaired sequences for a shuffled chain against 47 to
14,534 for a real one. Low confidence would then follow from missing evolutionary
information rather than a missing interface. The primary negatives are real IDRs in
non-cognate pairings. The scramble arm asks a different question and is reported separately.

**The comparator is the composition baseline, not 0.5.** Composition alone separates these
labels at AUC 0.851 in the homotypic arm — LC domains that self-associate genuinely *are*
aromatic-rich and weakly charged, which is biology and cannot be curated away. The primary
test is ΔAUC = AF3 − baseline with a paired bootstrap CI excluding zero. Stage 1 is exactly
the case this guards against: several metrics beat chance and none beats composition.

**Homotypic and heterotypic are analysed separately**, each against a baseline refitted
inside that arm (0.851 vs 0.471 raw). Different MSA situations; pooling lets chain identity
leak into the answer.

**ipTM stays the pre-registered primary even though it may be useless.** The 15 secondary
metrics exist because ipTM is likely floored for fuzzy pairs, but choosing the best one
after seeing data is p-hacking. The panel is reported with permutation p-values and
Benjamini–Hochberg correction; a secondary metric passing at q < 0.05 is stated as exactly
that — corrected, hypothesis-generating, needs a fresh set to confirm.

**Four positives are memorised** (P06/1KBH, P07/1FOS, P20/1AN2, P21/2ZTA). Pipeline
controls, excluded from the headline. **GCN4-p1 (P21) is the canary**: it passed at 0.73,
so a future failure there means the harness broke, not that IDRs got harder.

**Job names are sanitised to `[A-Za-z0-9_-]`.** The server rejects anything else — a dot in
`H1.0` broke the first job of batch 1 in the first version. If you add jobs, sanitise.

**Numbers in the report are generated, never typed.** `af3_report_stats.py` calls the
parser's own `auc`, `bootstrap_auc`, `perm_p`, `benjamini_hochberg` and
`delta_auc_vs_baseline`, plus `baseline_leakage.loo_logistic`, so the report cannot drift
from the pipeline. Prose in the PDF that quotes a number computes it from the JSON. **Do
not hand-edit a figure into that document** — change the data and rebuild.

**The length-adjusted analysis is a sensitivity check, not a pre-registered test**, and is
labelled that way in the report. Same for the negative-subtype, evidence-grade and
positive-subgroup groupings, and the paired scramble test.

## 7. Open questions and known weaknesses

Ordered by how likely they are to matter.

1. **Templates were on for all 95 chain entries**, by server default — the batch files are
   dialect v1 and never set `useStructureTemplate`, and the server upgraded them to v3 with
   it `true`. 6es7 templates both chains of P06 *and* N26 and N27; 1fos templates P07 *and*
   N15, N16, N17, N28. Until §4 step 2 runs, every memorisation statement is provisional.
2. **21 of 31 negatives rest on `inferred` or `absence_only` evidence.** "No paper reports
   this interaction" is not "these do not interact". Stage 1 makes this concrete rather than
   resolving it: the gold-standard subset inverts, but it is length- and foldability-
   confounded, so it cannot adjudicate. The structural fix is still to expand the
   `matched_internal` design — the one contrast that behaved — not to find more proteins
   nobody has reported binding.
3. **Two negatives are flagged `uncertain`**: H1.0 alone and ERD10 alone. Confirm or drop.
   N24 (ERD10 alone) is also the one job that ran on a server-assigned seed.
4. **Length is not matched by design.** Now quantified: length alone scores 0.985 on the
   gold split, 1.000 against the hard zippers, 0.676 on the headline. A second set should
   be length-matched at selection.
5. **A dimer is not a condensate.** LLPS is many-body and concentration-dependent; a 1:1
   co-fold cannot in principle report on it. The `llps` column is stratifying metadata,
   **never an outcome variable**. Stage 1 stratified by it and found nothing (0.604 with a
   reported LLPS phenotype vs 0.661 without) — do not let that become a claim either way.
6. **The server returns 5 samples, not 5 independent seeds.** Measured spread: 0.0141 mean
   within-job ipTM sd for fuzzy positives, 0.0146 for mutual folding — the plan expected
   fuzzy pairs to be more variable and they are not, because all five share a seed and a
   trunk embedding. Call it sample spread.
7. **The 464 ΔAUC intervals are uncorrected for multiplicity.** One non-control cell clears
   zero, which is roughly the chance expectation. Treat it as such.
8. **No PTMs.** FMRP/CAPRIN1 is explicitly phospho-dependent and so is tau LLPS. Two
   positives may be in the wrong functional state.
9. **No RNA.** G3BP1, CAPRIN1 and FMRP condensates are RNA-dependent. AF3 accepts RNA
   chains, so this is a scope choice and a natural extension.
10. **No monomer baseline.** Running each region alone (~30 jobs) would give a
    pLDDT/disorder reference to subtract. More attractive now that mean pLDDT turns out to
    be length-confounded (ρ = −0.60) and inverted in the headline test.
11. **hnRNPA2 boundary** (P22626, 194–353) is `approximate` — the literature construct is
    numbered on the A2 isoform, not B1.
12. **Species are mixed**: NCBD is mouse, GCN4 and Nsp1 yeast, ERD10 plant.

## 8. Practical notes

- **Budget**: 60 jobs/day across two accounts. §4 is 20 jobs, so one day on one account.
- **Seeds**: the batch files ship `"modelSeeds": []`, but 59 of the 60 jobs recorded seed
  `"1"` and N24 recorded a server-assigned `"1054611336"`. The `job_request.json` in each
  result folder is the record of what ran — the batch files are not.
- **Homodimers** are encoded as one `proteinChain` with `"count": 2`, not two chains.
- **Templates**: the switch is per chain, `"useStructureTemplate": false`, and it does not
  appear in a v1 batch file at all. The server rewrites the job to v3 and defaults it on.
- **The parser matches results to rows by job name**, normalised to alphanumerics, so the
  server's punctuation rewriting is handled. It prints anything it could not match — check
  that list before trusting the AUCs. Currently: 2 unmatched, both batch 4.
- **Rebuilding the report**: `python3 parse_af3_results.py data af3_idr_pairs_1.csv` then
  `python3 af3_report_stats.py` then `python3 af3_report_build.py`. The stats step takes
  about a minute on 16 cores; the parser holds all the full-data JSON in memory and peaks
  around 2.7 GB. Both report scripts have the repo path hard-coded — parameterise before
  moving them.
- **numpy** is required by the parser and the baseline; **reportlab** by the report builder,
  which also registers Arial from `/System/Library/Fonts/Supplemental` and so is macOS-
  specific as written.
- The sandbox where the set was built **cannot reach UniProt, MobiDB, FuzDB, MFIB or DIBS
  by curl** — only through a fetch tool that relays text through a model. Any future
  sequence work is safer done on a normal machine, where the FMR1 check was finally done.

## 9. Who to tell

- **The Mehdiabadi group** — this extends their benchmark into the regime they excluded,
  and the answer in that regime is a clean negative with a quantified trivial baseline.
- **Fuxreiter's group** for the fuzzy/mutual-folding labels (also the FuzPred contact from
  the other project; Vendruscolo, a FuzPred co-author, is in Cambridge Chemistry).
- **The Forman-Kay lab** for FMRP/CAPRIN1 and condensate IDR–IDR generally.
