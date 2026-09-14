# AF3 IDR–IDR benchmark — handoff

**Status:** dataset built and verified, nothing submitted to AlphaFold yet.
**Written:** 14 September 2026. **Companion doc:** `af3_idr_plan.md` (the scientific design).
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
| All sequences verified against UniProt length + mass | **done**, 26/26 exact |
| AlphaFold Server batch JSONs | **done**, 4 files, names validated |
| Composition-baseline leakage check | **done**, run it again if you change the set |
| Result parser + statistics | **done**, tested on synthetic output end to end |
| AlphaFold Server submissions | **not started** |

## 3. Files

| File | Role |
|---|---|
| `af3_idr_pairs.csv` | The dataset. One row per job: both sequences, all labels, `evidence_strength`, stage |
| `af3_idr_regions.csv` | 30 IDR regions: accession, boundaries, sequence, `boundary_confidence` |
| `af3_batch_01..04.json` | Upload-ready AlphaFold Server batch files (20/20/20/2) |
| `baseline_leakage.py` | Composition-only baseline. **Run before AlphaFold** |
| `parse_af3_results.py` | Result zips → scored table, ΔAUC tests, corrected metric panel |
| `af3_idr_plan.md` | Scientific design, pre-registered criteria, limitations |
| `af3_idr_handoff.md` | This file |

Build scripts (`build_dataset.py`, `verify.py`, `raw_fetch_1.fasta`) live in the session
workspace and are only needed if the set is rebuilt. Ask for them if so.

## 4. Do this next — day 1

1. **Run the leakage check.** `python3 baseline_leakage.py af3_idr_pairs.csv`.
   Write down the three LOO-logistic AUCs (currently 0.70 all / 0.53 heterotypic /
   0.85 homotypic). These are the numbers AF3 has to beat.
2. **Fix the FMRP sequence.** Copy UniProt **Q06787 residues 445–632** and paste it over
   `FMR1_CT` wherever it appears (rows P04 and its scramble). It is the only sequence in
   the set not confirmed by mass — see §6.
3. **Submit batches 1–3** (60 jobs) on account A, **batch 4** (2 jobs) on account B.
   Leave `"modelSeeds": []` as shipped.
4. Download every result zip into one folder as they finish. Do not rename them.
5. **Parse**: `python3 parse_af3_results.py <results_folder> af3_idr_pairs.csv`
   (add `--no-full-data` for a fast first look; the tau homodimer's full-data file is large).

## 5. Do this next — day 2

The scramble arm (`stage` = 2, 10 jobs) plus explicit-seed replicates of the four
memorised controls: `"modelSeeds": ["1","2","3"]` on P06, P07, P20, P21 only. That is
~22 jobs, leaving ~38 for whichever branch stage 1 opens.

Read the outcome table in `af3_idr_plan.md` before deciding the branch.

## 6. Decisions already made, and why — do not silently undo these

**Sequences are mass-verified, and one is not.** The only network route available in the
build environment relays text through a language model, which silently garbled the FMR1
sequence three separate times in three different places. Every sequence was therefore
checked against UniProt's reported length *and* average molecular mass. 26/26 pass.
FMR1's C-terminal region is a three-fetch consensus — replace it before it matters.
**If you add any sequence, verify it the same way.**

**Scrambles are not the primary negative.** A shuffled sequence has no homologs, so AF3
builds it no MSA; low confidence would then follow from missing evolutionary information
rather than a missing interface. The primary negatives are real IDRs in non-cognate
pairings. The scramble arm asks a different question and is reported separately.

**The comparator is the composition baseline, not 0.5.** Composition alone separates
these labels at AUC 0.85 in the homotypic arm — LC domains that self-associate genuinely
*are* aromatic-rich and weakly charged, which is biology and cannot be curated away. So
the primary test is ΔAUC = AF3 − baseline with a paired bootstrap CI excluding zero.
On simulated data an AF3 AUC of 0.84 against a baseline of 0.77 gives Δ = +0.07 with a
CI crossing zero: a **fail** that looks like a clean win if you compare against chance.

**Homotypic and heterotypic are analysed separately.** Different baselines (0.85 vs 0.53),
different MSA situations. Pooling lets chain identity leak into the answer.

**ipTM stays the pre-registered primary even though it may be useless.** The 15 secondary
metrics exist because ipTM is likely floored for fuzzy pairs, but choosing the best one
after seeing data is p-hacking. The panel is reported with permutation p-values and
Benjamini–Hochberg correction; a secondary metric passing at q < 0.05 is stated as exactly
that — corrected, hypothesis-generating, needs a fresh set to confirm.

**Four positives are memorised** (P06/1KBH, P07/1FOS, P20/1AN2, P21/2ZTA). Pipeline
controls, excluded from the headline. **GCN4-p1 (P21) is the canary**: if AF3 cannot
confidently dimerise it, stop and debug rather than concluding anything about IDRs.

**Job names are sanitised to `[A-Za-z0-9_-]`.** The server rejects anything else — a dot
in `H1.0` broke the first job of batch 1 in the first version. If you add jobs, sanitise.

## 7. Open questions and known weaknesses

Ordered by how likely they are to matter.

1. **20 of 31 negatives rest on `inferred` or `absence_only` evidence.** "No paper reports
   this interaction" is not "these do not interact", and IDRs are promiscuous. The
   `evidence_strength` column grades every row; test 3 in the parser uses only the 8
   `demonstrated` negatives. **If the headline and gold-standard results disagree, believe
   the gold-standard one.** The structural fix is to expand the `matched_internal`
   design — more chains appearing in both a positive and a negative, like ACTR, NCBD and
   c-Fos — not to find more proteins nobody has reported binding.
2. **Two negatives are flagged `uncertain`**: H1.0 alone (H1 coacervates with
   single-stranded nucleic acid; histone LLPS reported in cells) and ERD10 alone (plant
   condensate biology is active; dehydrin condensation under dehydration not ruled out).
   Confirm or drop both.
3. **A dimer is not a condensate.** LLPS is many-body and concentration-dependent; a 1:1
   co-fold cannot in principle report on it. The `llps` column is stratifying metadata,
   **never an outcome variable**. Do not let this slip in the write-up.
4. **The server returns 5 samples, not 5 independent seeds**, unless `modelSeeds` is set.
   Call it sample spread, not seed replication.
5. **No PTMs.** FMRP/CAPRIN1 is explicitly phospho-dependent and so is tau LLPS. Two
   positives may be in the wrong functional state.
6. **No RNA.** G3BP1, CAPRIN1 and FMRP condensates are RNA-dependent. AF3 accepts RNA
   chains, so this is a scope choice and a natural extension.
7. **No monomer baseline.** Running each region alone (~30 jobs) would give a
   pLDDT/disorder reference to subtract. Worth it if the pLDDT metrics turn out to matter.
8. **Length is not matched by design**, only checked post hoc (the parser flags any metric
   with |ρ| > 0.5 against total residue count).
9. **hnRNPA2 boundary** (P22626, 194–353) is `approximate` — the literature construct is
   numbered on the A2 isoform, not B1.
10. **Species are mixed**: NCBD is mouse, GCN4 and Nsp1 yeast, ERD10 plant.

## 8. Practical notes

- **Budget**: 60 jobs/day across two accounts. Stage 1 is 52 real-sequence jobs, so one day.
- **Seeds**: `"modelSeeds": []` for the main run; `["1"]` if you want it fixed — a **list
  of strings** of uint32 values, not integers. Same seed across all jobs so they compare.
- **Homodimers** are encoded as one `proteinChain` with `"count": 2`, not two chains.
- **The parser matches results to rows by job name**, normalised to alphanumerics, so the
  server's punctuation rewriting is handled. It prints anything it could not match — check
  that list before trusting the AUCs.
- **numpy is required** by the parser; `baseline_leakage.py` needs it too.
- The sandbox where this was built **cannot reach UniProt, MobiDB, FuzDB, MFIB or DIBS by
  curl** — only through a fetch tool that relays text through a model. Any future sequence
  work is safer done on a normal machine.

## 9. Who to tell

- **The Mehdiabadi group** once stage 1 has a result either way — this extends their
  benchmark into the regime they excluded.
- **Fuxreiter's group** for the fuzzy/mutual-folding labels (also the FuzPred contact from
  the other project; Vendruscolo, a FuzPred co-author, is in Cambridge Chemistry).
- **The Forman-Kay lab** for FMRP/CAPRIN1 and condensate IDR–IDR generally.
