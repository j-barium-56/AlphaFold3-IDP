# Can AF3 co-folding detect IDR–IDR interactions?

Dataset and staged plan. Written 14 September 2026, revised the same day after a
leakage check changed the primary criterion. Everything marked *verified* was checked by
running code or against a source. Where I have not checked something, I say so.

---

## The question

> Given two disordered regions and nothing else, does AlphaFold 3 co-folding
> distinguish pairs that really interact from pairs that do not?

Not "does AF3 get the structure right" — for most of these pairs there is no structure to
get right. The test is **discrimination**, scored on the confidence outputs AF3 already
produces.

## Why this is still open — *verified*

The closest paper is **Mehdiabadi et al., "AlphaFold3 and Intrinsically Disordered
Proteins: Reliable Monomer Prediction, Unpredictable Multimer Performance"**
(bioRxiv, Dec 2025; Protein Science 2026). It benchmarks AF3 vs AF2 on 90 CAPRI
**IDR–receptor** dimers (one disordered chain, one folded chain) and on **MFIB 2.0**
fold-upon-binding complexes. It does **not** cover IDR–IDR pairs where both chains are
disordered, disorder-to-disorder (fuzzy) complexes, homotypic self-association, or phase
separation, and it states it cannot "prospectively identify which disordered regions will
undergo binding-induced folding." Its headline numbers use DockQ, which needs a reference
structure — a metric that does not exist for most pairs below.

## The trap

The failure mode from the murmuration project was showing something the tool's authors
already concede. "AF3 is bad at disordered complexes" is that result. The version that
adds something is **stratified**, with the confounds separated out rather than left to
explain the answer.

---

## Stage 0 — the dataset

**62 jobs. 21 positives, 31 real-sequence negatives, 10 composition-matched scrambles.**

Every sequence was pulled from UniProt and *verified* against UniProt's own reported
length **and** average molecular mass — 26/26 exact. This matters: the network route
available here relays text through a language model, and it silently garbled one sequence
(FMR1) three separate times, in three different places. Mass-checking caught it.

### Positives (21)

| Class | n | Examples |
|---|---|---|
| Heterotypic, fuzzy | 5 | ProTα×H1.0; G3BP1 IDR1×RGG; FMRP×CAPRIN1; SRSF1 RS×U1-70K RS; α-syn×tau |
| Heterotypic, mutual folding | 2 | ACTR×NCBD (1KBH); c-Fos×c-Jun (1FOS) |
| Homotypic, fuzzy, LLPS | 12 | FUS LC, hnRNPA1 LC, hnRNPA2 LC, TDP-43 CTD, DDX4 N-IDR, EWSR1 LC, TAF15 LC, UBQLN2 C-term, Nsp1 FG, tau, α-syn, CAPRIN1 C-term |
| Homotypic, mutual folding | 2 | Max bHLHZip (1AN2); GCN4 leucine zipper (2ZTA) |

### Negatives (31 real-sequence + 10 scrambles)

| Type | n | Logic |
|---|---|---|
| `noncognate` | 10 | Two real IDRs, no reported interaction |
| `cross_kingdom` | 5 | Human IDR × plant dehydrin ERD10, or yeast region |
| `charge_repulsive` | 2 | Two strongly acidic IDRs |
| `hard_coiledcoil` | 3 | Fos×Max, Jun×Max, GCN4×Jun — zipper-competent, do **not** pair |
| `homo_negative` | 11 | IDRs that do **not** self-associate, run against themselves |
| `scramble` | 10 | Composition-matched shuffle of chain B — **separate arm** |

### How do we know the negatives are negative?

This is the weakest joint in any benchmark of this kind, so it is graded per row in an
`evidence_strength` column rather than assumed uniform:

| Grade | n | Meaning |
|---|---|---|
| `demonstrated` | 5 | The monomeric / non-interacting state is explicitly characterised |
| `demonstrated+matched_internal` | 3 | Same chain also appears in a positive, so the contrast is internal |
| `inferred` | 11 | Physical argument (charge repulsion, compartment) plus no reports |
| `absence_only` | 10 | Nothing more than "no paper reports this interaction" |
| `uncertain` | 2 | Literature exists that could contradict the label |

The three `matched_internal` rows are the sharpest tests in the set, because they do not
rest on a literature absence at all: **ACTR alone** and **NCBD alone** (each folds only
with the other) against P06, and **c-Fos alone** (no stable homodimer; O'Shea et al.,
*Science* 1989) against P07. Same chains, same MSAs, only the pairing differs. If AF3
scores ACTR×NCBD above ACTR×ACTR, that is unambiguous.

The two `uncertain` rows are flagged honestly: **H1.0 alone** (H1 coacervates with
single-stranded nucleic acid and histone LLPS is reported in cells — it should not
self-associate alone, but this was not confirmed) and **ERD10 alone** (plant condensate
biology is active and dehydrin condensation under dehydration stress is not ruled out).
Check both before they carry any weight, or drop them.

Test 3 in the script runs positives against the 8 `demonstrated` negatives only. It has
less power, but its labels are defensible; if the headline result and the gold-standard
result disagree, believe the gold-standard one.

### Two design decisions that came out of checking the set against itself

**Scrambles cannot be the primary negative.** A shuffled sequence has no homologs, so AF3
builds it no MSA. Low confidence then follows from missing evolutionary information, not
from a missing interface. The primary negatives are therefore real IDRs in non-cognate
pairings — both chains have genuine MSAs, only the pairing is wrong. The scramble arm is
kept, asking a different question: does the signal survive loss of evolutionary information?

**The comparator is not 0.5.** `baseline_leakage.py` fits a leave-one-out logistic model
on cheap composition features (NCPR, FCR, aromatic fraction, G/S content, hydropathy,
length, charge complementarity, homotypic flag) and asks how well *those alone* separate
the labels. Result, *verified by running it*:

| Subset | best single feature | LOO-logistic AUC |
|---|---|---|
| All real-sequence pairs | `polar_qn_absdiff` 0.72 | **0.70** |
| Heterotypic only | `polar_qn_absdiff` 0.79 | 0.53 |
| Homotypic only | `aromatic_mean` 0.86 | **0.85** |

The homotypic arm is the problem. LC domains that self-associate genuinely *are*
aromatic-rich and weakly charged, so composition predicts the label at AUC 0.85 with no
structural information whatsoever. I added three aromatic-containing monomeric IDRs
(p53 TAD, p21, 4E-BP1) as homotypic negatives, which pulled `charge_product` from 0.97
to 0.85 — but this is biology, not an imbalance that can be curated away.

So the pre-registered criterion changes. **AF3 must beat the composition baseline, not
chance.**

---

## Stage 1 — run it (day 1, 52 real-sequence jobs; scrambles on day 2)

`af3_batch_01..04.json` upload straight into AlphaFold Server (20 jobs each; homodimers
are one chain with `count: 2`). All job names are sanitised to `[A-Za-z0-9_-]` — the
server rejects anything else, which is what broke the first job of batch 1 in the first
version (the dot in `H1.0`).

`parse_af3_results.py` reads the result zips and writes the scored table plus every test
below.

### Metrics, fixed before looking at the data

**Primary**: ipTM of the top-ranked model.

**Secondary panel (16 metrics)**, because ipTM may be floored at ~0 for every disordered
pair and carry no information at all:

- from `summary_confidences.json`: `chain_pair_iptm`, `chain_pair_pae_min`, `ptm`,
  `ranking_score`, `fraction_disordered`, `has_clash`
- from `full_data.json`, sliced to the **interchain block**: mean / 5th-percentile PAE;
  **fraction of interchain token pairs below 10 Å and 5 Å PAE** (length-normalised, unlike
  the minima); interchain contact-probability sum, density, max and count above 0.5;
  interface pLDDT
- `contact_participation_ratio` — participation ratio of the interchain contact mass.
  Near 1 means one dominant contact; large means mass smeared over many weak ones. This
  is idea D from the mechanism doc ("is *many weak interactions* actually true?") answered
  directly off AF3's output, for free.

The panel is 16 wide, so choosing the best-separating metric post hoc is p-hacking. The
script reports the panel with permutation p-values and Benjamini–Hochberg correction. A
secondary metric that passes at q < 0.05 is a legitimate finding stated as exactly that:
corrected, hypothesis-generating, worth confirming on a fresh set.

### Pre-registered criteria

0. **Report the gold-standard subset alongside the headline.** Test 3 uses only
   negatives graded `demonstrated`; the rest of the negatives are softer evidence.
1. **Primary**: ΔAUC = AUC(AF3 ipTM) − AUC(composition baseline), on positives vs real
   negatives, paired bootstrap 95% CI **excluding zero**. Beating 0.5 is not the bar.
2. **Stratified — this is the actual question**:
   - **AUC_fold**: mutual-folding positives (P06, P07, P20, P21) vs the three
     `hard_coiledcoil` negatives. Expect high. If it is not, the pipeline is broken, not AF3.
   - **AUC_fuzzy**: the 17 fuzzy positives vs the remaining real negatives, again as ΔAUC
     over baseline. **This is the result.**
3. **Analyse homotypic and heterotypic separately.** They have different baselines (0.85
   vs 0.53) and different MSA situations, and pooling them lets chain identity leak into
   the answer.
4. **Length confound**: any metric with |Spearman ρ| > 0.5 against total residue count is
   flagged and reported partialled as well as raw.
5. **Memorisation**: P06 (1KBH), P07 (1FOS), P20 (1AN2), P21 (2ZTA) predate AF3 training.
   Pipeline controls, **not** evidence about prediction; excluded from the headline number.
   GCN4-p1 is the canary — if AF3 cannot confidently dimerise it, stop and debug.

### What each outcome means

- **ΔAUC_fuzzy > 0 with CI excluding zero** — the interesting result. AF3 carries
  information about interactions with no folded interface, beyond what composition gives.
  Next question: what is it reading? The scramble arm starts answering that.
- **ΔAUC_fuzzy ≈ 0, AUC_fold high** — the clean expected split, and still publishable: a
  quantified boundary showing co-folding confidence is informative down to mutual folding
  and adds nothing over composition for fuzzy binding. That is the benchmark the field
  does not have.
- **Both low** — debug before concluding.

## Stage 2 — day 2

Scramble arm (10 jobs) plus explicit-seed replicates of the four memorised controls.
That leaves ~45 jobs for whichever branch Stage 1 opens: if fuzzy discrimination exists,
test what drives it with charge-pattern-preserving vs full scrambles; if it does not,
find where the boundary sits by truncating the GCN4 zipper toward marginality and
locating the length at which ipTM collapses.

## Stage 3 — write-up

Short paper whose contribution is the labelled IDR–IDR set plus one plane
(interaction class × seed variance, or × ΔAUC). Natural readers: the Mehdiabadi group,
Fuxreiter (fuzzy labels), Forman-Kay (FMRP/CAPRIN1 and condensate IDR–IDR generally).

---

## What this design still does not cover

Listed so none of it arrives as a surprise from a referee.

1. **A dimer is not a condensate.** Phase separation is a many-body, concentration-
   dependent phenomenon; a 1:1 co-fold cannot in principle report on it. The `llps`
   column is descriptive metadata for stratifying, **never an outcome variable**.
2. **The server gives 5 samples, not 5 independent seeds**, unless `modelSeeds` is set
   explicitly. "Inter-seed spread" is really inter-sample spread — fine as a variance
   readout, but do not call it seed replication in the write-up.
3. **No monomer baseline.** Running each region alone would give a pLDDT/disorder
   reference to subtract. Costs ~30 jobs; worth it if the pLDDT metrics turn out to matter.
4. **No PTMs.** FMRP/CAPRIN1 is explicitly phospho-dependent, and so is tau LLPS. The
   unmodified sequences may be the wrong functional state for at least two positives.
5. **No RNA.** G3BP1, CAPRIN1 and FMRP condensates are RNA-dependent. AF3 accepts RNA
   chains, so this is a scope choice, not a limitation of the tool.
6. **Length is not matched by design**, only checked post hoc.
7. **Cross-kingdom negatives have no paired MSA**, so they share the scramble arm's
   confound in milder form. The non-cognate human–human pairs are the cleanest negatives;
   report negative subtypes separately if they behave differently.
8. **"No reported interaction" is not "does not interact."** IDRs are promiscuous, and
   20 of the 31 negatives rest on `inferred` or `absence_only` evidence. This is the
   single biggest threat to the benchmark's validity. The `evidence_strength` column
   exists so a reader can discount accordingly, and the fix if the result hinges on it is
   to expand the `matched_internal` design — more chains that appear in both a positive
   and a negative — rather than to find more proteins nobody has reported binding.
9. **No replicate jobs** to measure the server's own run-to-run variance.
10. **Species are mixed**: NCBD is mouse, GCN4 and Nsp1 are yeast, ERD10 is plant.
11. **FMRP sequence caveat**: FMR1 445–632 is the one sequence not confirmed by mass —
    the fetch route corrupted the protein's N-terminal half. The C-terminal region agreed
    across three independent fetches, but re-copy it from UniProt Q06787 to be clean.
12. **hnRNPA2 boundary** (P22626, 194–353) is marked `approximate`: the literature
    construct is numbered on the A2 isoform, not B1.

---

## Files

| File | What it is |
|---|---|
| `af3_idr_pairs.csv` | 62 jobs, one row each: both sequences, all labels, `evidence_strength`, stage |
| `af3_idr_regions.csv` | 30 IDR regions: accession, boundaries, sequence, boundary confidence |
| `af3_batch_01..04.json` | Upload-ready AlphaFold Server batch files |
| `baseline_leakage.py` | Composition-only baseline — **run before AlphaFold** |
| `parse_af3_results.py` | Result zips → scored table, ΔAUC tests, corrected metric panel |
| `af3_idr_plan.md` | This document |
