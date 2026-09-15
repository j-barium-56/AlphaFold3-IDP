# Can AF3 co-folding detect IDR–IDR interactions?

Dataset and staged plan. Written 14 September 2026, revised the same day after a
leakage check changed the primary criterion, and again that evening with the stage-1
result. Everything marked *verified* was checked by running code or against a source.
Where I have not checked something, I say so.

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

Two follow-ups closed that out, both from a machine with real network access. FMR1
445–632, the one sequence resting on a three-fetch consensus rather than a mass check, was
compared character-for-character against UniProt Q06787: **exact**, and it is what P04 and
S04 actually ran. And the three homotypic negatives added late — p53 TAD, p21, 4E-BP1 —
turned out never to have been recorded in `af3_idr_regions.csv` at all, so they had no
accession or boundaries on file and were never in the 26/26 tally. They are now, checked
the same way: p53 TAD is P04637 residues 1–61 exactly, p21 is P38936 full-length (average
mass 18,119.3 Da, matching), 4E-BP1 is Q13541 full-length (12,580.0 Da, matching).
**30/30 regions verified and all 30 documented.**

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

**Amended after stage 1: that rule cannot be applied as written.** Those 8 negatives are
also the shortest and most helix-competent rows in the set, so the subset is length- and
foldability-confounded — length alone separates it at AUC 0.985. The disagreement is real
but it is not evidence about AF3. See Stage 1 below, and settle the length handling before
test 3 is allowed to adjudicate anything.

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

## Stage 1 — how it was run

`af3_batch_01..04.json` upload straight into AlphaFold Server (20 jobs each; homodimers
are one chain with `count: 2`). All job names are sanitised to `[A-Za-z0-9_-]` — the
server rejects anything else, which is what broke the first job of batch 1 in the first
version (the dot in `H1.0`).

As actually run: batches 1–3 went in one sitting — 60 jobs, which includes 8 of the 10
scrambles rather than holding the whole scramble arm for day 2 — and batch 4 (S20, S21)
is still outstanding.

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

## Stage 1 — the result (*verified*, 14 September 2026)

60 of 62 jobs ran (batches 1–3; S20 and S21 still pending), five samples each. Every
number below is in `af3_idr_results_report.pdf`; every cell is in `af3_all_statistics.csv`.

| Test | AF3 ipTM AUC | composition | ΔAUC (95% CI) | |
|---|---|---|---|---|
| 1. positives vs 31 real negatives | 0.611 | 0.767 | −0.156 [−0.32, +0.01] | fail |
| 2a. mutual folding vs hard zipper | 0.417 | 0.250 | +0.167 [−0.50, +0.88] | fail |
| 2b. fuzzy positives — **the result** | 0.676 | 0.784 | −0.107 [−0.26, +0.05] | fail |
| 3. vs gold-standard negatives | 0.118 | 0.640 | −0.522 [−0.74, −0.29] | fail |
| Homotypic arm, arm baseline | 0.564 | 0.902 | −0.337 [−0.58, −0.12] | fail |
| Heterotypic arm, arm baseline | 0.485 | 0.500 | −0.015 [−0.32, +0.29] | fail |

**Which branch.** ΔAUC_fuzzy ≈ 0 *and* AUC_fold low, which by the table above means
"debug before concluding". The debugging is done, and the pipeline is not what is broken:
GCN4-p1 reaches ipTM 0.73, all four memorised controls sit at 0.70–0.80 (the top of the
whole set), and all 60 jobs ran the sequence they were supposed to, five samples each.
What fails test 2a is the model: the three `hard_coiledcoil` negatives score 0.70, 0.78
and 0.80 — **ipTM AUC 0.000** against the positives, meaning every one of them outranks
every positive. AF3 dimerises any zipper-competent pair, cognate or not.

**Test 3 cannot adjudicate, so criterion 0 cannot be applied as written.** The 8
`demonstrated` negatives are 94–236 residues; the 17 positives are 205–882. Total chain
length alone separates that exact split at AUC 0.985, and the metrics carrying the
inversion are the four the parser flags as length-confounded. Residualising on log length
moves ipTM from 0.118 to 0.449. Decide the length handling — matched subset or partialled —
and write it down *before* reading test 3 again.

**The secondary panel does not rescue it.** Interchain contact-probability sum is the one
metric that separates on its own (AUC 0.787, q = 0.008; 0.970 in the homotypic arm), and
it still does not beat composition: Δ = +0.020 [−0.10, +0.14]. Of the 464 cells with a
baseline, 15 have a Δ interval clear of zero and 14 of those are the memorised controls.
The fifteenth — fuzzy homotypic positives on contact-probability sum, Δ +0.091 [+0.011,
+0.204] — is uncorrected across 464 comparisons, and the pre-registered homotypic arm
does not pass. That is a hypothesis for a fresh set, not a result.

**The homotypic arm's failure has a structural explanation, not just a statistical one:
AF3 stacks disordered chains indiscriminately rather than recognising a specific
interface.** `contact_participation_ratio` (inverse participation ratio of the inter-chain
contact-probability mass — near 1 means one dominant contact patch, large means the mass
is spread thinly over many weak ones) is already in the 16-metric panel. On the two
`homo_fold` positives, whose complexes are real solved structures, it sits at 52–108: a
small, localised interface, exactly what a specific binding site looks like. On the 14
`homo_fuzzy` positives (FUS, hnRNPA1/2, TDP-43, EWSR1, TAF15, DDX4, tau, α-synuclein,
UBQLN2, Nsp1-FG, CAPRIN1) it sits at 700–2900 — contacts smeared across the whole chain.
So do the `homo_negative` controls that never appear in any self-association literature
(PTMA, H1.0, NPM1, ERD10, G3BP1, p53 TAD, p21, 4E-BP1): 1100–2800, statistically the same
range. The pre-registered homotypic test on this metric bears that out directly — AUC
0.511 [0.291, 0.736], chance — while composition still separates the same split at 0.896.
AF3 does not fail to detect fuzzy self-association *weakly*; it produces the same diffuse,
low-confidence, spread-out contact pattern for a real LLPS-forming LC domain and for an
acidic IDR nobody has ever reported self-associating. That pattern is consistent with a
generic packing prior for aromatic/low-complexity sequences — the same property that makes
composition alone predictive — rather than with sequence-specific recognition of either
kind. Any visual read of these structures as "the two chains are clearly interacting"
should be treated as this artefact, not as evidence of binding, unless the contact map is
also localised.

**What did work: the matched-internal contrasts.** They rest on no literature absence at
all, and AF3 gets all three right — ACTR×NCBD above both ACTR×ACTR and NCBD×NCBD, and
c-Fos×c-Jun above c-Fos×c-Fos; 3/3 on ipTM, 12 of 16 metrics for P06 and 15 of 16 for
P07. This is the design to expand. With the caveat immediately below.

**Templates were on, and nobody chose that.** The batch files are dialect version 1 and
never set `useStructureTemplate`; the server upgraded them to version 3 and set it `true`
on all 95 chain entries, returning 2–4 template hits per chain for all 60 jobs. 6es7 —
which contains both ACTR and NCBD, since it templates both chains of P06 — is also a
template for N26 and N27, and 1fos is a template in N15, N16, N17 and N28. AF3 applies
templates within a chain, so this is each partner's bound conformation rather than the
interface itself, but it is exactly the help a zipper needs and it lands on both the
memorised positives and the matched-internal negatives. Until the templates-off rerun,
every memorisation statement in this document is provisional.

## Stage 2 — what to run next

Revised after stage 1. Twenty jobs, one day on one account, in this order:

1. **Batch 4** (S20, S21) — finishes the scramble arm at 10 pairs. 2 jobs. Still pending.
2. **Templates off — done, run 15 Sep, N28 not yet submitted.** P06, P07, P20, P21, N15,
   N16, N17, N26, N27 re-ran with `"useStructureTemplate": false` on every chain
   (`af3_batch_05_notmpl.json`, rows T01–T09 in `af3_idr_pairs_1.csv`, stage 2). Result —
   ipTM with templates vs without:

   | Job | templated | no-template | Δ |
   |---|---|---|---|
   | P06 ACTR×NCBD (memorised) | 0.70 | 0.73 | +0.03 |
   | P07 Fos×Jun (memorised) | 0.74 | 0.72 | −0.02 |
   | P20 MAX×MAX (memorised) | 0.80 | 0.77 | −0.03 |
   | P21 GCN4×GCN4, canary | 0.73 | 0.68 | −0.05 |
   | N15 Fos×MAX (hard zipper) | 0.78 | 0.77 | −0.01 |
   | N16 Jun×MAX (hard zipper) | 0.70 | 0.74 | +0.04 |
   | N17 GCN4×Jun (hard zipper) | 0.80 | 0.77 | −0.03 |
   | N26 ACTR×ACTR (matched-internal neg.) | 0.52 | 0.48 | −0.04 |
   | N27 NCBD×NCBD (matched-internal neg.) | 0.46 | 0.18 | −0.28 |

   Eight of nine move by ≤0.05 either direction. **The memorised-control signal and the
   hard-zipper promiscuity are not template artefacts** — pulling the template does not
   collapse P06/P07/P20/P21, and it does not collapse the three hard-zipper negatives
   either, so test 2a's failure ("AF3 pairs any two zippers") is a property of the model's
   coiled-coil prior, not of template retrieval. The one large mover, N27, moves the
   *right* direction for a true negative — 0.46 to 0.18 — which if anything strengthens
   the matched-internal design rather than undermining it. Caveat #13 below is resolved by
   this rerun; every memorisation statement elsewhere in this document should be read as
   confirmed, not provisional. N28 (FOS×FOS homodimer negative) was never submitted
   without templates and is the one loose end here.
3. **Seed replicates.** `"modelSeeds": ["2","3"]` on the four memorised controls. Seed 1
   already ran on 59 of the 60 jobs, so only 2 and 3 add information. 8 jobs, and the
   first real measurement of run-to-run variance — the five samples per job share a seed.
4. **Confirm or drop** the two `uncertain` negatives (H1.0 alone, ERD10 alone) before they
   carry any weight in a second pass.

Then the branch. Fuzzy discrimination does not exist at this sample size, so the
informative direction is where the boundary sits: truncate the GCN4 zipper toward
marginality and locate the length at which ipTM collapses, and expand the
matched-internal design, which is the only contrast in the set that behaved.

## Stage 3 — write-up

Short paper whose contribution is the labelled IDR–IDR set plus one plane
(interaction class × seed variance, or × ΔAUC). Natural readers: the Mehdiabadi group,
Fuxreiter (fuzzy labels), Forman-Kay (FMRP/CAPRIN1 and condensate IDR–IDR generally).

Stage 1 reframes it slightly. The seed-variance plane is not available yet, so the paper
that exists today is the labelled set plus **a quantified negative against a trivial
baseline**: co-folding confidence adds nothing over sequence composition for fuzzy
IDR–IDR pairs, it cannot separate cognate from non-cognate zippers at all, and the
gold-standard subset that was supposed to arbitrate is confounded by length. The
methodological point — that an IDR benchmark has to beat composition, not chance — is
probably as useful to the field as the AF3 result itself.

---

## What this design still does not cover

Listed so none of it arrives as a surprise from a referee.

1. **A dimer is not a condensate.** Phase separation is a many-body, concentration-
   dependent phenomenon; a 1:1 co-fold cannot in principle report on it. The `llps`
   column is descriptive metadata for stratifying, **never an outcome variable**.
2. **The server gives 5 samples, not 5 independent seeds**, unless `modelSeeds` is set
   explicitly. "Inter-seed spread" is really inter-sample spread — fine as a variance
   readout, but do not call it seed replication in the write-up. As run: 59 of the 60 jobs
   recorded seed `"1"` and N24 recorded a server-assigned seed, so the measured spread
   (0.0141 mean within-job ipTM sd for fuzzy positives, 0.0146 for mutual folding) is
   sample spread across one seed. Run-to-run variance is still unmeasured.
3. **No monomer baseline.** Running each region alone would give a pLDDT/disorder
   reference to subtract. Costs ~30 jobs; worth it if the pLDDT metrics turn out to matter.
4. **No PTMs.** FMRP/CAPRIN1 is explicitly phospho-dependent, and so is tau LLPS. The
   unmodified sequences may be the wrong functional state for at least two positives.
5. **No RNA.** G3BP1, CAPRIN1 and FMRP condensates are RNA-dependent. AF3 accepts RNA
   chains, so this is a scope choice, not a limitation of the tool.
6. **Length is not matched by design**, only checked post hoc — and stage 1 shows it
   decides answers: length alone scores AUC 0.985 on the gold-standard split, 1.000
   against the hard zippers and 0.676 on the headline test. A second set should be
   length-matched at selection time rather than corrected afterwards.
7. **Cross-kingdom negatives have no paired MSA**, so they share the scramble arm's
   confound in milder form. The non-cognate human–human pairs are the cleanest negatives;
   report negative subtypes separately if they behave differently.
8. **"No reported interaction" is not "does not interact."** IDRs are promiscuous, and
   21 of the 31 negatives rest on `inferred` or `absence_only` evidence. This is the
   single biggest threat to the benchmark's validity. The `evidence_strength` column
   exists so a reader can discount accordingly, and the fix if the result hinges on it is
   to expand the `matched_internal` design — more chains that appear in both a positive
   and a negative — rather than to find more proteins nobody has reported binding.
9. **No replicate jobs** to measure the server's own run-to-run variance — scheduled as
   stage 2 step 3, and until it runs no variance claim in the write-up is supportable.
10. **Species are mixed**: NCBD is mouse, GCN4 and Nsp1 are yeast, ERD10 is plant.
11. **FMRP sequence caveat — resolved.** FMR1 445–632 was the one sequence not confirmed
    by mass, because the fetch route corrupted the protein's N-terminal half. It has since
    been compared character-for-character against UniProt Q06787 and is exact; P04 and S04
    ran the correct sequence.
12. **hnRNPA2 boundary** (P22626, 194–353) is marked `approximate`: the literature
    construct is numbered on the A2 isoform, not B1.
13. **Structure templates were on for all 95 chain entries** — the AlphaFold Server
    default, never a choice made here. **Resolved by the stage-2 templates-off rerun**
    (above): removing templates moves the four memorised controls and three hard-zipper
    negatives by ≤0.05 ipTM. The memorisation and zipper-promiscuity findings hold without
    templates and are no longer provisional.
14. **The ΔAUC intervals are not corrected for multiplicity.** The 16-metric panel carries
    Benjamini–Hochberg correction within each group, but the 464 Δ-vs-baseline intervals
    do not. Exactly one non-control cell clears zero, which is about what 464 uncorrected
    95% intervals would produce by chance.
15. **No length-matched negative arm**, which is the single change most likely to make
    test 3 interpretable on a second pass.

---

## Files

| File | What it is |
|---|---|
| `af3_idr_pairs_1.csv` | 62 jobs, one row each: both sequences, all labels, `evidence_strength`, stage. The `_1` is a download artefact; `parse_af3_results.py` defaults to `af3_idr_pairs.csv`, so pass the path explicitly or rename it. The copy in `~/Downloads` is an **older** version with no `evidence_strength` column — do not use it |
| `af3_idr_regions.csv` | 30 IDR regions: accession, boundaries, sequence, boundary confidence. p53 TAD, p21 and 4E-BP1 were back-filled after stage 1 — they were in the pair table but had never been recorded here. The histone id is `H1.0_full`; the pair table sanitises it to `H10_full` |
| `af3_batch_01..04.json` | Upload-ready AlphaFold Server batch files. Dialect version 1, `modelSeeds: []`, no template switch — none of which is what actually ran; the `job_request.json` inside each result folder is the record |
| `baseline_leakage.py` | Composition-only baseline — **run before AlphaFold** |
| `parse_af3_results.py` | Result zips → scored table, ΔAUC tests, corrected metric panel |
| `af3_report_stats.py` | Extends the parser to every grouping × metric → `af3_report_data.json` |
| `af3_report_build.py` | `af3_report_data.json` → `af3_idr_results_report.pdf` |
| `af3_idr_results_report.pdf` | Stage-1 results, 26 pages: verdicts, every grouping, confounds, both appendices |
| `af3_all_statistics.csv` | All 496 cells: AUC, CI, permutation p, BH q, ΔAUC and its CI |
| `af3_scored.csv` | Per-job scored table written by `parse_af3_results.py` |
| `af3_idr_plan.md` | This document |
| `af3_idr_handoff.md` | The operational companion — what to run next, and why not to undo things |
