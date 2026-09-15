# AF3 IDR–IDR benchmark — handoff

**Status:** stage 1 scored. 60 of 62 jobs run, parsed and reported; batch 4 (2 jobs)
pending. Stage 2 templates-off ablation is now fully closed (10 of 10 jobs — N28's
rerun landed 15 Sep evening as job T11). **Stage 2b landed 15 Sep evening: 13 new
jobs** expanding the matched-internal design into non-fuzzy SLiM/domain pairs, plus a
corrected rerun of the corrupted P22/T10 full-length p53–MDM2 pair. See §5b below.
**Written:** 14 September 2026, updated 15 September (stage-1 result, templates-off /
contact-localisation follow-ups, then again the same evening with stage 2b).
**Companion doc:** `af3_idr_plan.md` (the scientific design and the numbers in full;
not yet updated with stage 2b — this file is currently the more current record for it).
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
| AlphaFold Server submissions | **73 of 75** — batches 1–3, 5, 7(stage2b) done; batch 4 (S20, S21) still pending |
| Result parser + statistics | **done**, rerun 15 Sep evening over 85 scored jobs |
| Results report | **done**, `af3_idr_results_report.pdf` (30 pages), site rebuilt (`docs/`) |
| Stage 2 templates-off ablation | **done, 10/10** — N28/T11 closed the last gap, see §5b |
| Stage 2b (matched-internal SLiM/domain expansion) | **done, scored** — see §5b. Not yet folded into `af3_idr_plan.md`'s pre-registered tests, which still pool everything into the original 31 groupings |

## 3. Files

| File | Role |
|---|---|
| `af3_idr_pairs_1.csv` | The dataset. One row per job: both sequences, all labels, `evidence_strength`, stage |
| `af3_idr_regions.csv` | 30 IDR regions: accession, boundaries, sequence, `boundary_confidence`. p53 TAD, p21 and 4E-BP1 were back-filled after stage 1; the histone id is `H1.0_full` here and `H10_full` in the pair table |
| `af3_batch_01..04.json` | Upload-ready AlphaFold Server batch files (20/20/20/2) |
| `af3_batch_05_notmpl.json` | Templates-off ablation, 10 jobs — **all 10 run** (9 on 15 Sep, N28 followed the same evening as job T11) |
| `af3_batch_06_monomers.json` | **New, not yet submitted.** 30 single-chain jobs, one per region in `af3_idr_regions.csv`, names `M_<region_id>` sanitised. `parse_af3_results.py` won't pick these up as-is — it expects a pairs CSV row per job; a monomer needs its own small script reading `summary_metrics()` only (no `interchain_metrics`, since there's one chain) — write that when the results come back rather than before |
| `af3_batch_07_stage2b.json` | **Run 15 Sep evening, scored.** 13 jobs: corrected p53–MDM2 full-length rerun (P23), the minimal 1YCR interface (P24/N32), and 4 new SLiM-class matched-internal pairs (LXXLL/NR, SH3/PxxP, 14-3-3/client, KEN-box/APC-C) each with a motif-dead mutant or scramble negative, plus 2 cross-motif promiscuity checks. `useStructureTemplate: false` on every chain throughout, deliberately — several of these are real solved structures (1YCR, 1GWQ, 1ABO, 1A37, 4GGD) that would otherwise trivially retrieve their own answer as a template. See §5b |
| `baseline_leakage.py` | Composition-only baseline. **Run before AlphaFold** |
| `parse_af3_results.py` | Result zips → `af3_scored.csv`, ΔAUC tests, corrected metric panel |
| `af3_report_stats.py` | Every grouping × metric → `af3_report_data.json` (uses the parser's own functions) |
| `af3_report_build.py` | `af3_report_data.json` → the 26-page PDF |
| `af3_idr_results_report.pdf` | Stage-1 results |
| `af3_all_statistics.csv` | All 496 cells, machine-readable |
| `af3_scored.csv` | Per-job scored table |
| `af3_idr_plan.md` | Scientific design, pre-registered criteria, stage-1 result, limitations. **Not yet updated with stage 2b** — read §5b here instead until it is |
| `af3_idr_handoff.md` | This file |
| `data/` | Raw result folders, ~493 MB + stage 2b, **gitignored** — re-download from the server if lost. New folders: `09_15_batch7_stage2b/` (P23–P28, N32–N38) and an updated `09_15_batch5/` (now includes N28's no-template rerun) |
| `docs/`, `viewer/` | Static GitHub Pages site (`viewer/build_static.py` → `docs/`) and the live dev server (`viewer/serve.py`). Rebuilt 15 Sep evening — 85 of 87 jobs now have results in `docs/api/` |

Build scripts (`build_dataset.py`, `verify.py`, `raw_fetch_1.fasta`) live in the session
workspace and are only needed if the set is rebuilt. Ask for them if so.

## 4. Do this next — 20 jobs, one day

1. **Submit batch 4** (S20, S21). Finishes the scramble arm at 10 pairs. The paired test is
   currently 6/8 with a sign-test p of 0.29 — two more pairs will not fix the power, but it
   completes the arm as designed.
2. **Templates off — done, 10/10.** P06, P07, P20, P21, N15, N16, N17, N26, N27 re-ran
   15 Sep with `"useStructureTemplate": false` (`af3_batch_05_notmpl.json`, results in
   `data/09_15_batch5/`, rows T01–T09 in the pairs CSV). **Result: 8 of 9 moved ≤0.05
   ipTM; templates are not driving the memorisation or the hard-zipper promiscuity.**
   N28 (FOS×FOS homodimer) followed the same evening — its result landed in the same
   folder and needed a new row (**T11**, since T10 was already used for the p53/MDM2
   ablation) because the parser matches by exact normalised job name first and only
   falls back to prefix-matching when there's no exact hit; without T11 the notmpl
   result would have silently been ignored while the original templated N28 row stayed
   in every stat. **N28: templated ipTM 0.31 → no-template 0.32, Δ+0.01** — same
   conclusion as the other 9, templates aren't the explanation. Full table in
   `af3_idr_plan.md` §"Stage 2 — what to run next" (not yet updated with the N28/T11 row).
3. **Seed replicates.** `"modelSeeds": ["2","3"]` on P06, P07, P20, P21 only — 8 jobs.
   Seed 1 already ran on 59 of 60 jobs, so seeds 2 and 3 are the only new information.
   This is the first actual measurement of run-to-run variance. **Now the highest-value
   experiment left**, since templates-off closed out item 2.
4. **Confirm or drop** the two `uncertain` negatives: H1.0 alone (N22) and ERD10 alone
   (N24). They currently sit in every headline number on unconfirmed labels.
5. **Decide how test 3 handles length, and write it down before looking again.** Either a
   length-matched subset or the partialled version the parser already asks for. Deciding
   after seeing the numbers is the thing the pre-registration exists to prevent.
6. **Monomer baseline — batch drafted, not submitted.** `af3_batch_06_monomers.json` (new,
   generated from `af3_idr_regions.csv`) has all 30 regions as single-chain jobs, names
   sanitised. Gives a pLDDT/compactness reference to subtract per chain — useful now that
   mean pLDDT is length-confounded (ρ = −0.60) and inverted in the headline test. Upload
   when there's server budget; nothing else depends on it.
7. **Expand `matched_internal` with non-fuzzy, structured IDR–IDR pairs — done, see §5b.**
   Landed 15 Sep evening as stage 2b / `af3_batch_07_stage2b.json`: p53 TAD × MDM2
   N-domain (1YCR) plus 4 new SLiM classes (LXXLL/NR, SH3/PxxP, 14-3-3/client,
   KEN-box/APC-C), each with a literature-verified matched-internal mutant or scramble
   negative. The PUMA/BAD-BH3×MCL-1 and c-Myc-TAD×Max/Bin1 candidates below were **not**
   built this round — still open if a third expansion round happens:
   - PUMA BH3 × MCL-1 (2ROC) vs BAD BH3 × MCL-1 (non-binder) — matched-internal pair on the
     MCL-1 side.
   - c-Myc TAD × Bin1 MBD or c-Myc TAD × Max (1NKP) as a second bHLHZ-family positive,
     matched against the existing MAX/GCN4/Fos/Jun hard-zipper negatives.

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

## 5b. What stage 2b found (15 September 2026 evening)

13 new jobs, all scored, `useStructureTemplate: false` throughout. Numbers from
`af3_scored.csv`; job IDs P23–P28 / N32–N38 / T11.

**p53–MDM2 was corrupted, now fixed — and the fix raises a new question.**
P22/T10 (full-length p53 × full-length MDM2) were run on sequences that were never
verified against UniProt and turned out to be corrupted exactly the way FMR1 was: p53
had a 1-residue insertion plus a substitution block, MDM2 had a ~38-residue runaway
repeat in its C-terminal Pro-rich tail. Both rows are re-labelled `evidence_strength =
superseded` (sequences and original numbers kept for the record) and P23 reruns the
same full-length pairing with corrected sequences: **ipTM 0.26 (corrupted) → 0.31
(corrected)**, chain_pair_pae_min 1.08 → 1.00. Real signal either way, slightly
stronger once the input isn't garbled. P24 then runs the *actual* crystallised
minimal interface (1YCR: p53 TAD × MDM2 N-domain) — the `demonstrated`-grade version
of the same pair — and it scores **lower**, not higher: ipTM 0.20, chain_pair_pae_min
5.37. The matched-internal negative (N32, p53 TAD L22Q/W23S × MDM2 N-domain) does
score lower still (ipTM 0.17, pae_min 6.80, 0 median contacts vs P24's 4), so the
ordering is correct, but a real solved co-crystal structure scoring below the messy
full-length "for fun" version is not explained by anything in this pipeline and is
reported as an open finding, not resolved. Plausibly a length/MSA-depth effect (P24's
two chains are only 61+109 residues) rather than anything about MDM2 specifically —
worth checking against the monomer baseline (`af3_batch_06_monomers.json`) if that
ever gets submitted.

**Four new SLiM-class matched-internal pairs, mixed results — reported plainly, not
softened:**

| Positive | ipTM | Negative (mutant/scramble) | ipTM | Direction | Note |
|---|---|---|---|---|---|
| P25 ESR1 LBD × NCOA2 LXXLL (1GWQ) | 0.95 | N33 × NCOA2 LXXAA dead mutant | 0.76 | correct | gap is smaller than hoped — AF3 still gives a motif-killed mutant substantial confidence, same flavour as the hard-zipper promiscuity in stage 1 |
| P26 ABL1 SH3 × SH3BP1 PxxP (1ABO) | 0.87 | N34 × scrambled peptide | 0.28 | correct | cleanest margin in this batch, but scrambles lose their MSA (documented confound) so part of the gap may be missing evolutionary information, not pure specificity |
| P27 14-3-3ζ × Raf-1 pS259 (1A37) | 0.73 | N35 × Raf-1 S259A | 0.70 | **no separation** | expected null: neither chain is phosphorylated (AF3 can't model pSer), so this mainly confirms the PTM caveat rather than testing real recognition |
| P28 CDC20 WD40 × BUB1B KEN-box (4GGD) | 0.68 | N36 × KEN→AAA dead mutant | 0.55 | correct | moderate gap |

**A cross-motif promiscuity hit, and ipTM is misleading on it.** N37 (ABL1 SH3 domain
offered the NCOA2 LXXLL peptide instead of its own PxxP partner — a nonsensical
domain/motif mismatch) scores ipTM **0.75**, nearly as high as the real SH3-PxxP pair
(P26, 0.87) and well above several of the real matched-internal negatives above. But
`contact_prob_max` is only 0.46 and median contact count is **0** — no localized
interface actually forms; ipTM alone would have called this a confident hit. This is
the same "AF3 stacks/docks disordered chains indiscriminately" artefact
`af3_idr_plan.md` documents for the hard coiled-coils and the homotypic LC domains,
now reproduced on new material, and a concrete illustration of why ipTM can't be
trusted alone (the reason the 16-metric panel exists). N38 (ESR1 LBD offered the
BUB1B KEN-box) is a milder case: ipTM 0.47 reads more like a real negative, though
contact_prob_max 0.79 with 21 contacts is not nothing either.

**Net for the matched-internal design**: 4 of 5 new SLiM pairs (P24/N32, P25/N33,
P26/N34, P28/N36) order correctly on ipTM; P27/N35 is a predicted null, not a failure.
That is a better hit rate than the pre-registered headline tests ever got, but two of
the four gaps are small (P25/N33, P28/N36) and one cross-pair (N37) reproduces the
exact generic-promiscuity failure mode from stage 1 on brand-new domains. Read
together with §"What stage 1 found" below, not as a replacement for it — none of this
has been folded into the pre-registered ΔAUC tests yet, which still pool everything
into the original 31 groupings and don't know `hetero_domain` exists as a category
(see "New controlled vocabulary" note just below).

**New controlled vocabulary, not yet reflected in `af3_idr_plan.md`:**
`interaction_class = hetero_domain` and `binding_mode = coupled_folding` (a SLiM
binds a stably folded receptor domain — distinct from the existing `hetero_fold`
class, which is two disordered chains mutually folding each other, e.g. ACTR/NCBD);
`neg_type = motif_dead_mutant` for a point-mutant matched-internal negative (as
opposed to the existing `homo_negative` reuse of an unmutated chain). `evidence_
strength = superseded` is also new, used only for P22/T10.

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

**The parser matches an exact normalised job name before it tries a prefix match** —
this bit stage 2b. N28's original templated result and its later no-template rerun
share the same job_name prefix (`N28_FOS_bZIP__FOS_bZIP`), and `results.get(key)` finds
the exact match (the original templated folder) before the `startswith` fallback that
handles genuine suffix variants (like `_notmpl` for T01–T09, which have no exact-match
collision) ever runs. A no-template rerun of an *already-run* job silently keeps
scoring the old templated result unless you give the rerun its own row with a distinct
`job_name` — that's why N28's rerun got row **T11** rather than overwriting N28's row.
If you rerun any other already-scored job (seed replicates, further ablations), give
it a new job_id/job_name pair too, or check `af3_scored.csv` isn't quietly still using
the old folder.

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
   N15, N16, N17, N28. **Fully resolved**: all 10 templates-off reruns (T01–T09, plus N28's
   as T11) move ≤0.05 ipTM, so the memorisation and hard-zipper findings are not template
   artefacts.
2. **New (stage 2b): the demonstrated-grade minimal p53-TAD×MDM2-N-domain pair (P24, the
   actual 1YCR interface) scores *below* the corrected full-length pair (P23) — ipTM 0.20
   vs 0.31 — despite being the better-evidenced construct.** The matched-internal ordering
   inside the minimal-domain pair is still correct (P24 > N32), so this isn't a labelling
   error, but nothing in this pipeline explains why the actual solved interface is less
   confident than the messy full-length version. Candidate explanation: P24's two chains
   are short (61+109 residues) and may get thinner MSAs than the full-length 393+491
   residue chains — check against `af3_batch_06_monomers.json` baselines if that batch
   ever runs, and don't generalise "AF3 prefers full-length" from an n=1 comparison.
3. **21 of 31 negatives rest on `inferred` or `absence_only` evidence.** "No paper reports
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
