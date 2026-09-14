#!/usr/bin/env python3
"""
Comprehensive statistics for the AF3 IDR-IDR benchmark report.

Every number is computed with the project's OWN functions -- parse_af3_results.auc,
bootstrap_auc, perm_p, benjamini_hochberg, delta_auc_vs_baseline and
baseline_leakage.loo_logistic / pair_feats -- so the report cannot drift from the
pipeline. The pre-registered five subsets use the identical definitions and default
seeds as parse_af3_results.main, so those cells reproduce its printout exactly.

Extends the pipeline in four ways, each flagged in the output:
  * an exploratory local-PAE block (ipSAE, LIS) on every grouping, stored apart from
    the 16-metric cells and BH-corrected against the 16 + 2 family
  * every grouping of the set (pair type, negative type, evidence grade, positive
    subgroup, LLPS stratum), not just the five pre-registered subsets
  * arm-specific composition baselines for the homotypic / heterotypic split that
    plan criterion 3 requires and the parser does not compute
  * length-adjusted (residualised on log total residues) sensitivity analysis for
    the metrics the parser flags as length-confounded

Writes af3_report_data.json next to this script.
"""
import sys, os, csv, json, math, glob, statistics, collections, re

sys.dont_write_bytecode = True                    # never litter the repo with .pyc
REPO = '/Users/barry/Documents/GitHub/AlphaFold3-IDP'
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO)

import numpy as np
from multiprocessing import Pool
import parse_af3_results as P
import baseline_leakage as B

MEM = {'P06', 'P07', 'P20', 'P21'}                # memorised pipeline controls
METRICS = list(P.METRICS)                         # [(name, sign)] x 16
METRIC_NAMES = [m for m, _ in METRICS]
LOCAL = list(P.EXPLORATORY)                       # ipSAE, LIS: added after stage 1
LOCAL_NAMES = [m for m, _ in LOCAL]
PRIMARY = P.PRIMARY

# ----------------------------------------------------------------- load
rows = list(csv.DictReader(open(os.path.join(HERE, 'af3_scored.csv'))))
pairs = list(csv.DictReader(open(os.path.join(REPO, 'af3_idr_pairs_1.csv'))))
pr = {r['job_name']: r for r in pairs}
for r in rows:
    p = pr[r['job_name']]
    for k in ('stage', 'chain_A', 'chain_B', 'len_A', 'len_B', 'evidence'):
        r[k] = p[k]
by_name = {r['job_name']: r for r in rows}

def fv(r, m):
    """Signed metric value, exactly as parse_af3_results.vals() forms it."""
    v = r.get(m)
    return None if v in ('', None) else float(v)

SIGN = dict(METRICS + LOCAL)

# ----------------------------------------------------------------- composition baseline
def loo_scores(sub):
    """(LOO AUC, {job_name: score}) for a subset of real-sequence rows."""
    feats = [B.pair_feats(pr[r['job_name']]['seq_A'], pr[r['job_name']]['seq_B']) for r in sub]
    names_f = sorted(feats[0])
    X = np.array([[f[k] for k in names_f] for f in feats], float)
    y = np.array([1 if r['label'] == 'positive' else 0 for r in sub])
    a, s = B.loo_logistic(X, y, return_scores=True)
    return float(a), {r['job_name']: float(v) for r, v in zip(sub, s)}, names_f

real = [r for r in rows if r['label'] != 'negative_scramble']          # 52
homo_real = [r for r in real if r['pair_type'] == 'homotypic']         # 25
het_real = [r for r in real if r['pair_type'] == 'heterotypic']        # 27

pooled_auc, base_pooled, feat_names = loo_scores(real)
homo_auc, base_homo, _ = loo_scores(homo_real)
het_auc, base_het, _ = loo_scores(het_real)
BASES = {'pooled': base_pooled, 'homotypic': base_homo, 'heterotypic': base_het}

# single-feature leakage table (two-sided, as baseline_leakage reports it)
def single_feature_table(sub):
    feats = [B.pair_feats(pr[r['job_name']]['seq_A'], pr[r['job_name']]['seq_B']) for r in sub]
    names_f = sorted(feats[0])
    y = [1 if r['label'] == 'positive' else 0 for r in sub]
    out = []
    for j, n in enumerate(names_f):
        pos = [feats[i][n] for i in range(len(sub)) if y[i] == 1]
        neg = [feats[i][n] for i in range(len(sub)) if y[i] == 0]
        out.append({'feature': n, 'auc_two_sided': B.two_sided(P.auc(pos, neg)),
                    'auc_raw': P.auc(pos, neg)})
    return sorted(out, key=lambda d: -d['auc_two_sided'])

# ----------------------------------------------------------------- length adjustment
logL = {r['job_name']: math.log(float(r['total_residues'])) for r in rows}
resid = {}
for m in METRIC_NAMES + LOCAL_NAMES:
    x = np.array([logL[r['job_name']] for r in real])
    y = np.array([fv(r, m) for r in real], float)
    A = np.vstack([x, np.ones_like(x)]).T
    slope, icept = np.linalg.lstsq(A, y, rcond=None)[0]
    resid[m] = {r['job_name']: fv(r, m) - (slope * logL[r['job_name']] + icept) for r in rows}

# ----------------------------------------------------------------- groups
def P_nonmem(r):  return r['label'] == 'positive' and r['job_id'] not in MEM
def N_real(r):    return r['label'] == 'negative'
def N_gold(r):    return r['label'] == 'negative' and 'demonstrated' in r['evidence_strength']
def N_scr(r):     return r['label'] == 'negative_scramble'

GROUPS = []
def G(family, name, pp, nn, base='pooled', note=''):
    GROUPS.append({'family': family, 'name': name, 'pos': pp, 'neg': nn,
                   'base': base, 'note': note})

# --- family 1: the five pre-registered subsets (definitions copied from the parser)
G('Pre-registered', '1. positives vs real negatives', P_nonmem, N_real, 'pooled',
  'primary test; headline')
G('Pre-registered', '2a. mutual folding vs hard coiled-coil',
  lambda r: r['job_id'] in MEM, lambda r: r['neg_type'] == 'hard_coiledcoil', 'pooled',
  'pipeline expectation: high')
G('Pre-registered', '2b. fuzzy positives vs real negatives',
  lambda r: r['label'] == 'positive' and r['binding_mode'] == 'fuzzy',
  lambda r: r['label'] == 'negative' and r['neg_type'] != 'hard_coiledcoil', 'pooled',
  'plan: "this is the result"')
G('Pre-registered', '3. positives vs gold negatives', P_nonmem, N_gold, 'pooled',
  'defensible labels, less power')
G('Pre-registered', '4. positives vs scrambles', P_nonmem, N_scr, None,
  'MSA-confounded; separate arm')

# --- family 2: pair type, with arm-specific baselines (plan criterion 3)
G('By pair type', 'Homotypic only',
  lambda r: P_nonmem(r) and r['pair_type'] == 'homotypic',
  lambda r: N_real(r) and r['pair_type'] == 'homotypic', 'homotypic',
  'arm baseline 0.85')
G('By pair type', 'Heterotypic only',
  lambda r: P_nonmem(r) and r['pair_type'] == 'heterotypic',
  lambda r: N_real(r) and r['pair_type'] == 'heterotypic', 'heterotypic',
  'arm baseline 0.53')

# --- family 3: negative subtype (plan limitation 7 asks for these separately)
for nt in ('noncognate', 'cross_kingdom', 'charge_repulsive', 'hard_coiledcoil',
           'homo_negative', 'scramble'):
    G('By negative type', f'vs {nt}', P_nonmem,
      (lambda t: (lambda r: r['neg_type'] == t))(nt),
      None if nt == 'scramble' else 'pooled')

# --- family 4: evidence grade of the negatives
EV = [('demonstrated', lambda r: r['evidence_strength'] == 'demonstrated'),
      ('demonstrated+matched_internal', lambda r: r['evidence_strength'] == 'demonstrated+matched_internal'),
      ('gold (both demonstrated grades)', lambda r: 'demonstrated' in r['evidence_strength']),
      ('inferred', lambda r: r['evidence_strength'] == 'inferred'),
      ('absence_only', lambda r: r['evidence_strength'] == 'absence_only'),
      ('uncertain', lambda r: r['evidence_strength'] == 'uncertain'),
      ('weak (inferred + absence_only)', lambda r: r['evidence_strength'] in ('inferred', 'absence_only'))]
for nm, f in EV:
    G('By negative evidence', f'vs {nm}', P_nonmem,
      (lambda f: (lambda r: N_real(r) and f(r)))(f), 'pooled')

# --- family 5: positive subgroup, all against the same 31 real negatives
G('By positive subgroup', 'Fuzzy heterotypic positives',
  lambda r: P_nonmem(r) and r['pair_type'] == 'heterotypic', N_real, 'pooled')
G('By positive subgroup', 'Fuzzy homotypic positives',
  lambda r: P_nonmem(r) and r['pair_type'] == 'homotypic', N_real, 'pooled')
G('By positive subgroup', 'Mutual-folding (memorised) positives',
  lambda r: r['job_id'] in MEM, N_real, 'pooled', 'pipeline control, not evidence')
G('By positive subgroup', 'All 21 positives incl. memorised',
  lambda r: r['label'] == 'positive', N_real, 'pooled')
G('By positive subgroup', 'LLPS-reported positives',
  lambda r: P_nonmem(r) and r['llps'] == 'yes', N_real, 'pooled',
  'llps is stratifying metadata, never an outcome')
G('By positive subgroup', 'Positives without LLPS report',
  lambda r: P_nonmem(r) and r['llps'] != 'yes', N_real, 'pooled',
  'llps is stratifying metadata, never an outcome')

# --- family 6: length-adjusted sensitivity (uses residual metric values)
G('Length-adjusted', '1. positives vs real negatives (adj.)', P_nonmem, N_real, 'pooled')
G('Length-adjusted', '2b. fuzzy positives vs real negatives (adj.)',
  lambda r: r['label'] == 'positive' and r['binding_mode'] == 'fuzzy',
  lambda r: r['label'] == 'negative' and r['neg_type'] != 'hard_coiledcoil', 'pooled')
G('Length-adjusted', '3. positives vs gold negatives (adj.)', P_nonmem, N_gold, 'pooled')
G('Length-adjusted', 'Homotypic only (adj.)',
  lambda r: P_nonmem(r) and r['pair_type'] == 'homotypic',
  lambda r: N_real(r) and r['pair_type'] == 'homotypic', 'homotypic')
G('Length-adjusted', 'Heterotypic only (adj.)',
  lambda r: P_nonmem(r) and r['pair_type'] == 'heterotypic',
  lambda r: N_real(r) and r['pair_type'] == 'heterotypic', 'heterotypic')

# ----------------------------------------------------------------- build cells
def cell_task(g, metric, adjusted):
    src = resid[metric] if adjusted else None
    sign = SIGN[metric]
    def value(r):
        v = src[r['job_name']] if adjusted else fv(r, metric)
        return None if v is None else sign * v
    pos_rows = [r for r in rows if g['pos'](r)]
    neg_rows = [r for r in rows if g['neg'](r)]
    pos = [value(r) for r in pos_rows if value(r) is not None]
    neg = [value(r) for r in neg_rows if value(r) is not None]
    base = BASES.get(g['base']) if g['base'] else None
    bp = bn = None
    if base is not None:
        bp = [base[r['job_name']] for r in pos_rows if r['job_name'] in base]
        bn = [base[r['job_name']] for r in neg_rows if r['job_name'] in base]
        if not (len(bp) == len(pos) and len(bn) == len(neg) and bp and bn):
            bp = bn = None
    return {'group': g['name'], 'family': g['family'], 'metric': metric,
            'pos': pos, 'neg': neg, 'bp': bp, 'bn': bn}

def compute(t):
    pos, neg, bp, bn = t['pos'], t['neg'], t['bp'], t['bn']
    out = {k: t[k] for k in ('group', 'family', 'metric')}
    out['n_pos'], out['n_neg'] = len(pos), len(neg)
    if not pos or not neg:
        out['auc'] = None
        return out
    out['auc'] = P.auc(pos, neg)
    small = min(len(pos), len(neg)) < 3
    out['small_n'] = small
    lo, hi = P.bootstrap_auc(pos, neg)
    out['ci'] = [lo, hi]
    out['p'] = None if small else P.perm_p(pos, neg)
    if bp is not None:
        out['base_auc'] = P.auc(bp, bn)
        d, dlo, dhi = P.delta_auc_vs_baseline(pos, neg, bp, bn)
        out['delta'] = d
        out['delta_ci'] = [dlo, dhi]
        out['beats_baseline'] = bool(dlo > 0)
    return out

tasks, local_tasks = [], []
for g in GROUPS:
    adjusted = g['family'] == 'Length-adjusted'
    for m in METRIC_NAMES:
        tasks.append(cell_task(g, m, adjusted))
    for m in LOCAL_NAMES:
        local_tasks.append(cell_task(g, m, adjusted))

if __name__ == '__main__':
    with Pool() as pool:
        cells = pool.map(compute, tasks, chunksize=4)
        local_cells = pool.map(compute, local_tasks, chunksize=4)

    # BH correction across the 16 metrics within each group, as the parser does
    bygroup = collections.defaultdict(list)
    for c in cells:
        bygroup[c['group']].append(c)
    for gname, cs in bygroup.items():
        testable = [c for c in cs if c.get('p') is not None]
        if testable:
            q = P.benjamini_hochberg([c['p'] for c in testable])
            for c, qq in zip(testable, q):
                c['q'] = qq

    # Local PAE cells are corrected against the group's 16 panel p-values plus their own,
    # so metrics added after the fact pay for the bigger family; panel q is left as is.
    local_by_group = collections.defaultdict(list)
    for c in local_cells:
        local_by_group[c['group']].append(c)
    for gname, lcs in local_by_group.items():
        new = [c for c in lcs if c.get('p') is not None]
        if new:
            panel = [c['p'] for c in bygroup[gname] if c.get('p') is not None]
            q = P.benjamini_hochberg(panel + [c['p'] for c in new])
            for c, qq in zip(new, q[len(panel):]):
                c['q'] = qq

    # ---------------------------------------------------- reference predictors per group
    refs = []
    for g in GROUPS:
        pos_rows = [r for r in rows if g['pos'](r)]
        neg_rows = [r for r in rows if g['neg'](r)]
        rec = {'group': g['name'], 'family': g['family'], 'note': g['note'],
               'n_pos': len(pos_rows), 'n_neg': len(neg_rows),
               'base_key': g['base']}
        Lp = [float(r['total_residues']) for r in pos_rows]
        Ln = [float(r['total_residues']) for r in neg_rows]
        if Lp and Ln:
            rec['length_auc'] = P.auc(Lp, Ln)
            rec['len_pos_median'] = statistics.median(Lp)
            rec['len_neg_median'] = statistics.median(Ln)
            rec['len_pos_range'] = [min(Lp), max(Lp)]
            rec['len_neg_range'] = [min(Ln), max(Ln)]
        base = BASES.get(g['base']) if g['base'] else None
        if base:
            bp = [base[r['job_name']] for r in pos_rows if r['job_name'] in base]
            bn = [base[r['job_name']] for r in neg_rows if r['job_name'] in base]
            if bp and bn:
                rec['baseline_auc'] = P.auc(bp, bn)
        refs.append(rec)

    # ---------------------------------------------------- matched-internal (descriptive)
    matched = [{'positive': 'P06', 'negatives': ['N26', 'N27'],
                'why': 'ACTR and NCBD each fold only with the other'},
               {'positive': 'P07', 'negatives': ['N28'],
                'why': 'c-Fos forms no stable homodimer (O\'Shea 1989)'}]
    jid = {r['job_id']: r for r in rows}
    for m in matched:
        m['rows'] = []
        for j in [m['positive']] + m['negatives']:
            r = jid[j]
            m['rows'].append({'job_id': j, 'job_name': r['job_name'], 'label': r['label'],
                              **{k: fv(r, k) for k in METRIC_NAMES}})
        m['correct'] = {}
        for k in METRIC_NAMES:
            pv = fv(jid[m['positive']], k) * SIGN[k]
            m['correct'][k] = all(pv > fv(jid[n], k) * SIGN[k] for n in m['negatives'])

    # ---------------------------------------------------- scramble arm, paired
    scr_pairs = []
    for r in rows:
        if r['label'] != 'negative_scramble':
            continue
        parent_id = 'P' + r['job_id'][1:]
        if parent_id not in jid:
            continue
        p_row, s_row = jid[parent_id], r
        scr_pairs.append({
            'scramble': r['job_id'], 'parent': parent_id,
            'scramble_name': r['job_name'], 'parent_name': p_row['job_name'],
            'memorised_parent': parent_id in MEM,
            'delta': {k: (fv(p_row, k) - fv(s_row, k)) for k in METRIC_NAMES},
            'parent_vals': {k: fv(p_row, k) for k in METRIC_NAMES},
            'scramble_vals': {k: fv(s_row, k) for k in METRIC_NAMES}})
    scr_summary = {}
    for k in METRIC_NAMES:
        wins = sum(1 for d in scr_pairs if d['delta'][k] * SIGN[k] > 0)
        n = len(scr_pairs)
        # exact two-sided sign test
        p = min(1.0, 2 * sum(math.comb(n, i) for i in range(max(wins, n - wins), n + 1)) / 2 ** n)
        scr_summary[k] = {'wins': wins, 'n': n, 'sign_p': p,
                          'median_delta': statistics.median(d['delta'][k] for d in scr_pairs)}

    # ---------------------------------------------------- confounds & descriptives
    conf = []
    for m in METRIC_NAMES:
        v = [(fv(r, m), float(r['total_residues'])) for r in rows if fv(r, m) is not None]
        rho = P.spearman([a for a, _ in v], [b for _, b in v])
        conf.append({'metric': m, 'rho_length': rho, 'flagged': abs(rho) > 0.5})

    # per-class metric distributions
    CLASSES = [('Positive, fuzzy', lambda r: P_nonmem(r)),
               ('Positive, mutual folding', lambda r: r['job_id'] in MEM),
               ('Negative, noncognate', lambda r: r['neg_type'] == 'noncognate'),
               ('Negative, cross-kingdom', lambda r: r['neg_type'] == 'cross_kingdom'),
               ('Negative, charge-repulsive', lambda r: r['neg_type'] == 'charge_repulsive'),
               ('Negative, hard coiled-coil', lambda r: r['neg_type'] == 'hard_coiledcoil'),
               ('Negative, homotypic', lambda r: r['neg_type'] == 'homo_negative'),
               ('Scramble', N_scr)]
    dists = []
    for nm, f in CLASSES:
        sel = [r for r in rows if f(r)]
        d = {'class': nm, 'n': len(sel),
             'jobs': [r['job_id'] for r in sel]}
        for m in METRIC_NAMES + ['total_residues']:
            vals = [fv(r, m) for r in sel if fv(r, m) is not None]
            if vals:
                d[m] = {'median': statistics.median(vals), 'min': min(vals), 'max': max(vals),
                        'values': vals}
        dists.append(d)

    # raw summary-confidence extras the scored table drops (has_clash, chain_ptm, recycles)
    extras = {}
    for d in sorted(glob.glob(f'{REPO}/data/*/*/')):
        name = os.path.basename(d.rstrip('/'))
        S = [json.load(open(x)) for x in sorted(glob.glob(d + '*summary_confidences_*.json'))]
        if not S:
            continue
        key = P.norm(name)
        match = [r for r in rows if P.norm(r['job_name']) == key]
        if not match:
            continue
        extras[match[0]['job_id']] = {
            'has_clash_any': any(bool(s.get('has_clash')) for s in S),
            'chain_ptm_mean': float(np.mean([np.mean(s['chain_ptm']) for s in S if s.get('chain_ptm')])),
            'num_recycles': sorted({s.get('num_recycles') for s in S if s.get('num_recycles') is not None}),
            'iptm_samples': [s.get('iptm') for s in S],
            'iptm_sample_sd': float(np.std([s['iptm'] for s in S])) if len(S) > 1 else 0.0}

    # MSA depth + templates per job, recomputed here so the report is data-driven
    msa_tmpl = {}
    for d in sorted(glob.glob(f'{REPO}/data/*/*/')):
        name = os.path.basename(d.rstrip('/'))
        key = P.norm(name)
        match = [r for r in rows if P.norm(r['job_name']) == key]
        if not match:
            continue
        tm = collections.defaultdict(list)
        for f in sorted(glob.glob(d + 'templates/*template_hit_*_chains_*.cif')):
            mm = re.search(r'_hit_(\d+)_chains_([a-z_]+)\.cif$', os.path.basename(f))
            if mm:
                tm[mm.group(2)].append(open(f).readline().strip().replace('data_', ''))
        ms = {}
        for f in sorted(glob.glob(d + 'msas/*.a3m')):
            mm = re.search(r'_(paired|unpaired)_msa_chains_([a-z_]+)\.a3m$', os.path.basename(f))
            if mm:
                ms[f'{mm.group(2)}/{mm.group(1)}'] = sum(
                    1 for line in open(f, errors='ignore') if line.startswith('>'))
        jr = glob.glob(d + '*job_request.json')
        seeds, use_tmpl = None, None
        if jr:
            j = json.load(open(jr[0]))
            j = j[0] if isinstance(j, list) else j
            seeds = j.get('modelSeeds')
            use_tmpl = [c.get('proteinChain', {}).get('useStructureTemplate') for c in j['sequences']]
        msa_tmpl[match[0]['job_id']] = {'templates': dict(tm), 'msa': ms,
                                        'seeds': seeds, 'use_template': use_tmpl}

    # ---------------------------------------------------- per-job table
    table = []
    for r in sorted(rows, key=lambda r: r['job_id']):
        rec = {k: r[k] for k in ('job_id', 'job_name', 'label', 'pair_type', 'interaction_class',
                                 'binding_mode', 'llps', 'in_pdb', 'neg_type',
                                 'evidence_strength', 'stage', 'total_residues',
                                 'chain_A', 'chain_B', 'len_A', 'len_B', 'n_seeds')}
        rec.update({m: fv(r, m) for m in METRIC_NAMES + LOCAL_NAMES})
        rec['iptm_sd'] = fv(r, 'iptm_sd')
        rec['baseline_pooled'] = base_pooled.get(r['job_name'])
        rec['baseline_arm'] = (base_homo if r['pair_type'] == 'homotypic' else base_het).get(r['job_name'])
        rec.update(extras.get(r['job_id'], {}))
        table.append(rec)

    # ---------------------------------------------------- local PAE block (exploratory)
    # ipSAE and LIS were added after stage 1 was scored. Everything about them lives
    # under 'local_pae' so no existing key, count or page changes meaning.
    local_conf = []
    for m in LOCAL_NAMES:
        v = [(fv(r, m), float(r['total_residues'])) for r in rows if fv(r, m) is not None]
        rho = P.spearman([a for a, _ in v], [b for _, b in v])
        local_conf.append({'metric': m, 'rho_length': rho, 'flagged': abs(rho) > 0.5})

    local_scr_pairs = [{'scramble': s['scramble'], 'parent': s['parent'],
                        'memorised_parent': s['memorised_parent'],
                        'parent_vals': {k: fv(jid[s['parent']], k) for k in LOCAL_NAMES},
                        'scramble_vals': {k: fv(jid[s['scramble']], k) for k in LOCAL_NAMES}}
                       for s in scr_pairs]
    local_scr_summary = {}
    for k in LOCAL_NAMES:
        diffs = [s['parent_vals'][k] - s['scramble_vals'][k] for s in local_scr_pairs]
        # parent and scramble often both score exactly 0, so ties are dropped (sign test)
        untied = [x for x in diffs if x != 0]
        wins, n = sum(1 for x in untied if x * SIGN[k] > 0), len(untied)
        p = (min(1.0, 2 * sum(math.comb(n, i) for i in range(max(wins, n - wins), n + 1)) / 2 ** n)
             if n else 1.0)
        local_scr_summary[k] = {'wins': wins, 'n': n, 'ties': len(diffs) - n, 'sign_p': p,
                                'median_delta': statistics.median(diffs)}

    local_matched = []
    for m in matched:
        js = [m['positive']] + m['negatives']
        mv = {j: {k: fv(jid[j], k) for k in LOCAL_NAMES} for j in js}
        local_matched.append({
            'positive': m['positive'], 'negatives': m['negatives'], 'why': m['why'], 'vals': mv,
            'correct': {k: all(mv[m['positive']][k] * SIGN[k] > mv[neg][k] * SIGN[k]
                               for neg in m['negatives']) for k in LOCAL_NAMES}})

    out = {
        'meta': {
            'generated': '2026-09-14',
            'n_jobs_expected': len(pairs), 'n_jobs_scored': len(rows),
            'missing': [p['job_name'] for p in pairs if p['job_name'] not in by_name],
            'primary_metric': PRIMARY,
            'metrics': [{'name': m, 'sign': s} for m, s in METRICS],
            'memorised': sorted(MEM),
            'baseline_features': feat_names,
            'baseline_loo': {'pooled_raw': pooled_auc, 'homotypic_raw': homo_auc,
                             'heterotypic_raw': het_auc,
                             'pooled_two_sided': B.two_sided(pooled_auc),
                             'homotypic_two_sided': B.two_sided(homo_auc),
                             'heterotypic_two_sided': B.two_sided(het_auc)},
            'baseline_n': {'pooled': len(real), 'homotypic': len(homo_real),
                           'heterotypic': len(het_real)},
        },
        'single_features': {
            'all': single_feature_table(real)[:8],
            'homotypic': single_feature_table(homo_real)[:8],
            'heterotypic': single_feature_table(het_real)[:8]},
        'cells': cells,
        'refs': refs,
        'matched_internal': matched,
        'scramble_pairs': scr_pairs,
        'scramble_summary': scr_summary,
        'confounds': conf,
        'distributions': dists,
        'msa_templates': msa_tmpl,
        'table': table,
        'local_pae': {
            'metrics': [{'name': m, 'sign': s} for m, s in LOCAL],
            'cutoffs_angstrom': {'ipsae': P.IPSAE_CUTOFF, 'lis': P.LIS_CUTOFF},
            'q_family': 'BH over the 16 panel metrics plus ipSAE and LIS, within each group',
            'cells': local_cells,
            'confounds': local_conf,
            'scramble_pairs': local_scr_pairs,
            'scramble_summary': local_scr_summary,
            'matched_internal': local_matched,
        },
    }
    with open(os.path.join(HERE, 'af3_report_data.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)

    # ------------------------------------------------ console check of key numbers
    def show(gname, metric='iptm', pool_=None):
        c = [x for x in (cells if pool_ is None else pool_)
             if x['group'] == gname and x['metric'] == metric]
        if not c:
            return
        c = c[0]
        s = f"  {gname:46s} AUC {c['auc']:.3f}"
        if c.get('ci'):
            s += f" CI [{c['ci'][0]:.3f},{c['ci'][1]:.3f}]"
        if 'delta' in c:
            s += f"  base {c['base_auc']:.3f}  d {c['delta']:+.3f} [{c['delta_ci'][0]:+.3f},{c['delta_ci'][1]:+.3f}]"
            s += '  PASS' if c['beats_baseline'] else '  FAIL'
        print(s + f"  n={c['n_pos']}/{c['n_neg']}")

    print(f"baseline LOO raw: pooled {pooled_auc:.3f}  homo {homo_auc:.3f}  het {het_auc:.3f}")
    print('\nPRE-REGISTERED (ipTM) -- must match parse_af3_results output:')
    for g in GROUPS:
        if g['family'] == 'Pre-registered':
            show(g['name'])
    print('\nBY PAIR TYPE (arm baselines):')
    for g in GROUPS:
        if g['family'] == 'By pair type':
            show(g['name'])
    print('\nBEST SECONDARY (contact_prob_sum) vs baseline:')
    for g in GROUPS:
        if g['family'] in ('Pre-registered', 'By pair type'):
            show(g['name'], 'contact_prob_sum')
    print('\nLENGTH-ONLY AUC per group:')
    for r in refs:
        if 'length_auc' in r:
            print(f"  {r['group']:46s} length {r['length_auc']:.3f}"
                  + (f"  baseline {r['baseline_auc']:.3f}" if 'baseline_auc' in r else ''))
    print('\nmatched-internal correct orderings (ipTM):',
          {m['positive']: m['correct']['iptm'] for m in matched})
    print('scramble arm, parent beats scramble (ipTM):',
          scr_summary['iptm']['wins'], '/', scr_summary['iptm']['n'],
          f"sign p={scr_summary['iptm']['sign_p']:.3f}")
    print('has_clash anywhere:', sum(1 for v in extras.values() if v['has_clash_any']), '/', len(extras))
    print('\nLOCAL PAE, exploratory (q = BH over the 16 panel metrics + these 2):')
    for g in GROUPS:
        if g['family'] in ('Pre-registered', 'By pair type'):
            for m in LOCAL_NAMES:
                show(g['name'], m, local_cells)
    print('local PAE rho vs length:', {c['metric']: round(c['rho_length'], 3) for c in local_conf})
    print('local PAE, parent beats scramble:',
          {k: f"{v['wins']}/{v['n']} (+{v['ties']} ties) p={v['sign_p']:.3f}"
           for k, v in local_scr_summary.items()})
    print('local PAE matched-internal correct:',
          {m['positive']: m['correct'] for m in local_matched})
    print(f"\nwrote af3_report_data.json  ({len(cells)} + {len(local_cells)} local cells, "
          f"{len(GROUPS)} groups)")
