#!/usr/bin/env python3
"""
Parse AlphaFold Server results for the IDR-IDR co-folding benchmark.

    python3 parse_af3_results.py <results_dir> [af3_idr_pairs.csv] [--no-full-data]

<results_dir> holds the zips (or unzipped folders) AlphaFold Server gives you. Job
folders are matched back to the pair table by the job_name you uploaded.

Reads BOTH confidence files per seed:

  *summary_confidences*.json   iptm, ptm, ranking_score, fraction_disordered,
                               has_clash, chain_pair_iptm, chain_pair_pae_min,
                               chain_ptm, chain_iptm, chain_ids
  *full_data*.json             pae [ntok,ntok], contact_probs [ntok,ntok],
                               atom_plddts [natom], token_chain_ids, atom_chain_ids

and derives the interchain-block metrics that matter when ipTM is floored.
Pass --no-full-data to skip the large files (summary metrics only, much faster).

Writes af3_scored.csv and prints the pre-registered tests plus the secondary panel
with Benjamini-Hochberg correction.

Requires numpy.
"""
import csv, json, os, re, sys, zipfile, statistics
from collections import defaultdict

import numpy as np

CONTACT_P = 0.5      # "a contact" for counting purposes
PAE_CONF = 10.0      # Angstrom; interchain token pairs below this are "confident"
PAE_TIGHT = 5.0

# ---------------------------------------------------------------- file discovery

def _walk_members(results_dir):
    """Yield (job_key, filename, loader) for every json in every job, zipped or not."""
    for root, _dirs, files in os.walk(results_dir):
        for fn in files:
            path = os.path.join(root, fn)
            if fn.endswith('.zip'):
                try:
                    z = zipfile.ZipFile(path)
                except zipfile.BadZipFile:
                    print(f'  ! not a zip: {path}', file=sys.stderr)
                    continue
                key = os.path.splitext(fn)[0]
                for member in z.namelist():
                    if member.endswith('.json'):
                        yield key, member, (lambda zz=z, m=member: json.load(zz.open(m)))
            elif fn.endswith('.json'):
                yield os.path.basename(root), fn, (lambda p=path: json.load(open(p)))


def norm(s):
    """AF Server rewrites punctuation in job names; compare on alphanumerics only."""
    return re.sub(r'[^a-z0-9]', '', s.lower())


def seed_of(name):
    """Group a summary file and its full_data sibling by their shared seed/model index."""
    m = re.search(r'(?:seed[-_]?)?(\d+)\.json$', os.path.basename(name))
    return m.group(1) if m else os.path.basename(name)


def collect(results_dir, want_full):
    """job_key -> seed -> {'summary': dict, 'full': dict}"""
    out = defaultdict(lambda: defaultdict(dict))
    for key, member, load in _walk_members(results_dir):
        base = os.path.basename(member)
        if 'summary_confidences' in base:
            out[norm(key)][seed_of(base)]['summary'] = load()
        elif want_full and 'full_data' in base:
            out[norm(key)][seed_of(base)]['full'] = load()
    return out


# ---------------------------------------------------------------- metrics

def _offdiag(mat):
    a = np.asarray(mat, dtype=float)
    iu = ~np.eye(a.shape[0], dtype=bool)
    return a[iu]


def summary_metrics(d):
    """Scalars available without the big files."""
    m = {
        'iptm': d.get('iptm'),
        'ptm': d.get('ptm'),
        'ranking_score': d.get('ranking_score'),
        'fraction_disordered': d.get('fraction_disordered'),
        'has_clash': float(bool(d.get('has_clash'))) if d.get('has_clash') is not None else None,
    }
    if d.get('chain_pair_iptm'):
        m['chain_pair_iptm'] = float(_offdiag(d['chain_pair_iptm']).mean())
    if d.get('chain_pair_pae_min'):
        m['chain_pair_pae_min'] = float(_offdiag(d['chain_pair_pae_min']).min())
    if d.get('chain_ptm'):
        m['chain_ptm_mean'] = float(np.mean(d['chain_ptm']))
    return m


def interchain_metrics(full):
    """
    Everything derived from the token x token blocks. These are the metrics that can
    still separate pairs when ipTM is floored near zero for every disordered complex.
    """
    chains = np.asarray(full['token_chain_ids'])
    uniq = list(dict.fromkeys(chains.tolist()))
    if len(uniq) < 2:
        return {}
    A, B = chains == uniq[0], chains == uniq[1]

    out = {}

    pae = np.asarray(full['pae'], dtype=float)
    block = np.concatenate([pae[np.ix_(A, B)].ravel(), pae[np.ix_(B, A)].ravel()])
    out['pae_inter_mean'] = float(block.mean())
    out['pae_inter_min'] = float(block.min())
    out['pae_inter_p05'] = float(np.percentile(block, 5))
    # fractions are length-normalised; min is a max-statistic and grows with chain length
    out['pae_inter_frac_lt10'] = float((block < PAE_CONF).mean())
    out['pae_inter_frac_lt5'] = float((block < PAE_TIGHT).mean())

    cp = np.asarray(full['contact_probs'], dtype=float)
    cblock = cp[np.ix_(A, B)]
    out['contact_prob_sum'] = float(cblock.sum())
    out['contact_prob_max'] = float(cblock.max())
    out['n_contacts_p50'] = int((cblock > CONTACT_P).sum())
    # normalised so it does not simply track chain length
    out['contact_prob_density'] = float(cblock.sum() / cblock.size)

    # "few strong vs many weak": participation ratio of the interchain contact mass.
    # Near 1 = one dominant contact; large = mass spread over many weak contacts.
    p = cblock.ravel()
    s = p.sum()
    out['contact_participation_ratio'] = float(s ** 2 / np.square(p).sum()) if s > 0 else 0.0

    if 'atom_plddts' in full and 'atom_chain_ids' in full:
        plddt = np.asarray(full['atom_plddts'], dtype=float)
        ach = np.asarray(full['atom_chain_ids'])
        out['plddt_mean'] = float(plddt.mean())
        out['plddt_chainA_mean'] = float(plddt[ach == uniq[0]].mean())
        out['plddt_chainB_mean'] = float(plddt[ach == uniq[1]].mean())
    return out


# ordering matters: sign says which direction means "more interaction"
METRICS = [
    ('iptm', +1), ('chain_pair_iptm', +1), ('ranking_score', +1), ('ptm', +1),
    ('fraction_disordered', -1),
    ('chain_pair_pae_min', -1), ('pae_inter_mean', -1), ('pae_inter_p05', -1),
    ('pae_inter_frac_lt10', +1), ('pae_inter_frac_lt5', +1),
    ('contact_prob_sum', +1), ('contact_prob_density', +1), ('contact_prob_max', +1),
    ('n_contacts_p50', +1), ('contact_participation_ratio', +1),
    ('plddt_mean', +1),
]
PRIMARY = 'iptm'


# ---------------------------------------------------------------- statistics

def auc(pos, neg):
    """Mann-Whitney AUC with tie correction."""
    if not pos or not neg:
        return float('nan')
    allv = sorted(pos + neg)
    ranks, i = {}, 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1] == allv[i]:
            j += 1
        ranks[allv[i]] = (i + j) / 2 + 1
        i = j + 1
    n1, n2 = len(pos), len(neg)
    return (sum(ranks[v] for v in pos) - n1 * (n1 + 1) / 2) / (n1 * n2)


def bootstrap_auc(pos, neg, n=10000, seed=0):
    import random
    rng = random.Random(seed)
    vals = sorted(auc([rng.choice(pos) for _ in pos], [rng.choice(neg) for _ in neg])
                  for _ in range(n))
    return vals[int(0.025 * n)], vals[int(0.975 * n)]


def perm_p(pos, neg, n=10000, seed=1):
    """Two-sided permutation p-value for AUC != 0.5."""
    import random
    rng = random.Random(seed)
    obs = abs(auc(pos, neg) - 0.5)
    pool, n1 = pos + neg, len(pos)
    hits = 0
    for _ in range(n):
        rng.shuffle(pool)
        if abs(auc(pool[:n1], pool[n1:]) - 0.5) >= obs:
            hits += 1
    return (hits + 1) / (n + 1)


def benjamini_hochberg(pvals):
    """Return BH-adjusted p-values, same order as input."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj, prev = [0.0] * m, 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        prev = min(prev, pvals[i] * m / k)
        adj[i] = prev
    return adj


def delta_auc_vs_baseline(af3_pos, af3_neg, base_pos, base_neg, n=5000, seed=3):
    """
    Paired bootstrap on AUC(AF3) - AUC(composition baseline), resampling the SAME jobs
    for both score vectors. Baseline scores are out-of-sample (leave-one-out), so this
    asks whether AF3 adds information beyond what sequence composition already gives.
    """
    import random
    rng = random.Random(seed)
    obs = auc(af3_pos, af3_neg) - auc(base_pos, base_neg)
    diffs = []
    for _ in range(n):
        ip = [rng.randrange(len(af3_pos)) for _ in af3_pos]
        ineg = [rng.randrange(len(af3_neg)) for _ in af3_neg]
        diffs.append(auc([af3_pos[i] for i in ip], [af3_neg[i] for i in ineg])
                     - auc([base_pos[i] for i in ip], [base_neg[i] for i in ineg]))
    diffs.sort()
    return obs, diffs[int(0.025 * n)], diffs[int(0.975 * n)]


def spearman(x, y):
    rx = np.argsort(np.argsort(np.asarray(x, dtype=float)))
    ry = np.argsort(np.argsort(np.asarray(y, dtype=float)))
    return float(np.corrcoef(rx, ry)[0, 1])


# ---------------------------------------------------------------- main

def main(results_dir, pairs_csv='af3_idr_pairs.csv', want_full=True):
    pairs = list(csv.DictReader(open(pairs_csv)))
    results = collect(results_dir, want_full)
    print(f'{len(results)} result jobs found for {len(pairs)} expected\n')

    rows, missing = [], []
    for p in pairs:
        key = norm(p['job_name'])
        seeds = results.get(key)
        if seeds is None:
            hits = [k for k in results if k.startswith(key)]
            seeds = results[hits[0]] if len(hits) == 1 else None
        if not seeds:
            missing.append(p['job_name'])
            continue

        per_seed = []
        for _sd, files in sorted(seeds.items()):
            if 'summary' not in files:
                continue
            m = summary_metrics(files['summary'])
            if 'full' in files:
                try:
                    m.update(interchain_metrics(files['full']))
                except (KeyError, ValueError) as e:
                    print(f"  ! {p['job_name']}: full_data unusable ({e})", file=sys.stderr)
            per_seed.append(m)
        if not per_seed:
            missing.append(p['job_name'])
            continue

        row = {k: p[k] for k in ('job_id', 'job_name', 'label', 'pair_type',
                                 'interaction_class', 'binding_mode', 'llps',
                                 'in_pdb', 'neg_type', 'evidence_strength',
                                 'total_residues')}
        row['n_seeds'] = len(per_seed)
        row['has_full_data'] = int(any('pae_inter_mean' in m for m in per_seed))
        # rank seeds by the model's own ranking_score so "best" means AF3's best
        best_i = max(range(len(per_seed)),
                     key=lambda i: per_seed[i].get('ranking_score') or -1e9)
        for name, sign in METRICS:
            vals = [m[name] for m in per_seed if m.get(name) is not None]
            if not vals:
                row[name] = row[f'{name}_sd'] = ''
                continue
            row[name] = round(per_seed[best_i].get(name, vals[0]), 5)
            row[f'{name}_sd'] = round(statistics.pstdev(vals), 5) if len(vals) > 1 else ''
        rows.append(row)

    if missing:
        print(f'MISSING RESULTS for {len(missing)} jobs:')
        for m in missing:
            print('   ', m)
        print()
    if not rows:
        print('nothing to score'); return

    with open('af3_scored.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    nfull = sum(r['has_full_data'] for r in rows)
    print(f'wrote af3_scored.csv  ({len(rows)} jobs, {nfull} with full_data metrics)\n')

    # ------------------------------------------------ subsets
    memorised = {'P06', 'P07', 'P20', 'P21'}
    SUBSETS = {
        '1. positives vs real negatives':
            (lambda r: r['label'] == 'positive' and r['job_id'] not in memorised,
             lambda r: r['label'] == 'negative'),
        '2a. mutual-folding vs hard coiled-coil':
            (lambda r: r['job_id'] in memorised,
             lambda r: r['neg_type'] == 'hard_coiledcoil'),
        '2b. FUZZY positives vs real negatives':
            (lambda r: r['label'] == 'positive' and r['binding_mode'] == 'fuzzy',
             lambda r: r['label'] == 'negative' and r['neg_type'] != 'hard_coiledcoil'),
        '3. positives vs GOLD negatives only':
            (lambda r: r['label'] == 'positive' and r['job_id'] not in memorised,
             lambda r: r['label'] == 'negative'
                       and 'demonstrated' in r.get('evidence_strength', '')),
        '4. positives vs scrambles (MSA-confounded)':
            (lambda r: r['label'] == 'positive' and r['job_id'] not in memorised,
             lambda r: r['label'] == 'negative_scramble'),
    }

    def vals(metric, pred, sign):
        return [sign * float(r[metric]) for r in rows if r.get(metric) not in ('', None) and pred(r)]

    # ------------------------------------------------ composition baseline per job
    # The comparator is NOT 0.5. Composition alone separates these labels (LC domains
    # that self-associate really are aromatic-rich and uncharged), so the question is
    # whether AF3 beats that, not whether it beats chance.
    base = {}
    try:
        from baseline_leakage import pair_feats, loo_logistic
        pr = {r['job_name']: r for r in pairs}
        sub = [r for r in rows if r['label'] != 'negative_scramble' and r['job_name'] in pr]
        if len(sub) > 10:
            names_f = sorted(pair_feats(pr[sub[0]['job_name']]['seq_A'],
                                        pr[sub[0]['job_name']]['seq_B']))
            Xb = np.array([[pair_feats(pr[r['job_name']]['seq_A'],
                                       pr[r['job_name']]['seq_B'])[k] for k in names_f]
                           for r in sub], float)
            yb = np.array([1 if r['label'] == 'positive' else 0 for r in sub])
            _a, scores = loo_logistic(Xb, yb, return_scores=True)
            base = {r['job_name']: s for r, s in zip(sub, scores)}
    except Exception as e:                                    # never block the AF3 numbers
        print(f'(composition baseline unavailable: {e})\n')

    # ------------------------------------------------ pre-registered, primary only
    print(f'=== pre-registered, primary metric = {PRIMARY} ===')
    print('   comparator is the composition baseline, not 0.5 - see delta column\n')
    sign = dict(METRICS)[PRIMARY]
    for name, (pp, nn) in SUBSETS.items():
        pos, neg = vals(PRIMARY, pp, sign), vals(PRIMARY, nn, sign)
        if not pos or not neg:
            print(f'{name:44s} not computable'); continue
        a = auc(pos, neg); lo, hi = bootstrap_auc(pos, neg)
        line = f'{name:44s} AUC {a:.3f}  CI [{lo:.3f}, {hi:.3f}]  n={len(pos)}/{len(neg)}'
        bp = [base[r['job_name']] for r in rows if pp(r) and r['job_name'] in base]
        bn = [base[r['job_name']] for r in rows if nn(r) and r['job_name'] in base]
        if len(bp) == len(pos) and len(bn) == len(neg) and bp and bn:
            d, dlo, dhi = delta_auc_vs_baseline(pos, neg, bp, bn)
            ok = dlo > 0
            line += (f'\n{"":44s} baseline {auc(bp, bn):.3f}   delta {d:+.3f} '
                     f'CI [{dlo:+.3f}, {dhi:+.3f}]' + ('  PASS' if ok else '  FAIL'))
        print(line)

    # ------------------------------------------------ secondary panel, BH-corrected
    print('\n=== secondary metric panel (BH-corrected across metrics, per subset) ===')
    print('   ipTM is expected to be floored for fuzzy pairs; these are the metrics')
    print('   that can still separate. Corrected because the panel is 16 wide.\n')
    for name, (pp, nn) in SUBSETS.items():
        out = []
        for metric, sgn in METRICS:
            pos, neg = vals(metric, pp, sgn), vals(metric, nn, sgn)
            if len(pos) < 3 or len(neg) < 3:
                continue
            out.append((metric, auc(pos, neg), perm_p(pos, neg), len(pos), len(neg)))
        if not out:
            continue
        adj = benjamini_hochberg([o[2] for o in out])
        print(f'  {name}')
        for (metric, a, p, n1, n2), q in sorted(zip(out, adj), key=lambda t: -abs(t[0][1] - 0.5)):
            star = ' *' if q < 0.05 else ''
            print(f'    {metric:28s} AUC {a:.3f}   p={p:.4f}  q={q:.4f}  n={n1}/{n2}{star}')
        print()

    # ------------------------------------------------ confounds
    print('=== confound checks ===')
    L = [float(r['total_residues']) for r in rows]
    for metric, _s in METRICS:
        v = [(float(r[metric]), float(r['total_residues']))
             for r in rows if r.get(metric) not in ('', None)]
        if len(v) < 5:
            continue
        rho = spearman([a for a, _ in v], [b for _, b in v])
        if abs(rho) > 0.5:
            print(f'  {metric:28s} rho(metric, total residues) = {rho:+.3f}'
                  f'   <-- length-confounded, report partialled')

    print('\n=== inter-seed spread of ipTM (fuzzy should exceed folded) ===')
    for mode in ('fuzzy', 'mutual_folding'):
        sds = [float(r['iptm_sd']) for r in rows
               if r.get('iptm_sd') not in ('', None)
               and r['binding_mode'] == mode and r['label'] == 'positive']
        if sds:
            print(f'  {mode:16s} mean sd {statistics.mean(sds):.4f}  (n={len(sds)})')


if __name__ == '__main__':
    argv = [a for a in sys.argv[1:] if not a.startswith('--')]
    if not argv:
        print(__doc__); sys.exit(1)
    main(*argv[:2], want_full='--no-full-data' not in sys.argv)
