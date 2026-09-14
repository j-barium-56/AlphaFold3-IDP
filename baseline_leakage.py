#!/usr/bin/env python3
"""
Trivial-baseline / leakage check for the IDR-IDR benchmark. Run this BEFORE AlphaFold.

If a sequence-composition feature already separates positives from negatives, then any
AUC AlphaFold gets is uninterpretable - it could be reading composition rather than
interaction. This computes, for every cheap feature, how well it alone separates the
labels, and does the same for a leave-one-out logistic model over all of them.

    python3 baseline_leakage.py [af3_idr_pairs.csv]

Rule of thumb: single-feature AUC above ~0.75, or LOO-logistic AUC above ~0.8, means the
set needs rebalancing before it can support a claim about AF3.
"""
import csv, sys
import numpy as np

KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'E':-3.5,'Q':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
POS, NEG, AROM = set('KR'), set('DE'), set('FWY')


def feats(s):
    n = len(s)
    f = lambda st: sum(c in st for c in s) / n
    fpos, fneg = f(POS), f(NEG)
    return {
        'len': n,
        'ncpr': fpos - fneg,
        'fcr': fpos + fneg,
        'kappa_proxy': abs(fpos - fneg) / (fpos + fneg) if (fpos + fneg) else 0.0,
        'aromatic': f(AROM),
        'gly_ser': f(set('GS')),
        'polar_qn': f(set('QN')),
        'pro': f(set('P')),
        'hydropathy': sum(KD.get(c, 0) for c in s) / n,
        'lowcomplexity': 1 - len(set(s)) / 20,
    }


def pair_feats(a, b):
    fa, fb = feats(a), feats(b)
    out = {}
    for k in fa:
        out[f'{k}_mean'] = (fa[k] + fb[k]) / 2
        out[f'{k}_absdiff'] = abs(fa[k] - fb[k])
    out['total_len'] = len(a) + len(b)
    out['len_ratio'] = max(len(a), len(b)) / min(len(a), len(b))
    out['charge_product'] = fa['ncpr'] * fb['ncpr']        # negative = complementary
    out['is_homotypic'] = 1.0 if a == b else 0.0
    return out


def auc(pos, neg):
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


def two_sided(a):
    """A feature that anti-correlates is just as much leakage as one that correlates."""
    return max(a, 1 - a)


def loo_logistic(X, y, return_scores=False):
    """Leave-one-out logistic regression AUC; plain gradient descent, no sklearn."""
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    X = np.hstack([X, np.ones((len(X), 1))])
    scores = np.zeros(len(y))
    for i in range(len(y)):
        m = np.ones(len(y), bool); m[i] = False
        Xt, yt = X[m], y[m]
        w = np.zeros(Xt.shape[1])
        for _ in range(2000):
            p = 1 / (1 + np.exp(-Xt @ w))
            w -= 0.05 * (Xt.T @ (p - yt) / len(yt) + 0.02 * w)   # L2, keeps it honest
        scores[i] = X[i] @ w
    a = auc(list(scores[y == 1]), list(scores[y == 0]))
    return (a, scores) if return_scores else a


def main(pairs_csv='af3_idr_pairs.csv'):
    rows = [r for r in csv.DictReader(open(pairs_csv)) if r['label'] != 'negative_scramble']
    y = np.array([1 if r['label'] == 'positive' else 0 for r in rows])
    F = [pair_feats(r['seq_A'], r['seq_B']) for r in rows]
    names = sorted(F[0])
    X = np.array([[f[k] for k in names] for f in F], float)

    def report(title, mask):
        yy, XX = y[mask], X[mask]
        if len(set(yy.tolist())) < 2 or min((yy == 1).sum(), (yy == 0).sum()) < 3:
            print(f'\n{title}: too few of one class to test'); return
        print(f'\n{title}  ({(yy==1).sum()} positive / {(yy==0).sum()} negative)')
        scored = sorted(
            ((n, two_sided(auc(list(XX[yy == 1, j]), list(XX[yy == 0, j]))))
             for j, n in enumerate(names)), key=lambda t: -t[1])
        for n, a in scored[:6]:
            flag = '   <-- LEAKAGE' if a > 0.75 else ''
            print(f'    {n:22s} AUC {a:.3f}{flag}')
        m = loo_logistic(XX, yy)
        flag = '   <-- LEAKAGE' if two_sided(m) > 0.80 else ''
        print(f'    {"[LOO logistic, all]":22s} AUC {two_sided(m):.3f}{flag}')

    allm = np.ones(len(y), bool)
    homo = np.array([r['pair_type'] == 'homotypic' for r in rows])
    print('Trivial sequence-only baseline. High AUC here means the AF3 result would be')
    print('uninterpretable, because composition alone already predicts the label.')
    report('ALL real-sequence pairs', allm)
    report('HETEROTYPIC only', ~homo)
    report('HOMOTYPIC only', homo)


if __name__ == '__main__':
    main(*sys.argv[1:2])
