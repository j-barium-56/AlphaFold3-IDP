#!/usr/bin/env python3
"""
Build a static GitHub Pages site for the AF3 IDR viewer, equivalent to what
viewer/serve.py serves live, so nothing needs to run at request time.

    python3 viewer/build_static.py [results_dir]      # writes ./docs

Reuses serve.py's load() for the job/model matching so there is one place that
decides which folder belongs to which job. Each full_data_N.json is trimmed to
the four fields the viewer actually reads (pae, contact_probs, token_chain_ids,
token_res_ids) with floats rounded, which is most of the size win: PAE/contact
matrices compress ~10x under gzip, and GitHub Pages serves gzip automatically.
"""
import json, os, shutil, sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)

from serve import load                             # noqa: E402  (path set up above)
import af3_workbook_build as W                      # noqa: E402

MATRIX_FIELDS = ('pae', 'contact_probs', 'token_chain_ids', 'token_res_ids')
ROUND = {'pae': 2, 'contact_probs': 3}


def round_nested(x, nd):
    if isinstance(x, float):
        return round(x, nd)
    if isinstance(x, list):
        return [round_nested(v, nd) for v in x]
    return x


def trim(path):
    d = W.load_json(path)
    out = {k: d[k] for k in MATRIX_FIELDS}
    for k, nd in ROUND.items():
        out[k] = round_nested(out[k], nd)
    return out


def build(results_dir, out_dir):
    jobs, matrices = load(results_dir)
    matrix_root = os.path.join(out_dir, 'api', 'matrix')
    n = 0
    for job_id, models in matrices.items():
        if not models:
            continue
        job_dir = os.path.join(matrix_root, job_id)
        os.makedirs(job_dir, exist_ok=True)
        for model, path in models.items():
            with open(os.path.join(job_dir, f'{model}.json'), 'w') as f:
                json.dump(trim(path), f, separators=(',', ':'))
            n += 1
    os.makedirs(os.path.join(out_dir, 'api'), exist_ok=True)
    with open(os.path.join(out_dir, 'api', 'jobs.json'), 'w') as f:
        json.dump(jobs, f, separators=(',', ':'))
    shutil.copyfile(os.path.join(HERE, 'index.html'), os.path.join(out_dir, 'index.html'))
    print(f'{sum(bool(j["samples"]) for j in jobs)} of {len(jobs)} jobs have results; '
          f'wrote {n} matrix files to {out_dir}')


if __name__ == '__main__':
    results_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, 'data')
    build(results_dir, os.path.join(REPO, 'docs'))
