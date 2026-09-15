#!/usr/bin/env python3
"""
A very small local viewer for the AF3 results: the job table and, for each server
sample, the PAE and contact-probability matrices.

    python3 viewer/serve.py [results_dir] [--port 8765]      # then open the printed URL

Reads af3_idr_pairs_1.csv, af3_scored.csv and the result folders (default: data/) at
start-up and hands each full_data_N.json to the page when it is opened, so there is
nothing to build and nothing to go stale. Listens on 127.0.0.1 only.
"""
import argparse, json, os, re, sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import unquote, urlparse

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import af3_workbook_build as W                    # same job <-> folder matching as the workbook


def load(results_dir):
    """(jobs for the page, {job_id: {sample: full_data path}})"""
    scored = {r['job_id']: r for r in W.read_csv('af3_scored.csv')}
    dirs = W.find_job_dirs(results_dir) if os.path.isdir(results_dir) else {}
    jobs, matrices = [], {}
    for p in W.read_csv('af3_idr_pairs_1.csv'):
        job_dir = W.match_dir(p['job_name'], dirs)
        samples, full = [], {}
        if job_dir:
            full = W.model_files(job_dir, 'full_data')
            for i, path in W.model_files(job_dir, 'summary_confidences').items():
                s = W.load_json(path)
                s.pop('chain_ids', None)          # one entry per token; the matrix file has them
                samples.append({'model': i, 'has_matrix': i in full, **s})
        matrices[p['job_id']] = full
        top = W.top_model([s['ranking_score'] for s in samples]) if samples else None
        jobs.append({**p, 'scored': scored.get(p['job_id']), 'samples': samples,
                     'top_model': None if top is None else samples[top]['model'],
                     'memorised_control': p['job_id'] in W.MEM})
    return jobs, matrices


def handler(jobs_body, matrices):
    index = os.path.join(HERE, 'index.html')

    class Handler(BaseHTTPRequestHandler):
        def reply(self, code, body, ctype):
            self.send_response(code)
            self.send_header('Content-Type', ctype)
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            url = urlparse(self.path)
            if url.path in ('/', '/index.html'):
                with open(index, 'rb') as f:          # read per request, so edits show on reload
                    return self.reply(200, f.read(), 'text/html; charset=utf-8')
            if url.path == '/api/jobs.json':
                return self.reply(200, jobs_body, 'application/json')
            m = re.fullmatch(r'/api/matrix/([^/]+)/(\d+)\.json', url.path)   # same layout build_static.py writes
            if m:
                try:
                    path = matrices[unquote(m.group(1))][int(m.group(2))]
                except (KeyError, ValueError):
                    return self.reply(404, b'no such job or sample', 'text/plain')
                with open(path, 'rb') as f:
                    return self.reply(200, f.read(), 'application/json')
            self.reply(404, b'not found', 'text/plain')

        def log_message(self, *_args):                # keep the terminal quiet
            pass

    return Handler


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Local PAE / contact-map viewer for the AF3 results.')
    ap.add_argument('results_dir', nargs='?', default=os.path.join(REPO, 'data'))
    ap.add_argument('--port', type=int, default=8765)
    args = ap.parse_args()
    jobs, matrices = load(args.results_dir)
    print(f'{sum(bool(j["samples"]) for j in jobs)} of {len(jobs)} jobs have results in {args.results_dir}')
    server = ThreadingHTTPServer(('127.0.0.1', args.port), handler(json.dumps(jobs).encode(), matrices))
    print(f'open http://127.0.0.1:{args.port}   (ctrl-c stops it)')
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
