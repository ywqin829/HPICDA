#!/usr/bin/env python3
"""Fetch ORCS curated gene-level score tables for Track B screens.

For each ORCS screen ID: GET the screen page (session + hidden fields),
then POST to scripts/datatableTools.php (tool=serverSideRows, type=scores)
to retrieve the full curated score table. Saves raw JSON + parsed TSV.

Memory-safe, polite (1.5 s delay), resumable (skips existing TSV).
"""
import json, re, sys, time, urllib.request, urllib.parse, http.cookiejar
from pathlib import Path

UA = 'Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0'
BASE = 'https://orcs.thebiogrid.org'
OUT = Path('/workspace/HPICDA/03_functional_screens/trackB_orcs')
OUT.mkdir(parents=True, exist_ok=True)

SCREEN_IDS = ['2067']

cj = http.cookiejar.CookieJar()
opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
opener.addheaders = [('User-Agent', UA)]

def get_hidden(html, name, default=''):
    m = re.search(r"id='%s' value='([^']*)'" % name, html)
    return m.group(1) if m else default

def fetch_screen(sid):
    # 1) screen page for session + hidden fields
    html = opener.open(f'{BASE}/Screen/{sid}', timeout=60).read().decode('utf-8', 'replace')
    sig_ind = get_hidden(html, 'significanceIndicator')
    sig_det = get_hidden(html, 'significanceDetails').replace('&gt;', '>').replace('&lt;', '<')
    score_cols = get_hidden(html, 'scoreCols', '1')
    m = re.search(r'Full Dataset Size</strong>: ([0-9,]+)', html)
    full_size = m.group(1) if m else ''
    m = re.search(r'Number of Hits</strong>: ([0-9,]+)', html)
    n_hits = m.group(1) if m else ''
    time.sleep(1.5)
    # 2) header
    hd = {'tool': 'serverSideHeader', 'type': 'scores', 'screenID': sid,
          'isHit': 'true', 'significanceIndicator': sig_ind,
          'significanceDetails': sig_det, 'scoreCols': score_cols}
    hd['expData'] = json.dumps(hd)
    req = urllib.request.Request(BASE + '/scripts/datatableTools.php',
                                 data=urllib.parse.urlencode({'expData': hd['expData']}).encode())
    header = json.loads(opener.open(req, timeout=120).read().decode())
    time.sleep(1.5)
    # 3) full table
    d = {'draw': 1, 'start': 0, 'length': 40000, 'tool': 'serverSideRows',
         'type': 'scores', 'screenID': sid, 'isHit': 'true',
         'significanceIndicator': sig_ind, 'significanceDetails': sig_det,
         'scoreCols': score_cols, 'totalRecords': '', 'checkedBoxes': '{}'}
    d['expData'] = json.dumps(d)
    req = urllib.request.Request(BASE + '/scripts/datatableTools.php',
                                 data=urllib.parse.urlencode({'expData': d['expData']}).encode())
    out = json.loads(opener.open(req, timeout=180).read().decode())
    rows = out.get('data', [])
    # parse
    col_names = []
    if isinstance(header, list):
        for c in header:
            t = c.get('title', '') if isinstance(c, dict) else str(c)
            col_names.append(re.sub(r'<[^>]+>', '', t).strip())
    parsed = []
    for r in rows:
        gm = re.search(r"title='([^']+)'", r[0])
        gene = gm.group(1) if gm else re.sub(r'<[^>]+>', '', r[0])
        scores = [x for x in r[4:] if not str(x).startswith('<div')]
        scores = [re.sub(r'<[^>]+>', '', str(s)).strip() for s in scores]
        parsed.append([gene] + scores)
    n_score_cols = len(parsed[0]) - 1 if parsed else 0
    score_headers = col_names[4:4 + n_score_cols] if len(col_names) >= 4 + n_score_cols else [f'score{i+1}' for i in range(n_score_cols)]
    tsv = OUT / f'ORCS-{sid}.gene_scores.tsv'
    with open(tsv, 'w') as f:
        f.write('gene\t' + '\t'.join(score_headers) + '\n')
        for p in parsed:
            f.write('\t'.join(p) + '\n')
    meta = {'screen_id': sid, 'full_dataset_size': full_size, 'n_hits': n_hits,
            'n_curated_rows': len(parsed), 'score_headers': score_headers,
            'significance': sig_det}
    with open(OUT / f'ORCS-{sid}.meta.json', 'w') as f:
        json.dump(meta, f, indent=2)
    return meta

results = []
for sid in SCREEN_IDS:
    tsv = OUT / f'ORCS-{sid}.gene_scores.tsv'
    if tsv.exists() and tsv.stat().st_size > 100:
        print(f'SKIP {sid} (exists)', flush=True)
        continue
    try:
        meta = fetch_screen(sid)
        print(f"OK {sid}: {meta['n_curated_rows']} rows, hits={meta['n_hits']}, full={meta['full_dataset_size']}, cols={meta['score_headers']}", flush=True)
        results.append(meta)
        time.sleep(1.5)
    except Exception as e:
        print(f'FAIL {sid}: {e}', flush=True)
        results.append({'screen_id': sid, 'error': str(e)})

with open(OUT / 'fetch_summary.json', 'w') as f:
    json.dump(results, f, indent=2)
print('DONE', len(results))
