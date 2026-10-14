#!/usr/bin/env python3
"""HPICDA Stage 2 — Li 2020 (PMID 31919360, ORCS-1207) Track A processing.

Data: Edinburgh DataShare doi:10.7488/ds/7818. Files contain log2 RPM
per-guide values (NOT raw counts): log2(count*1e6/total + 1).
Transformation: pseudo-counts = round((2^l2rpm - 1) * 10) at assumed 10M
depth; MAGeCK median normalization is scale-robust. Documented approximation.

Screen type (erratum #5): FACS-HA viral-antigen screen (16 hpi), NOT survival.
  screen1 primary (a priori): P5_6 (low/very-low HA) vs P4 (mid HA)
  screen2 primary (a priori): P6 (very-low HA) vs P4 (mid HA)
  enriched in low-HA fraction = KO_REDUCES_PATHOGEN (dependency factor).
"""
import csv, math, os, subprocess, sys

BASE = "/workspace/HPICDA/03_functional_screens/raw/Li2020_Edinburgh"
S1 = os.path.join(BASE, "93635e29-12b5-45ee-812a-d5b9036e90c2.bin")  # screen1.txt
S2 = os.path.join(BASE, "9dde4e84-9827-4c4e-903c-dadeaf90dd0a.bin")  # screen2.txt

def to_pseudo(v):
    try:
        x = float(v)
    except ValueError:
        return "0"
    c = (2.0 ** x - 1.0) * 10.0
    return str(max(0, int(round(c))))

def convert(src, dst, treat_cols, ctrl_cols):
    with open(src) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: header.index(c) for c in treat_cols + ctrl_cols}
        with open(dst, "w") as o:
            o.write("sgRNA\tGene\t" + "\t".join(treat_cols + ctrl_cols) + "\n")
            n = 0
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < len(header):
                    continue
                o.write(p[0] + "\t" + p[1] + "\t" +
                        "\t".join(to_pseudo(p[idx[c]]) for c in treat_cols + ctrl_cols) + "\n")
                n += 1
    return n

env = dict(os.environ)
env["PATH"] = "/workspace/.conda/envs/crispr/bin:" + env["PATH"]

jobs = [
    ("ORCS-1207-S1", S1, ["P5_6"], ["P4"]),
    ("ORCS-1207-S2", S2, ["P6"], ["P4"]),
]
for sid, src, treat, ctrl in jobs:
    cnt = os.path.join(BASE, sid + "_counts.txt")
    n = convert(src, cnt, treat, ctrl)
    print(sid, "guides:", n, file=sys.stderr)
    cmd = ["mageck", "test", "-k", cnt, "-t", ",".join(treat), "-c", ",".join(ctrl),
           "-n", os.path.join(BASE, sid), "--norm-method", "median"]
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-1500:], file=sys.stderr); print(res.stderr[-1500:], file=sys.stderr)
        sys.exit(f"mageck failed for {sid}")
print("LI2020_DONE", file=sys.stderr)
