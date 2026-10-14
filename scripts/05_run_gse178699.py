#!/usr/bin/env python3
"""HPICDA Stage 2 — GSE178699 (Chan 2022, PMID 36597481) Track A processing.

Builds MAGeCK input from GSE178699_covid_screens_readCounts.txt.gz and runs
the 3 discovery contrasts with raw counts in this series:
  ORCS-2066  Calu-3   T43.SARS(MOI 0.0005) vs T25.mock-control
  ORCS-2070  HEK293A_T T12.SARS(MOI 2)     vs T12.mock-control
  ORCS-2071  HuH-7    T15.SARS2            vs T14.mock-control

Excluded / not present (leakage + provenance discipline):
  VERO.*            non-human (blueprint exclusion)
  CACO2.*           not in the frozen 52-screen registry corpus
  ORCS-2067         Calu-3 0.05 MOI / 58-day arm NOT in this GEO series
                    (only 35 samples GSM5396228-262) -> Track B via ORCS
  ORCS-2068/2069    UM-UC-4 holdout arms NOT in this GEO series (verified)
Reference-arm decision (HG2): time-matched mock-expanded control is primary
reference (isolates infection-specific selection; T0 retained as sensitivity).
"""
import csv, gzip, subprocess, sys, os

BASE = "/workspace/HPICDA/03_functional_screens/raw/GSE178699"
SRC = os.path.join(BASE, "GSE178699_covid_screens_readCounts.txt.gz")

SCREENS = {
    "ORCS-2066": (["CALU3.T43A.lowMOI.SARS", "CALU3.T43B.lowMOI.SARS", "CALU3.T43C.lowMOI.SARS"],
                  ["CALU3.T25A", "CALU3.T25B", "CALU3.T25C"]),
    "ORCS-2070": (["HEK293A_T.T12A.MOI2", "HEK293A_T.T12B.MOI2", "HEK293A_T.T12C.MOI2"],
                  ["HEK293A_T.T12A", "HEK293A_T.T12B", "HEK293A_T.T12C"]),
    "ORCS-2071": (["HuH7.T15A.SARS2", "HuH7.T15B.SARS2", "HuH7.T15C.SARS2"],
                  ["HuH7.T14A", "HuH7.T14B", "HuH7.T14C"]),
}

# load full matrix once
with gzip.open(SRC, "rt") as f:
    r = csv.reader(f, delimiter="\t")
    header = next(r)
    rows = list(r)
print("samples in matrix:", header[3:], file=sys.stderr)

env = dict(os.environ)
env["PATH"] = "/workspace/.conda/envs/crispr/bin:" + env["PATH"]

for sid, (treat, ctrl) in SCREENS.items():
    cols = treat + ctrl
    idx = []
    for c in cols:
        if c not in header:
            sys.exit(f"MISSING COLUMN: {c}")
        idx.append(header.index(c))
    out = os.path.join(BASE, f"{sid}_counts.txt")
    with open(out, "w") as o:
        o.write("sgRNA\tGene\t" + "\t".join(cols) + "\n")
        for row in rows:
            o.write(row[1] + "\t" + row[2] + "\t" + "\t".join(row[i] for i in idx) + "\n")
    cmd = ["mageck", "test", "-k", out, "-t", ",".join(treat), "-c", ",".join(ctrl),
           "-n", os.path.join(BASE, sid), "--norm-method", "median"]
    print("RUN:", " ".join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-2000:], file=sys.stderr)
        print(res.stderr[-2000:], file=sys.stderr)
        sys.exit(f"mageck failed for {sid}")
print("ALL_SCREENS_DONE", file=sys.stderr)
