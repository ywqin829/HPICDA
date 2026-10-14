#!/usr/bin/env python3
"""HPICDA Stage 2 — GSE175666 (Rebendenne 2022 Nat Genet, PMID 35879413).

Track A from deposited per-guide reads (genomewide_reads sheets; float
normalized values rounded to integers — documented approximation).
  ORCS-1862 Calu-3 KO  (Hsu lib): infected D4/D10 vs D0 controls
  ORCS-1863 Calu-3 Act (Hsu lib): D5 infected vs D0 controls (CRISPRa;
        direction sign-inverted to LOF at Stage 5 per Stage 1 adjudication)
  ORCS-1864 Caco-2-ACE2 KO (Brunello): resistant vs initial, screens 1+2.1
        (gene via guide_mapping sheet join)
"""
import csv, os, subprocess, sys

BASE = "/workspace/HPICDA/03_functional_screens/raw/GSE175666"

def rnd(v):
    try:
        return str(max(0, int(round(float(v)))))
    except ValueError:
        return "0"

def build(src, dst, treat, ctrl, guide_col=0, gene_col=1, map_file=None):
    gmap = None
    if map_file:
        gmap = {}
        with open(map_file) as f:
            r = csv.reader(f, delimiter="\t")
            next(r)
            for row in r:
                if len(row) >= 2:
                    gmap[row[0]] = row[1]
    with open(src) as f:
        header = f.readline().rstrip("\n").split("\t")
        ti = [header.index(c) for c in treat]
        ci = [header.index(c) for c in ctrl]
        with open(dst, "w") as o:
            o.write("sgRNA\tGene\t" + "\t".join(treat + ctrl) + "\n")
            n = 0
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < len(header):
                    continue
                guide = p[guide_col]
                gene = gmap.get(guide, "") if gmap else p[gene_col]
                if not gene:
                    continue
                o.write(guide + "\t" + gene + "\t" +
                        "\t".join(rnd(p[i]) for i in ti + ci) + "\n")
                n += 1
    return n

env = dict(os.environ)
env["PATH"] = "/workspace/.conda/envs/crispr/bin:" + env["PATH"]

jobs = [
    ("ORCS-1862", f"{BASE}/GSE175666_Calu3_KO_Hsu_v1_reads.tsv", None,
     ["E-SS_D4_MOI01_r1_run1", "Hv3_4_av-SS_R2_D4_run3", "Iv4-And_R2_D10_run3"],
     ["B-SS_D0_r1_run1", "D-SS_D0_R2_run1", "F_av-And_R2_D0_Pt2_run2"]),
    ("ORCS-1863", f"{BASE}/GSE175666_Calu3_Act_Hsu_v1_reads.tsv", None,
     ["L-SAM_D5_r1", "M-SAM_D5_r2"],
     ["J-SAM_D0_r1", "K-SAM_D0_r2"]),
    ("ORCS-1864", f"{BASE}/GSE175666_Supplementary_Data_5_Caco-2-ACE2_KO_Brunello_Goujon_reads.tsv",
     f"{BASE}/caco2_guide_map.tsv",
     ["Caco2 screen 1 resistant population", "Caco2 screen 2.1 resistant population"],
     ["Caco2 screen 1 initial population", "Caco2 screen 2.1 initial population"]),
]
for sid, src, mapf, treat, ctrl in jobs:
    cnt = os.path.join(BASE, sid + "_counts.txt")
    n = build(src, cnt, treat, ctrl, map_file=mapf)
    print(sid, "guides:", n, file=sys.stderr)
    cmd = ["mageck", "test", "-k", cnt, "-t", ",".join(treat), "-c", ",".join(ctrl),
           "-n", os.path.join(BASE, sid), "--norm-method", "median"]
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-1500:], file=sys.stderr); print(res.stderr[-1500:], file=sys.stderr)
        sys.exit(f"mageck failed for {sid}")
print("GSE175666_DONE", file=sys.stderr)
