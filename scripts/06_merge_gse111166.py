#!/usr/bin/env python3
"""HPICDA Stage 2 — GSE111166 (Han 2018, PMID 29642015) merge, memory-safe.

CRITICAL BUG FIX: PREa/POSTa files and BM-* files use DIFFERENT UID systems
for the same guide sequence (e.g. SLC35A1 guide ATAGCGCCAAACCCTAATAA =
HGLibA_52234 in PREa but HGLibA_44947 in BM files). Merge key MUST be the
guide sequence, never the UID.

Design (verified against paper PMC5939577 + ORCS screen pages):
  sequential screen (expansion between rounds): Rd1..Rd5, 2 reps each
    BM-17/18=Rd1, BM-13/14=Rd2, BM-9/10=Rd3, BM-5/6=Rd4, BM-1/2=Rd5
  preliminary consecutive screen (no expansion): POSTa = Prelim Rd5
  PREa = pre-selection A549-GeCKO LibA library (reference)
  ORCS-1148 = Rd2 vs PRE (798 hits p<0.05)
  ORCS-1149 = Rd5 vs PRE (501 hits p<0.05)
  ORCS-1150 = Prelim Rd5 (POSTa) vs PRE (119 hits p<0.1)
"""
import csv, gzip, os, sys, zipfile, re
import xml.etree.ElementTree as ET

BASE = "/workspace/HPICDA/03_functional_screens/raw/GSE111166"

def read_gz_csv(path, col_names, seq_idx, cnt_idx):
    d = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split(",")
            if len(p) <= max(seq_idx, cnt_idx):
                continue
            d[p[seq_idx]] = p[cnt_idx]
    return d

def read_bm(path):
    # header: lib_index, sgRNA, Gene, BM-N  (tab or comma? inspect)
    d = {}
    with gzip.open(path, "rt") as f:
        header = f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                p = line.rstrip("\n").split(",")
            d[p[1]] = p[3]
    return d

def read_xlsx_counts(path, val_col_name):
    """Stream-parse xlsx sheet1: columns seq, gene_id, UID, <val_col_name>."""
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for ev, el in ET.iterparse(z.open("xl/sharedStrings.xml"), events=("end",)):
            if el.tag == ns + "si":
                shared.append("".join(t.text or "" for t in el.iter(ns + "t")))
                el.clear()
    d = {}
    header = None
    for ev, row in ET.iterparse(z.open("xl/worksheets/sheet1.xml"), events=("end",)):
        if row.tag != ns + "row":
            continue
        vals = []
        for c in row.iter(ns + "c"):
            v = c.find(ns + "v")
            txt = v.text if v is not None else ""
            if c.get("t") == "s" and txt != "":
                txt = shared[int(txt)]
            vals.append(txt)
        if header is None:
            header = vals
        else:
            if len(vals) >= 4:
                d[vals[0]] = vals[3]
        row.clear()
    return d

def main():
    pre = read_gz_csv(os.path.join(BASE, "GSM3024555_LibA_sgRNA_PREa.csv.gz"),
                      None, 0, 3)
    # gene map from PREa
    gene_of = {}
    with gzip.open(os.path.join(BASE, "GSM3024555_LibA_sgRNA_PREa.csv.gz"), "rt") as f:
        for line in f:
            p = line.rstrip("\n").split(",")
            if len(p) >= 4 and p[2].startswith("HGLibA"):
                gene_of[p[0]] = p[1]
    bm = {
        "Rd1_rep1": "GSM3024552_LibA_sgRNA_BM-17.txt.gz", "Rd1_rep2": "GSM3024553_LibA_sgRNA_BM-18.txt.gz",
        "Rd2_rep1": "GSM3024550_LibA_sgRNA_BM-13.txt.gz", "Rd2_rep2": "GSM3024551_LibA_sgRNA_BM-14.txt.gz",
        "Rd3_rep1": "GSM3024548_LibA_sgRNA_BM-9.txt.gz",  "Rd3_rep2": "GSM3024549_LibA_sgRNA_BM-10.txt.gz",
        "Rd4_rep1": "GSM3024546_LibA_sgRNA_BM-5.txt.gz",  "Rd4_rep2": "GSM3024547_LibA_sgRNA_BM-6.txt.gz",
        "Rd5_rep1": "GSM3024544_LibA_sgRNA_BM-1.txt.gz",  "Rd5_rep2": "GSM3024545_LibA_sgRNA_BM-2.txt.gz",
    }
    counts = {name: read_bm(os.path.join(BASE, fn)) for name, fn in bm.items()}
    posta = read_xlsx_counts(os.path.join(BASE, "GSM3024554_LibA_sgRNA_POSTa.xlsx"), "POSTa")
    counts["POSTa"] = posta

    cols = ["PREa"] + [f"Rd{r}_rep{rep}" for r in range(1, 6) for rep in (1, 2)] + ["POSTa"]
    all_counts = {"PREa": pre, **counts}
    out = os.path.join(BASE, "counts_mageck.txt")
    n0 = 0
    with open(out, "w") as o:
        o.write("sgRNA\tGene\t" + "\t".join(cols) + "\n")
        for seq, gene in gene_of.items():
            row = [all_counts[c].get(seq, "0") for c in cols]
            o.write(seq + "\t" + gene + "\t" + "\t".join(row) + "\n")
            n0 += 1
    # sanity
    print("guides written:", n0)
    for g in ["SLC35A1", "CIC"]:
        with open(out) as f:
            for line in f:
                if "\t" + g + "\t" in line:
                    print(line.rstrip("\n"))

if __name__ == "__main__":
    main()
