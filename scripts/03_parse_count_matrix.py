#!/usr/bin/env python3
"""HPICDA Stage 2 — Track A count-matrix parser (production rule).

Parses a CRISPR screen raw count matrix into MAGeCK-ready format
(sgRNA<TAB>Gene<TAB>sample counts...).

Parsing rules (documented per blueprint provenance discipline):
  R1. Guide/gene separator: last underscore + digits = global construct ID
      (verified for GSE210689: number is a global construct index, NOT a
      per-gene guide index; e.g. A1BG_1..4, GPR156_26222..26225).
  R2. Multi-mapping guides (comma-separated Construct IDs, e.g.
      PRR20A_53067,PRR20B_53071,...): DROPPED. A guide that cannot be
      uniquely attributed to one gene must not contribute to any single
      gene's RRA evidence. Counted and reported.
  R3. NO_CURRENT_* and MATCH_LIMIT_EXCEEDED_* (unannotated guides):
      DROPPED, counted and reported.
  R4. Gene symbols are NOT case-normalized; library annotation is used
      as-is and harmonized at Stage 5 against the ontology/DepMap universe.

Usage:
  python3 03_parse_count_matrix.py --input RAW.txt[.gz] \
      --id-col "Construct IDs" --barcode-col "Construct Barcode" \
      --output counts_mageck.txt --report parse_report.json
"""
import argparse, gzip, json, re, sys
from collections import Counter

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--id-col", default="Construct IDs")
    ap.add_argument("--barcode-col", default="Construct Barcode")
    ap.add_argument("--output", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    pat = re.compile(r"^([A-Za-z0-9.\-]+)_(\d+)$")
    n_total = 0
    n_kept = 0
    n_multi = 0
    n_unannotated = 0
    n_malformed = 0
    genes = Counter()

    with open_maybe_gz(args.input) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            id_idx = header.index(args.id_col)
            bc_idx = header.index(args.barcode_col)
        except ValueError:
            sys.exit(f"ERROR: required columns not found in header: {header}")
        sample_cols = [i for i in range(len(header)) if i not in (id_idx, bc_idx)]
        sample_names = [header[i] for i in sample_cols]

        with open(args.output, "w") as out:
            out.write("sgRNA\tGene\t" + "\t".join(sample_names) + "\n")
            for line in f:
                n_total += 1
                parts = line.rstrip("\n").split("\t")
                cid = parts[id_idx]
                bc = parts[bc_idx]
                if cid.startswith("NO_CURRENT") or cid.startswith("MATCH_LIMIT_EXCEEDED"):
                    n_unannotated += 1
                    continue
                if "," in cid:
                    n_multi += 1
                    continue
                m = pat.match(cid)
                if not m:
                    n_malformed += 1
                    continue
                gene = m.group(1)
                counts = [parts[i] for i in sample_cols]
                out.write(bc + "\t" + gene + "\t" + "\t".join(counts) + "\n")
                genes[gene] += 1
                n_kept += 1

    report = {
        "input": args.input,
        "output": args.output,
        "total_guides": n_total,
        "kept_guides": n_kept,
        "dropped_multimapping": n_multi,
        "dropped_unannotated": n_unannotated,
        "dropped_malformed": n_malformed,
        "n_genes": len(genes),
        "median_guides_per_gene": sorted(genes.values())[len(genes)//2] if genes else 0,
        "samples": sample_names,
    }
    with open(args.report, "w") as r:
        json.dump(report, r, indent=2)
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
