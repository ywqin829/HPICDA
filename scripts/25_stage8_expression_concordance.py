#!/usr/bin/env python3
"""Stage 8: candidate-gene expression concordance in human infection datasets.

Per human_layer_rules.yaml: pre-registered candidate set only; direction
consistency evidence; upgrade-only; no de novo discovery.

Usage:
  python3 25_stage8_expression_concordance.py inspect   # print group sizes only
  python3 25_stage8_expression_concordance.py run       # full analysis
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path("/workspace/HPICDA")
RAW = ROOT / "09_transcriptomic_context/raw"
OUT = ROOT / "09_transcriptomic_context"

# ---------------- metadata parsing ----------------
def read_series_matrix_meta(path):
    """Return DataFrame: gsm, title, chars (list)."""
    with gzip.open(path, "rt", errors="replace") as f:
        lines = [l.rstrip("\n") for l in f if l.startswith("!Sample_")]
    d = {}
    for l in lines:
        parts = l.split("\t")
        key = parts[0]
        vals = [p.strip('"') for p in parts[1:]]
        d.setdefault(key, []).append(vals)
    gsm = d.get("!Sample_geo_accession", [[]])[0]
    title = d.get("!Sample_title", [[]])[0]
    chars = d.get("!Sample_characteristics_ch1", [])
    char_df = {}
    for row in chars:
        for i, cell in enumerate(row):
            if ":" in cell:
                k, v = cell.split(":", 1)
                char_df.setdefault(k.strip().lower(), {})[i] = v.strip()
    m = pd.DataFrame({"gsm": gsm, "title": title})
    for k, v in char_df.items():
        m[k] = pd.Series(v)
    return m

def groups_from(meta, case_fn, ctrl_fn):
    case = meta[meta.apply(case_fn, axis=1)]["gsm"].tolist()
    ctrl = meta[meta.apply(ctrl_fn, axis=1)]["gsm"].tolist()
    return case, ctrl

# ---------------- expression loading ----------------
def load_series_matrix_expr(path):
    with gzip.open(path, "rt", errors="replace") as f:
        lines = f.readlines()
    hdr = None
    for i, l in enumerate(lines):
        if l.startswith("!series_matrix_table_begin"):
            hdr = i + 1
            break
    df = pd.read_csv(path, sep="\t", skiprows=hdr, compression="gzip",
                     index_col=0, comment="!")
    df = df[df.index.notna()]
    return df

def load_table(path, sep=None):
    kw = dict(compression="gzip", index_col=0)
    if sep and len(sep) > 1:  # regex separator needs python engine
        kw["engine"] = "python"
    else:
        kw["low_memory"] = False
    df = pd.read_csv(path, sep=sep or "\t", **kw)
    df.columns = [str(c).strip('"') for c in df.columns]
    return df

# ---------------- per-dataset configs ----------------
# group functions operate on metadata row
CFG = {
    "GSE157103": dict(col_style="covid_nc", pathogen="SARS-CoV-2", expr="GSE157103_tpm.tsv.gz", kind="tpm",
        meta="GSE157103_series_matrix.txt.gz",
        case=lambda r: str(r["title"]).startswith("COVID"),
        ctrl=lambda r: str(r["title"]).startswith("NONCOVID")),
    "GSE171110": dict(col_style="title_prefix", pathogen="SARS-CoV-2", expr="GSE171110_rawcounts.txt.gz", kind="counts",
        meta="GSE171110_series_matrix.txt.gz",
        case=lambda r: "covid" in str(r.get("title", "")).lower() and "healthy" not in str(r.get("title", "")).lower(),
        ctrl=lambda r: "healthy" in str(r.get("title", "")).lower() or "control" in str(r.get("title", "")).lower()),
    "GSE155635": dict( file_sep=" ", symbol_col="gene_symbol", col_style="title", pathogen="Influenza", expr="GSE155635_norm.txt.gz", kind="norm",
        meta="GSE155635_series_matrix.txt.gz",
        case=lambda r: str(r.get("status", "")).lower() == "infected",
        ctrl=lambda r: str(r.get("status", "")).lower() == "healthy_control"),
    "GSE158592": dict( file_sep=" ", symbol_col="gene_symbol", col_style="title_first", pathogen="Influenza", expr="GSE158592_norm.txt.gz", kind="norm",
        meta="GSE158592_series_matrix.txt.gz",
        case=lambda r: str(r.get("status", "")).lower() == "infected",
        ctrl=lambda r: str(r.get("status", "")).lower() == "healthy_control"),
    "GSE101702": dict( probe_gpl="GPL21185", pathogen="Influenza", expr="GSE101702_series_matrix.txt.gz", kind="array",
        meta="GSE101702_series_matrix.txt.gz",
        case=lambda r: str(r.get("diagnosis", "")).lower() == "influenza",
        ctrl=lambda r: str(r.get("diagnosis", "")).lower() == "healthy"),
    "GSE38900": dict( probe_gpl="GPL10558", pathogen="RSV", expr="GSE38900-GPL10558_series_matrix.txt.gz", kind="array",
        meta="GSE38900-GPL10558_series_matrix.txt.gz",
        case=lambda r: bool(re.search(r"\bRSV\b", str(r.get("title", "")), re.I)),
        ctrl=lambda r: bool(re.search(r"healthy control", str(r.get("title", "")), re.I))),
    "GSE77087": dict( probe_gpl="GPL10558", pathogen="RSV", expr="GSE77087_series_matrix.txt.gz", kind="array",
        meta="GSE77087_series_matrix.txt.gz",
        # disease==RSV covers Inpatient+Outpatient (81); contrast = RSV vs Healthy
        case=lambda r: str(r.get("disease", "")).lower() == "rsv",
        ctrl=lambda r: str(r.get("disease", "")).lower() == "healthy"),
    "GSE83456": dict( probe_gpl="GPL10558", pathogen="Mtb", expr="GSE83456_series_matrix.txt.gz", kind="array",
        meta="GSE83456_series_matrix.txt.gz",
        # primary contrast PTB vs HC; Sarcoid/EPTB excluded per selection audit
        case=lambda r: str(r.get("disease state", "")).lower() == "ptb",
        ctrl=lambda r: str(r.get("disease state", "")).lower() == "hc"),
    "GSE94438": dict( symbol_col="symbol", col_style="code", pathogen="Mtb", expr="GSE94438_rawcounts.csv.gz", kind="counts_csv",
        meta="GSE94438_series_matrix.txt.gz",
        case=lambda r: str(r.get("group", "")).lower() == "case (tb)",
        ctrl=lambda r: str(r.get("group", "")).lower() == "control"),
}

def normalize_expr(df, kind):
    df = df.apply(pd.to_numeric, errors="coerce")
    if kind in ("tpm", "norm"):
        x = np.log2(df + 1.0)
    elif kind in ("counts", "counts_csv"):
        lib = df.sum(axis=0)
        cpm = df.div(lib, axis=1) * 1e6
        x = np.log2(cpm + 1.0)
    else:  # array: check if already log
        mx = np.nanmax(df.to_numpy().ravel())
        x = df if mx < 100 else np.log2(df.clip(lower=0.1))
    return x

# ---------------- probe->symbol maps (array platforms) ----------------
_GPL_CACHE = {}
def gpl_map(gpl):
    if gpl in _GPL_CACHE:
        return _GPL_CACHE[gpl]
    if gpl == "GPL10558":
        t = pd.read_csv(RAW / "GPL10558.annot.gz", sep="\t", skiprows=28, low_memory=False)
        t = t[~t["ID"].astype(str).str.startswith("!")]
        m = dict(zip(t["ID"], t["Gene symbol"].fillna("")))
    elif gpl == "GPL21185":
        t = pd.read_csv(RAW / "GPL21185_table.tsv", sep="\t", skiprows=1, low_memory=False)
        t = t[~t["ID"].astype(str).str.startswith("!")]
        m = dict(zip(t["ID"], t["GENE_SYMBOL"].fillna("")))
    else:
        raise ValueError(gpl)
    _GPL_CACHE[gpl] = m
    return m

def map_columns_to_gsm(expr_cols, meta, style):
    """Return {expr_col: gsm}."""
    out = {}
    if style == "title":
        lut = dict(zip(meta["title"], meta["gsm"]))
        out = {c: lut[c] for c in expr_cols if c in lut}
    elif style == "title_first":
        lut = {}
        for t, g in zip(meta["title"], meta["gsm"]):
            lut.setdefault(str(t).split()[0], g)
        out = {c: lut[c] for c in expr_cols if c in lut}
    elif style == "title_prefix":
        # expr col is a prefix of the title, e.g. col "507-V" vs title "507-V Healthy"
        pairs = sorted(((str(t), g) for t, g in zip(meta["title"], meta["gsm"])),
                       key=lambda x: -len(x[0]))
        for c in expr_cols:
            for t, g in pairs:
                if t.startswith(c) or str(t).split()[0] == c:
                    out[c] = g
                    break
    elif style == "code":
        lut = {"X" + str(c): g for c, g in zip(meta["code"], meta["gsm"])}
        out = {c: lut[c] for c in expr_cols if c in lut}
    elif style == "covid_nc":
        # GSE157103: COVID_NN_... -> C{N}; NONCOVID_NN_... -> NC{N}
        lut = {}
        for t, g in zip(meta["title"], meta["gsm"]):
            m = re.match(r"(NONCOVID|COVID)_(\d+)_", str(t))
            if m:
                key = ("NC" if m.group(1) == "NONCOVID" else "C") + str(int(m.group(2)))
                lut[key] = g
        out = {c: lut[c] for c in expr_cols if c in lut}
    return out

def main(mode):
    cand_path = ROOT / "02_metadata/preregistered_candidate_set_v1.tsv"
    if mode == "inspect" and not cand_path.exists():
        cand_genes = set()  # group verification only; candidate set not yet frozen
    else:
        cand = pd.read_csv(cand_path, sep="\t")
        cand_genes = set(cand["gene"])
    results = []
    for gse, c in CFG.items():
        meta = read_series_matrix_meta(RAW / c["meta"])
        case, ctrl = groups_from(meta, c["case"], c["ctrl"])
        print(f"{gse} ({c['pathogen']}): case={len(case)} ctrl={len(ctrl)}")
        if mode == "inspect":
            tt = meta[meta['gsm'].isin(case)]['title'].head(3).tolist()
            tc = meta[meta['gsm'].isin(ctrl)]['title'].head(3).tolist()
            print(f"   case e.g. {tt}")
            print(f"   ctrl e.g. {tc}")
            continue
        if len(case) < 5 or len(ctrl) < 5:
            print(f"   SKIP: insufficient groups")
            continue
        expr = load_series_matrix_expr(RAW / c["expr"]) if c["kind"] == "array" \
            else load_table(RAW / c["expr"], sep=c.get("file_sep") or ("," if c["kind"] == "counts_csv" else "\t"))
        # row index -> gene symbols
        if c.get("symbol_col"):
            expr.index = expr[c["symbol_col"]].fillna("").astype(str).str.strip('"')
            expr = expr.drop(columns=[col for col in
                ["Row.names","ENSEMBL_gene_ID","ensembl_transcript_id","gene_symbol",
                 "entrezgene_id","uniprotswissprot","description","refseq_mrna",
                 "transcript_length","symbol"] if col in expr.columns])
        elif c.get("probe_gpl"):
            pm = gpl_map(c["probe_gpl"])
            expr.index = [str(pm.get(str(i), "")) for i in expr.index]
        # columns -> gsm
        if c.get("col_style"):
            cm = map_columns_to_gsm(list(expr.columns), meta, c["col_style"])
            expr = expr.rename(columns=cm)
            expr = expr.loc[:, [col for col in expr.columns if str(col).startswith("GSM")]]
        common = [g for g in case + ctrl if g in expr.columns]
        if len(common) < 10:
            print(f"   SKIP: only {len(common)} columns aligned to GSM")
            continue
        x = normalize_expr(expr[common], c["kind"])
        case_c = [g for g in case if g in x.columns]
        ctrl_c = [g for g in ctrl if g in x.columns]
        print(f"   expr matrix {x.shape}, case={len(case_c)} ctrl={len(ctrl_c)}")
        # gene index harmonization: uppercase symbol match
        x.index = [str(i).split("///")[0].split(";")[0].strip().upper() for i in x.index]
        x = x[(x.index != "") & (x.index != "NAN")]
        x = x[~x.index.duplicated(keep="first")]
        for gene in cand_genes:
            g = gene.upper()
            if g not in x.index:
                continue
            a = x.loc[g, case_c].dropna()
            b = x.loc[g, ctrl_c].dropna()
            if len(a) < 5 or len(b) < 5:
                continue
            lfc = a.mean() - b.mean()
            try:
                p = mannwhitneyu(a, b, alternative="two-sided").pvalue
            except ValueError:
                p = np.nan
            results.append({"gse": gse, "pathogen": c["pathogen"], "gene": gene,
                            "log2fc": lfc, "p": p,
                            "direction": "up" if lfc > 0 else "down",
                            "n_case": len(a), "n_ctrl": len(b)})
    if mode == "run":
        res = pd.DataFrame(results)
        res["fdr"] = res.groupby("gse")["p"].transform(
            lambda p: pd.Series(bh(p.to_numpy()), index=p.index))
        res.to_csv(OUT / "stage8_candidate_expression_concordance.tsv", sep="\t", index=False)
        print("\nwrote", len(res), "gene-dataset tests")
        print(res.groupby(["gse"])["gene"].nunique().to_string())

def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    q = np.empty(n)
    q[o] = np.clip(r, 0, 1)
    return q

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "inspect")
