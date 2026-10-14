#!/usr/bin/env python3
"""HPICDA Stage 5: per-screen conditional dependency calls (blueprint §14-15).

Harmonizes Track A (MAGeCK reanalysis) + Track B (ORCS author stats) into a
common LOF-directional gene-score layer, then assigns per-screen calls:
CDEP_PATHOGEN_EXPLOITED / CDEP_HOST_PROTECTIVE / BASELINE_FITNESS_CONFOUNDED /
TOXIN_SPECIFIC / DIRECTION_UNRESOLVED / NO_FUNCTIONAL_SIGNAL.

Direction semantics (HG1 + AM2-10):
  survival screens: enriched = KO_PROTECTS_HOST -> CDEP_PATHOGEN_EXPLOITED (INHIBIT)
  burden screens:   enriched = KO_REDUCES_PATHOGEN -> CDEP_PATHOGEN_EXPLOITED (INHIBIT)
  depleted (both)   -> CDEP_HOST_PROTECTIVE (PRESERVE/ACTIVATE)
  CRISPRa arms sign-inverted to LOF equivalents (ORCS-1175, 1706, 1708, 1863, 2020)
  Mtb PATHOGEN_STATE: signal recorded, intervention DIRECTION_UNRESOLVED (Stage 12)
Baseline confound: depleted call + background Chronos < -1 -> BASELINE_FITNESS_CONFOUNDED.
"""
import csv, json, math, sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, "/workspace/HPICDA/runner")
from checkpoint import write_checkpoint

ROOT = Path("/workspace/HPICDA")
FS = ROOT / "03_functional_screens"
OUT = ROOT / "05_conditional_dependency"
OUT.mkdir(exist_ok=True)
NOW = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")

CRISPRA = {"ORCS-1175", "ORCS-1706", "ORCS-1708", "ORCS-1863", "ORCS-2020"}
BOTTLENECK = {"ORCS-1150"}  # RRA primary, LFC sign unreliable (HG2 rule)

# ---------- registry ----------
with open(FS / "screen_estimand_registry.tsv") as f:
    reg = {r["screen_id"]: r for r in csv.DictReader(f, delimiter="\t")}

# effective track after amendments
def eff_track(sid, reg_track):
    if sid.startswith("ORCS-182"):
        return "B"
    if sid == "LIT-Pa-GSE154751":
        return "A"
    if sid == "LIT-Pa-SRP514827":
        return "UNAVAILABLE"
    if sid == "ORCS-2067":
        return "B"
    return reg_track

# ---------- baseline matrix ----------
baseline = {}  # background -> {gene: chronos}
with open(ROOT / "04_depmap_baseline" / "baseline_chronos_by_background.tsv") as f:
    rdr = csv.reader(f, delimiter="\t")
    header = next(rdr)
    bgs = header[1:]
    for bg in bgs:
        baseline[bg] = {}
    for row in rdr:
        g = row[0]
        for bg, v in zip(bgs, row[1:]):
            if v != "NA":
                baseline[bg][g] = float(v)

def get_baseline(cell_line, gene):
    # registry cell_line -> background key mapping
    m = {"A-549": "A-549", "Calu-3": "Calu-3", "THP-1": "THP-1", "HuH-7": "HuH-7",
         "Huh-7.5": "Huh-7.5", "Huh-7.5.1": "Huh-7.5.1", "HEK293T/ACE2": "HEK-293",
         "HEK293-ACE2-TMPRSS2": "HEK-293", "Caco-2-ACE2": "C2BBe1", "HAP-1": "HAP-1",
         "HCT 116": "HCT 116", "UM-UC-4": "UM-UC-4", "HT-29": "HT-29", "K-562": "K-562",
         "A-375": "A-375", "SH-SY5Y": "SH-SY5Y", "PLC/PRF/5": "PLC/PRF/5",
         "U-937": "U-937", "HeLa": "HeLa", "HEp-2": "HEp-2", "Vero-E6": "Vero-E6"}
    bg = m.get(cell_line)
    if bg and bg in baseline:
        return baseline[bg].get(gene)
    return None

# ---------- screen -> data file mapping ----------
# Track A: screen_id -> list of (subcontrast, mageck gene_summary path)
TRACKA = {
    "ORCS-1148": [("ORCS-1148", FS / "raw/GSE111166/ORCS-1148.gene_summary.txt")],
    "ORCS-1149": [("ORCS-1149", FS / "raw/GSE111166/ORCS-1149.gene_summary.txt")],
    "ORCS-1150": [("ORCS-1150", FS / "raw/GSE111166/ORCS-1150.gene_summary.txt")],
    "ORCS-1207": [("ORCS-1207-S1", FS / "raw/Li2020_Edinburgh/ORCS-1207-S1.gene_summary.txt"),
                  ("ORCS-1207-S2", FS / "raw/Li2020_Edinburgh/ORCS-1207-S2.gene_summary.txt")],
    "ORCS-1862": [("ORCS-1862", FS / "raw/GSE175666/ORCS-1862.gene_summary.txt")],
    "ORCS-1863": [("ORCS-1863", FS / "raw/GSE175666/ORCS-1863.gene_summary.txt")],
    "ORCS-1864": [("ORCS-1864", FS / "raw/GSE175666/ORCS-1864.gene_summary.txt")],
    "ORCS-2066": [("ORCS-2066", FS / "raw/GSE178699/ORCS-2066.gene_summary.txt")],
    "ORCS-2070": [("ORCS-2070", FS / "raw/GSE178699/ORCS-2070.gene_summary.txt")],
    "ORCS-2071": [("ORCS-2071", FS / "raw/GSE178699/ORCS-2071.gene_summary.txt")],
    "LIT-Mtb-GSE210689": [("active_replication", FS / "raw/GSE210689/mtb_pro_vs_ui.gene_summary.txt"),
                          ("stressed_nonreplicating", FS / "raw/GSE210689/mtb_lat_vs_ui.gene_summary.txt")],
    "LIT-RSV-PRJNA1281304-viability": [("libA", FS / "raw/PRJNA1281304/RSV_libA_viability.gene_summary.txt"),
                                       ("libB", FS / "raw/PRJNA1281304/RSV_libB_viability.gene_summary.txt")],
    "LIT-RSV-PRJNA1281304-early": [("libA_GFPNEG", FS / "raw/PRJNA1281304/RSV_libA_GFPNEG.gene_summary.txt"),
                                   ("libB_GFPNEG", FS / "raw/PRJNA1281304/RSV_libB_GFPNEG.gene_summary.txt"),
                                   ("libA_GFPLOW", FS / "raw/PRJNA1281304/RSV_libA_GFPLOW.gene_summary.txt"),
                                   ("libB_GFPLOW", FS / "raw/PRJNA1281304/RSV_libB_GFPLOW.gene_summary.txt")],
    "LIT-Pa-GSE154751": [("Pa_ExoU", FS / "raw/GSE154751/Pa_ExoU.gene_summary.txt")],
}

TRACKB_IDS = ["1175", "1220", "1385", "1386", "1431", "1433", "1581", "1583", "1584",
              "1704", "1706", "1708", "1821", "1822", "1823", "1824", "1825", "1826",
              "1827", "1955", "1956", "1997", "1998", "2020", "2067"]

def fnum(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return float("nan")

# ---------- Track A parser: MAGeCK gene_summary ----------
def parse_mageck(path, invert=False):
    """Return {gene: (direction_lof, score, fdr, lfc)} using winning side.
    score = -log10(permutation p-value) of winning side (bigger = more significant),
    consistent with Track B -log10 scale."""
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pf, nf = fnum(r["pos|fdr"]), fnum(r["neg|fdr"])
            pp, np_ = fnum(r["pos|p-value"]), fnum(r["neg|p-value"])
            if pf <= nf:
                d, pv, fd, lfc = "enriched", pp, pf, fnum(r["pos|lfc"])
            else:
                d, pv, fd, lfc = "depleted", np_, nf, fnum(r["neg|lfc"])
            s = -math.log10(pv) if pv > 0 else 300.0
            if invert:
                d = "depleted" if d == "enriched" else "enriched"
                lfc = -lfc if lfc == lfc else lfc
            out[r["id"]] = (d, s, fd, lfc)
    return out

# ---------- Track B parser: ORCS tables ----------
def parse_orcs(sid):
    """Return {gene: (direction_lof, score, fdr, lfc)}; score = -log10(p) or equivalent."""
    path = FS / "trackB_orcs" / f"ORCS-{sid}.gene_scores.tsv"
    out = {}
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        cols = rdr.fieldnames[1:]
        rows = list(rdr)
    n = len(rows)
    invert = f"ORCS-{sid}" in CRISPRA
    for r in rows:
        g = r["gene"]
        vals = {c: r[c] for c in cols}
        d = s = fd = lfc = None
        if sid in ("1220", "2020"):  # MaGeCK Score, p, FDR, Log2FC
            lfc = fnum(vals.get("Score.4Log2FC"))
            p = fnum(vals.get("Score.2p-Value"))
            fd = fnum(vals.get("Score.3FDR"))
            d = "enriched" if lfc > 0 else "depleted"
            s = -math.log10(p) if p > 0 else 300
        elif sid in ("1385", "1386"):  # Z-score, p, FDR
            z = fnum(vals.get("Score.1Z-score"))
            p = fnum(vals.get("Score.2p-Value"))
            fd = fnum(vals.get("Score.3FDR"))
            d = "enriched" if z > 0 else "depleted"
            s = abs(z)
            lfc = z
        elif sid in ("1431", "1433"):  # RRA rank, FDR, Log2, RRA score
            lfc = fnum(vals.get("Score.3Log2"))
            fd = fnum(vals.get("Score.2FDR"))
            rr = fnum(vals.get("Score.1RRA rank"))
            d = "enriched" if lfc > 0 else "depleted"
            s = -math.log10(rr / n) if rr > 0 else 0  # rank-based score
        elif sid in ("1581", "1583", "1584"):  # Log2FC, p
            lfc = fnum(vals.get("Score.1Log2FC"))
            p = fnum(vals.get("Score.2p-Value"))
            d = "enriched" if lfc > 0 else "depleted"
            s = -math.log10(p) if p > 0 else 300
            fd = float("nan")
        elif sid in ("1704", "1706"):  # MAGeCK pos score, FDR, Rank (one-sided enrichment)
            ps = fnum(vals.get("Score.1MAGeCK pos score"))
            fd = fnum(vals.get("Score.2FDR"))
            d = "enriched"
            s = -math.log10(ps) if ps > 0 else 300
        elif sid == "1708":  # MAGeCK neg score (one-sided depletion)
            ns = fnum(vals.get("Score.1MAGeCK neg score"))
            fd = fnum(vals.get("Score.2FDR"))
            d = "depleted"
            s = -math.log10(ns) if ns > 0 else 300
        elif sid in ("1821", "1822", "1823", "1824", "1825", "1826", "1827"):  # FDR, Log2FC
            lfc = fnum(vals.get("Score.2Log2FC"))
            fd = fnum(vals.get("Score.1FDR"))
            d = "enriched" if lfc > 0 else "depleted"
            s = -math.log10(fd) if fd and fd > 0 else (0 if fd == 0 else float("nan"))
        elif sid in ("1955", "1956", "1997", "1998"):  # Log2FC, FDR, FDR
            lfc = fnum(vals.get("Score.1Log2FC"))
            fd = fnum(vals.get("Score.2FDR"))
            d = "enriched" if lfc > 0 else "depleted"
            s = -math.log10(fd) if fd and fd > 0 else (0 if fd == 0 else float("nan"))
        elif sid == "2067":  # Differential score (qGI), FDR
            ds = fnum(vals.get("Score.1Differential score"))
            fd = fnum(vals.get("Score.2FDR"))
            d = "enriched" if ds > 0 else "depleted"
            s = abs(ds)
            lfc = ds
        elif sid == "1175":  # Score.1 Log10 (CRISPRa enrichment score)
            sc = fnum(vals.get("Score.1Log10"))
            d = "enriched"  # author one-sided enrichment
            s = sc
            fd = float("nan")
        else:
            continue
        if invert:
            d = "depleted" if d == "enriched" else "enriched"
            if lfc is not None and lfc == lfc:
                lfc = -lfc
        if d is not None:
            out[g] = (d, s, fd, lfc if lfc is not None else float("nan"))
    return out

def author_hit(sid, gene_row_vals, direction, fdr, score, lfc):
    """Author-defined hit thresholds per ORCS screen (from screen pages)."""
    if sid == "1175":
        return score > 4.0
    if sid in ("1220",):
        return fnum(gene_row_vals.get("Score.1MaGeCK Score")) < 0.005
    if sid in ("2020",):
        return fnum(gene_row_vals.get("Score.1MaGeCK Score")) < 0.00024
    if sid in ("1385", "1386"):
        return fdr == fdr and fdr < 0.25
    if sid in ("1431", "1433"):
        return fnum(gene_row_vals.get("Score.1RRA rank")) <= 1000
    if sid == "1581":
        return abs(lfc) >= 1.0 if lfc == lfc else False
    if sid == "1583":
        p = fnum(gene_row_vals.get("Score.2p-Value"))
        return p < 0.05 and lfc == lfc and lfc > 0
    if sid == "1584":
        return abs(lfc) >= 1.0 if lfc == lfc else False
    if sid in ("1704", "1706", "1708"):
        return fdr == fdr and fdr < 0.05
    if sid in ("1821", "1822", "1823", "1824", "1825", "1826", "1827"):
        return fdr == fdr and fdr < 0.25
    if sid in ("1955", "1956", "1997", "1998"):
        f2 = fnum(gene_row_vals.get("Score.2FDR"))
        f3 = fnum(gene_row_vals.get("Score.3FDR"))
        if lfc != lfc:
            return False
        if lfc > 0:
            return f2 == f2 and f2 < 0.45
        return f3 == f3 and f3 < 0.075
    if sid == "2067":
        return fdr == fdr and fdr < 0.1
    return False

# ---------- call assignment ----------
def assign_call(sid, regrow, direction, fdr, score, rank_pct, bl, sig=None, sig_hi=None):
    """Return (call, confidence, intervention)."""
    endpoint = regrow["endpoint_family"]
    if sid == "ORCS-1207":
        endpoint = "PATHOGEN_BURDEN"  # AM2-10
    if sig is None:
        sig_hi = (fdr == fdr and fdr < 0.05)
        sig = (fdr == fdr and fdr < 0.25)
    if endpoint == "OTHER_DEFINED":
        return ("NO_FUNCTIONAL_SIGNAL", "na", "none")  # SRP514827 unavailable
    if endpoint == "PATHOGEN_STATE":
        # Mtb state-reporter: intervention direction deferred to Stage 12 (HG1).
        # Exploratory tier: nominal RRA score >= -log10(0.01) ~ 2.0 (underpowered study,
        # feeds Stage 7 module analysis only, never gene-level claims).
        exploratory = (score == score and score >= 2.0)
        if sig:
            return ("SIGNAL_STATE_ALTERING", conf, "DIRECTION_UNRESOLVED")
        if exploratory:
            return ("SIGNAL_STATE_ALTERING", "exploratory", "DIRECTION_UNRESOLVED")
        return ("NO_FUNCTIONAL_SIGNAL", "na", "none")
    if not sig:
        return ("NO_FUNCTIONAL_SIGNAL", "na", "none")
    conf = "high" if sig_hi else "standard"
    if endpoint == "TOXIN_SENSITIVITY":
        return ("TOXIN_SPECIFIC", conf, "context_only")
    # HOST_FITNESS / PATHOGEN_BURDEN
    if direction == "enriched":
        # KO protects host / reduces pathogen
        note = ""
        if bl is not None and bl < -1:
            note = "_despite_essential"  # enriched despite baseline essentiality = stronger
            conf = conf + "+"
        return ("CDEP_PATHOGEN_EXPLOITED" + ("" if not note else ""), conf, "INHIBIT")
    else:
        # KO sensitizes host / increases pathogen; check baseline confound
        if bl is not None and bl < -1:
            return ("BASELINE_FITNESS_CONFOUNDED", "na", "none")
        return ("CDEP_HOST_PROTECTIVE", conf, "PRESERVE_ACTIVATE")

# ---------- main ----------
rows_out = []
screen_stats = defaultdict(lambda: defaultdict(int))

# Track A
for sid, subs in TRACKA.items():
    regrow = reg.get(sid)
    if not regrow:
        continue
    invert = sid in CRISPRA
    for sub, path in subs:
        if not path.exists():
            print(f"MISSING {path}")
            continue
        data = parse_mageck(path, invert=invert)
        # rank pct by score
        ranked = sorted(data.items(), key=lambda kv: -(kv[1][1] if kv[1][1] == kv[1][1] else -1))
        rp = {g: (i + 1) / len(ranked) for i, (g, _) in enumerate(ranked)}
        unit = f"{sid}|{sub}"
        for g, (d, s, fd, lfc) in data.items():
            bl = get_baseline(regrow["cell_line"], g)
            call, conf, interv = assign_call(sid, regrow, d, fd, s, rp[g], bl)
            screen_stats[unit][call] += 1
            rows_out.append([unit, sid, sub, regrow["pathogen"], regrow["cell_line"],
                             regrow["endpoint_family"], "A", g, d,
                             f"{s:.4g}" if s == s else "NA",
                             f"{fd:.4g}" if fd == fd else "NA",
                             f"{lfc:.4f}" if lfc == lfc else "NA",
                             f"{rp[g]:.5f}",
                             f"{bl:.3f}" if bl is not None else "NA",
                             call, conf, interv])

# Track B
for sidnum in TRACKB_IDS:
    sid = f"ORCS-{sidnum}"
    regrow = reg.get(sid)
    if not regrow:
        continue
    data = parse_orcs(sidnum)
    # reload raw row values for author-threshold evaluation
    raw_vals = {}
    with open(FS / "trackB_orcs" / f"ORCS-{sidnum}.gene_scores.tsv") as fh:
        for rr in csv.DictReader(fh, delimiter="\t"):
            raw_vals[rr["gene"]] = rr
    ranked = sorted(data.items(), key=lambda kv: -(kv[1][1] if kv[1][1] == kv[1][1] else -1))
    rp = {g: (i + 1) / len(ranked) for i, (g, _) in enumerate(ranked)}
    for g, (d, s, fd, lfc) in data.items():
        bl = get_baseline(regrow["cell_line"], g)
        ah = author_hit(sidnum, raw_vals.get(g, {}), d, fd, s, lfc)
        sig = ah
        sig_hi = bool(ah and fd == fd and fd < 0.05)
        call, conf, interv = assign_call(sid, regrow, d, fd, s, rp[g], bl, sig=sig, sig_hi=sig_hi)
        screen_stats[sid][call] += 1
        rows_out.append([sid, sid, "", regrow["pathogen"], regrow["cell_line"],
                         regrow["endpoint_family"], "B", g, d,
                         f"{s:.4g}" if s == s else "NA",
                         f"{fd:.4g}" if fd == fd else "NA",
                         f"{lfc:.4f}" if lfc == lfc else "NA",
                         f"{rp[g]:.5f}",
                         f"{bl:.3f}" if bl is not None else "NA",
                         call, conf, interv])

header = ["unit_id", "screen_id", "subcontrast", "pathogen", "cell_line",
          "endpoint_family", "track", "gene", "direction_lof", "score", "fdr",
          "lfc", "rank_pct", "baseline_chronos", "call", "confidence",
          "intervention_implication"]
out_path = OUT / "per_screen_calls.tsv"
with open(out_path, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(header)
    w.writerows(rows_out)
print(f"wrote {out_path}: {len(rows_out)} rows, {len(screen_stats)} screen-units")

# screen-level summary
sum_path = OUT / "per_screen_call_counts.tsv"
with open(sum_path, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["unit_id", "CDEP_PATHOGEN_EXPLOITED", "CDEP_HOST_PROTECTIVE",
                "BASELINE_FITNESS_CONFOUNDED", "TOXIN_SPECIFIC",
                "SIGNAL_STATE_ALTERING", "NO_FUNCTIONAL_SIGNAL"])
    for unit in sorted(screen_stats):
        st = screen_stats[unit]
        w.writerow([unit, st.get("CDEP_PATHOGEN_EXPLOITED", 0),
                    st.get("CDEP_HOST_PROTECTIVE", 0),
                    st.get("BASELINE_FITNESS_CONFOUNDED", 0),
                    st.get("TOXIN_SPECIFIC", 0),
                    st.get("SIGNAL_STATE_ALTERING", 0),
                    st.get("NO_FUNCTIONAL_SIGNAL", 0)])
print(f"wrote {sum_path}")

with open(OUT / "stage5_screen_stats.json", "w") as f:
    json.dump({k: dict(v) for k, v in screen_stats.items()}, f, indent=2)
print("done")
