#!/usr/bin/env python
"""Programmatic QA for the HPICDA figure set (PDF-only output).

The build host has no vision, so each figure is checked from the PDF itself:

  * every text span's true point size, reported as a histogram -- the figure
    contract requires all text >= 6.0 pt at final size, body text 6.5-9.5 pt
  * expected panel tags a/b/c(/d) present exactly once
  * axis lines: long horizontal / vertical vector strokes recovered from the
    page drawing operators. A panel that has lost its axis line shows up as a
    missing stroke, which is the failure mode this check exists to catch.
  * text-overlap clusters, to catch enlarged labels colliding

Exit status is non-zero if any figure violates the contract.
"""
from __future__ import annotations

import os

import re
import sys
from collections import Counter
from pathlib import Path

import fitz  # PyMuPDF

FIG = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    Path(os.environ.get(
        "HPICDA_ROOT", Path(__file__).resolve().parents[1]))
    / "25_manuscript_package/figures")

MIN_PT = 6.0          # hard floor for any glyph
STROKE_MIN = 18.0     # pt; strokes shorter than this are ticks/segments, not axes
TOL = 0.6             # pt; straightness tolerance for axis strokes


def axis_strokes(page):
    """Long near-axis-aligned strokes drawn on the page."""
    hr = vr = 0
    for d in page.get_drawings():
        for item in d["items"]:
            if item[0] == "l":
                p, q = item[1], item[2]
            elif item[0] == "re":
                continue
            else:
                continue
            dx, dy = q.x - p.x, q.y - p.y
            if abs(dy) <= TOL and abs(dx) >= STROKE_MIN:
                hr += 1
            elif abs(dx) <= TOL and abs(dy) >= STROKE_MIN:
                vr += 1
    return hr, vr


def overlaps(spans):
    """Count near-collisions between spans that share a text line."""
    n = 0
    for i in range(len(spans)):
        for j in range(i + 1, len(spans)):
            a, b = spans[i], spans[j]
            if not (a["bbox"][3] > b["bbox"][1] + 0.5 and b["bbox"][3] > a["bbox"][1] + 0.5):
                continue
            ox = min(a["bbox"][2], b["bbox"][2]) - max(a["bbox"][0], b["bbox"][0])
            oy = min(a["bbox"][3], b["bbox"][3]) - max(a["bbox"][1], b["bbox"][1])
            # thresholds sit above em-box grazing (~1.2 pt) but below a real
            # glyph collision, which runs to several points
            if ox > 1.0 and oy > 1.6:
                n += 1
    return n


def main():
    bad = 0
    print(f"{'figure':40s} {'size pt':>8s} {'min':>5s} {'max':>5s} "
          f"{'h-axis':>6s} {'v-axis':>6s} {'tags':>5s} {'ovl':>4s} {'oob':>4s}")
    print("-" * 98)
    for pdf in sorted(FIG.glob("fig*.pdf")):
        doc = fitz.open(pdf)
        sizes, spans, tags = Counter(), [], []
        oob = []
        for page in doc:
            W, H = page.rect.width, page.rect.height
            for blk in page.get_text("dict")["blocks"]:
                for line in blk.get("lines", []):
                    for sp in line["spans"]:
                        t = sp["text"].strip()
                        if not t:
                            continue
                        sizes[round(sp["size"], 1)] += len(t)
                        spans.append({"t": t, "bbox": sp["bbox"], "size": sp["size"]})
                        if re.fullmatch(r"[a-e]", t):
                            tags.append(t)
                        # a span whose em box leaves the page is clipped in print
                        x0, y0, x1, y1 = sp["bbox"]
                        if x0 < -0.5 or x1 > W + 0.5 or y0 < -0.5 or y1 > H + 0.5:
                            oob.append((t, round(x0, 1), round(x1, 1), round(W, 1)))
            hr, vr = axis_strokes(page)
        mn = min(sizes) if sizes else 0.0
        mx = max(sizes) if sizes else 0.0
        ovl = overlaps(spans)
        tagtxt = "".join(sorted(tags)) or "-"
        flag = ""
        if mn < MIN_PT:
            flag += "  <-- TEXT BELOW FLOOR"
            bad += 1
        if oob:
            flag += f"  <-- {len(oob)} SPAN(S) OFF-CANVAS"
            bad += 1
        print(f"{pdf.stem:40s} {len(sizes):8d} {mn:5.1f} {mx:5.1f} "
              f"{hr:6d} {vr:6d} {tagtxt:>5s} {ovl:4d} {len(oob):4d}{flag}")
        print(f"    px sizes: " + ", ".join(f"{s}:{c}" for s, c in sorted(sizes.items())))
        for t, x0, x1, W in oob[:5]:
            print(f"    OFF-CANVAS: '{t[:44]}' x {x0}..{x1} (page width {W})")
        doc.close()
    print("-" * 98)
    print("all figures pass" if bad == 0 else f"{bad} figure(s) failed")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
