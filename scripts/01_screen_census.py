#!/usr/bin/env python3
"""HPICDA Stage 0 — Screen census via BioGRID ORCS 2.0.18 search endpoint.

Collects SCREEN-LEVEL METADATA ONLY (screen id, publication, organism, cell
line, phenotype, condition, library, analysis). Per the blueprint leakage
rules, this script never touches gene-level scores (.screen.tab.txt content).

Usage:
    python scripts/01_screen_census.py --out-dir 02_metadata
"""
from __future__ import annotations

import argparse
import html
import json
import re
import time
from pathlib import Path

import urllib.request
import urllib.parse
import http.cookiejar

BASE = "https://orcs.thebiogrid.org"

# Census keyword set: 8 core pathogens + catch-alls. Each keyword is queried
# against organism=9606 (human) and organism=all (to catch mis-annotation).
CORE_PATHOGEN_KEYWORDS = {
    "influenza": "Influenza A virus",
    "coronavirus": "SARS-CoV-2",
    "SARS-CoV-2": "SARS-CoV-2",
    "respiratory syncytial": "RSV",
    "tuberculosis": "Mycobacterium tuberculosis",
    "mycobacterium": "Mycobacterium tuberculosis",
    "pseudomonas": "Pseudomonas aeruginosa",
    "klebsiella": "Klebsiella pneumoniae",
    "staphylococcus": "Staphylococcus aureus",
    "streptococcus": "Streptococcus pneumoniae",
    "pneumococcus": "Streptococcus pneumoniae",
}
CATCHALL_KEYWORDS = ["bacteria", "toxin", "pneumonia", "macrophage infection"]


def strip_tags(x: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", " ", x)).strip()


def collapse_ws(x: str) -> str:
    return re.sub(r"\s+", " ", x).strip()


class OrcsClient:
    def __init__(self, delay: float = 1.5):
        self.delay = delay
        self.cj = http.cookiejar.CookieJar()
        self.opener = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(self.cj)
        )
        self.opener.addheaders = [("User-Agent", "HPICDA-stage0-census/1.0")]

    def _get(self, url: str) -> bytes:
        time.sleep(self.delay)
        with self.opener.open(url, timeout=90) as r:
            return r.read()

    def _post(self, url: str, data: dict) -> bytes:
        time.sleep(self.delay)
        body = urllib.parse.urlencode(data).encode()
        req = urllib.request.Request(url, data=body, method="POST")
        with self.opener.open(req, timeout=90) as r:
            return r.read()

    def search(self, keyword: str, organism: str = "9606") -> list[dict]:
        """Session-based search: GET primes the session, POST fetches rows."""
        q = urllib.parse.urlencode(
            {"searchType": "key", "search": keyword, "organism": organism}
        )
        self._get(f"{BASE}/Search/?{q}")
        exp = {
            "tool": "searchRows",
            "draw": 1,
            "start": 0,
            "length": 1000,
            "totalRecords": "0",
            "checkedBoxes": [],
        }
        raw = self._post(
            f"{BASE}/scripts/datatableTools.php", {"expData": json.dumps(exp)}
        )
        payload = json.loads(raw)
        rows = []
        for r in payload.get("data", []):
            rows.append(self._parse_row(r, keyword, organism))
        return rows

    @staticmethod
    def _parse_row(r: list, keyword: str, organism: str) -> dict:
        cells = [collapse_ws(strip_tags(str(c))) for c in r]
        raw = [str(c) for c in r]
        screen_id = ""
        m = re.search(r"value='(\d+)'", raw[0]) if raw else None
        if m:
            screen_id = m.group(1)
        screen_name = ""
        pmid = ""
        for c in raw:
            m2 = re.search(r"/Screen/(\d+)' title='([^']*)'", c)
            if m2:
                screen_id = screen_id or m2.group(1)
                screen_name = html.unescape(m2.group(2))
            m3 = re.search(r"PMID(\d+)", c)
            if m3:
                pmid = m3.group(1)
        # Column order observed: checkbox, Dataset, Screen, Organism,
        # Cell Type, Cell Line, Phenotype, Condition, Library, Enzyme,
        # Analysis, Score details, ...
        return {
            "orcs_screen_id": screen_id,
            "query_keyword": keyword,
            "query_organism": organism,
            "dataset_ref": cells[1] if len(cells) > 1 else "",
            "screen_name": screen_name or (cells[2] if len(cells) > 2 else ""),
            "pmid": pmid,
            "organism": cells[3] if len(cells) > 3 else "",
            "cell_type": cells[4] if len(cells) > 4 else "",
            "cell_line": cells[5] if len(cells) > 5 else "",
            "phenotype": cells[6] if len(cells) > 6 else "",
            "condition": cells[7] if len(cells) > 7 else "",
            "library_summary": cells[8] if len(cells) > 8 else "",
            "enzyme": cells[9] if len(cells) > 9 else "",
            "analysis_summary": cells[10] if len(cells) > 10 else "",
            "score_summary": cells[11] if len(cells) > 11 else "",
        }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="02_metadata")
    ap.add_argument("--delay", type=float, default=1.5)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    client = OrcsClient(delay=args.delay)

    all_rows: list[dict] = []
    keywords = list(CORE_PATHOGEN_KEYWORDS) + CATCHALL_KEYWORDS
    for kw in keywords:
        for org in ("9606", "all"):
            try:
                rows = client.search(kw, org)
            except Exception as e:  # noqa: BLE001 - log and continue census
                rows = []
                print(f"WARN query failed kw={kw!r} org={org}: {e}")
            print(f"kw={kw!r:24s} org={org:5s} -> {len(rows)} rows")
            all_rows.extend(rows)

    # Deduplicate by screen id, keeping the union of query keywords.
    by_id: dict[str, dict] = {}
    for row in all_rows:
        sid = row["orcs_screen_id"]
        if not sid:
            continue
        if sid in by_id:
            prev = by_id[sid]
            kws = set(prev["query_keyword"].split(";")) | {row["query_keyword"]}
            prev["query_keyword"] = ";".join(sorted(kws))
        else:
            by_id[sid] = row

    rows = sorted(by_id.values(), key=lambda r: int(r["orcs_screen_id"]))
    raw_path = out / "orcs_census_raw.json"
    raw_path.write_text(json.dumps(rows, indent=2, ensure_ascii=False))

    cols = [
        "orcs_screen_id", "pmid", "dataset_ref", "screen_name", "organism",
        "cell_type", "cell_line", "phenotype", "condition", "library_summary",
        "enzyme", "analysis_summary", "score_summary", "query_keyword",
        "query_organism",
    ]
    tsv_path = out / "orcs_census_screens.tsv"
    with tsv_path.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(c, "") for c in cols) + "\n")

    print(f"\nUnique screens: {len(rows)}")
    print(f"Wrote {raw_path} and {tsv_path}")


if __name__ == "__main__":
    main()
