#!/usr/bin/env python3
"""Build the reader-facing HPICDA supplementary workbook.

The scientific TSVs remain the source of truth. This script packages the
manageable tables into one indexed XLSX. The 852,216-row S5 call table remains
as its compressed TSV because embedding it would produce an unnecessarily
large and fragile Excel file.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

from openpyxl import Workbook


ROOT = Path(
    os.environ.get("HPICDA_ROOT", Path(__file__).resolve().parents[1])
).resolve()
PACKAGE = ROOT / "25_manuscript_package/revision_v2_2026-10-21"
TABLE_SOURCE = ROOT / "25_manuscript_package/tables"
INDEX = TABLE_SOURCE / "TABLES_INDEX.json"
OUTPUT = PACKAGE / "HPICDA_Supplementary_Tables_S1-S21b.xlsx"


def clean_cell(value: str) -> str:
    value = value.replace("\x00", "")
    if len(value) > 32767:
        return value[:32700] + " [TRUNCATED; SEE SOURCE TSV]"
    return value


def main() -> None:
    entries = json.loads(INDEX.read_text(encoding="utf-8"))
    workbook = Workbook(write_only=True)

    contents = workbook.create_sheet("Contents")
    contents.append(
        [
            "Table",
            "Workbook sheet or external file",
            "Description",
            "Rows",
            "Columns",
        ]
    )

    included: list[dict[str, object]] = []
    external: list[dict[str, object]] = []

    for entry in entries:
        source = TABLE_SOURCE / str(entry["file"])
        if not source.exists():
            raise FileNotFoundError(source)
        if source.suffix == ".gz":
            external.append(entry)
            contents.append(
                [
                    entry["id"],
                    "Zenodo data archive: 10.5281/zenodo.23550165",
                    entry["caption"],
                    entry["rows"],
                    entry["cols"],
                ]
            )
            continue

        sheet = workbook.create_sheet(str(entry["id"])[:31])
        with source.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            rows = 0
            for raw_row in reader:
                row = [clean_cell(value) for value in raw_row]
                sheet.append(row)
                rows += 1
        contents.append(
            [
                entry["id"],
                entry["id"],
                entry["caption"],
                entry["rows"],
                entry["cols"],
            ]
        )
        included.append(entry)

    readme = workbook.create_sheet("README", 0)
    for row in [
        ["HPICDA supplementary tables — revision 2026-10-21"],
        [
            "This workbook contains the manageable supplementary tables. "
            "The full S5 per-screen call table contains 852,216 rows and is "
            "supplied through the Zenodo data archive "
            "(DOI 10.5281/zenodo.23550165), because its compressed file exceeds "
            "the 15 MB initial-submission limit. S5b gives the per-screen summary."
        ],
        [
            "S21a and S21b were added in this revision for the opportunity-adjusted "
            "and equal-opportunity benchmark analyses."
        ],
        [
            "The underlying TSV files remain the source of truth and are retained "
            "in the submission package."
        ],
    ]:
        readme.append(row)

    workbook.save(OUTPUT)
    print(f"Wrote {OUTPUT}")
    print(f"Included worksheet tables: {len(included)}")
    print(f"External compressed tables: {len(external)}")


if __name__ == "__main__":
    main()
