#!/usr/bin/env python3
"""HPICDA — memory-safe streaming xlsx sheet extractor (read-only).

Usage: python3 09_xlsx_stream.py FILE.xlsx SHEET_INDEX OUT.tsv [MAXROWS]
Sheet index is 1-based over xl/worksheets/sheetN.xml. Prints header + dims.
"""
import sys, zipfile, re
import xml.etree.ElementTree as ET

def main():
    path, sheet_idx, out_path = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    maxrows = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for ev, el in ET.iterparse(z.open("xl/sharedStrings.xml"), events=("end",)):
            if el.tag == ns + "si":
                shared.append("".join(t.text or "" for t in el.iter(ns + "t")))
                el.clear()
    nrows = 0
    ncols = 0
    with open(out_path, "w") as out:
        for ev, row in ET.iterparse(z.open(f"xl/worksheets/sheet{sheet_idx}.xml"), events=("end",)):
            if row.tag != ns + "row":
                continue
            # map cells by column reference to handle sparse rows
            cells = {}
            maxc = -1
            for c in row.iter(ns + "c"):
                ref = c.get("r", "")
                m = re.match(r"([A-Z]+)(\d+)", ref)
                if m:
                    col = 0
                    for ch in m.group(1):
                        col = col * 26 + (ord(ch) - 64)
                    col -= 1
                else:
                    col = maxc + 1
                v = c.find(ns + "v")
                txt = v.text if v is not None else ""
                if c.get("t") == "s" and txt != "":
                    txt = shared[int(txt)]
                elif c.get("t") == "inlineStr":
                    is_el = c.find(ns + "is")
                    txt = "".join(t.text or "" for t in is_el.iter(ns + "t")) if is_el is not None else ""
                cells[col] = txt
                maxc = max(maxc, col)
            if cells:
                width = max(cells) + 1
                ncols = max(ncols, width)
                out.write("\t".join(cells.get(i, "") for i in range(width)) + "\n")
                nrows += 1
                if maxrows and nrows >= maxrows:
                    break
            row.clear()
    print(f"rows={nrows} cols={ncols} -> {out_path}")

if __name__ == "__main__":
    main()
