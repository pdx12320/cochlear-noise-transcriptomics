#!/usr/bin/env python3
import csv
import os
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "pr": "http://schemas.openxmlformats.org/package/2006/relationships"}


def col_index(ref):
    letters = "".join(ch for ch in ref if ch.isalpha())
    out = 0
    for ch in letters:
        out = out * 26 + ord(ch.upper()) - 64
    return out - 1


def rows_from_sheet(zf, target, shared):
    root = ET.fromstring(zf.read(target))
    for row in root.findall("m:sheetData/m:row", NS):
        vals = []
        for cell in row.findall("m:c", NS):
            idx = col_index(cell.attrib.get("r", "A1"))
            while len(vals) <= idx:
                vals.append("")
            ctype = cell.attrib.get("t")
            value = cell.find("m:v", NS)
            if ctype == "inlineStr":
                inline = cell.find("m:is", NS)
                text = "" if inline is None else "".join(t.text or "" for t in inline.findall(".//m:t", NS))
            elif value is None:
                text = ""
            elif ctype == "s":
                text = shared[int(value.text)]
            else:
                text = value.text or ""
            vals[idx] = text
        yield vals


src, dest = sys.argv[1:3]
with zipfile.ZipFile(src) as zf:
    shared = []
    if "xl/sharedStrings.xml" in zf.namelist():
        root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
        for si in root.findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.findall(".//m:t", NS)))
    wb = ET.fromstring(zf.read("xl/workbook.xml"))
    relroot = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rels = {x.attrib["Id"]: x.attrib["Target"] for x in relroot.findall("pr:Relationship", NS)}
    sheets = []
    for sheet in wb.findall("m:sheets/m:sheet", NS):
        target = rels[sheet.attrib[f"{{{NS['r']}}}id"]].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        sheets.append((sheet.attrib["name"], target))

    records = []
    for sheet_name, target in sheets[1:]:
        rows = list(rows_from_sheet(zf, target, shared))
        header_index = next(i for i, row in enumerate(rows) if "p_val" in row and "gene" in row)
        header = [x.strip() for x in rows[header_index]]
        for row in rows[header_index + 1:]:
            row = row + [""] * (len(header) - len(row))
            item = dict(zip(header, row))
            if not item.get("gene"):
                continue
            item["author_sheet"] = sheet_name
            records.append(item)

fields = ["author_sheet", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"]
os.makedirs(os.path.dirname(dest), exist_ok=True)
with open(dest, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(records)
print(f"WROTE\t{dest}\t{len(records)}")
