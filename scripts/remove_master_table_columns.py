#!/usr/bin/env python3

import csv
import json
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "backend" / "data"
DROP_COLUMNS = [
    "苯乙酮_活性口袋命中Zn",
    "Zn第一配位壳最近空隙_A",
    "ZN_弯曲通道半径_A",
    "ZN_弯曲通道长度_A",
]


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload):
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def clean_header(header):
    out = []
    seen = set()
    for col in header or []:
        key = str(col or "").strip()
        if not key or key in seen or key in DROP_COLUMNS:
            continue
        seen.add(key)
        out.append(key)
    return out


def clean_row(row):
    if not isinstance(row, dict):
        return {}
    return {
        key: value
        for key, value in row.items()
        if str(key or "").strip() and str(key or "").strip() not in DROP_COLUMNS
    }


def cleanup_master_full(path: Path):
    obj = read_json(path)
    old_rows = list(obj.get("rows") or [])
    old_header = list(obj.get("header") or [])
    header = clean_header(old_header)
    rows = []
    changed = header != old_header
    for src in old_rows:
        row = clean_row(src)
        if len(row) != len(src):
            changed = True
        rows.append({col: row.get(col, "") for col in header})
    if not changed:
        return False
    obj["header"] = header
    obj["rows"] = rows
    write_json(path, obj)
    return True


def cleanup_csv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        old_fieldnames = list(reader.fieldnames or [])
    fieldnames = clean_header(old_fieldnames)
    changed = fieldnames != old_fieldnames
    clean_rows = []
    for src in rows:
        row = clean_row(src)
        if len(row) != len(src):
            changed = True
        clean_rows.append({col: row.get(col, "") for col in fieldnames})
    if not changed:
        return False
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in clean_rows:
            writer.writerow(row)
    return True


def cleanup_data_json(path: Path):
    obj = read_json(path)
    items = list(obj.get("items") or [])
    changed = False
    clean_items = []
    for src in items:
        row = clean_row(src)
        if len(row) != len(src):
            changed = True
        clean_items.append(row)
    if not changed:
        return False
    obj["items"] = clean_items
    write_json(path, obj)
    return True


def iter_paths(pattern: str, singleton: str):
    paths = {DATA_DIR / singleton}
    paths.update(DATA_DIR.glob(pattern))
    return sorted(path for path in paths if path.exists())


def main():
    changed_full = sum(cleanup_master_full(path) for path in iter_paths("*_master_full.json", "master_full.json"))
    changed_csv = sum(cleanup_csv(path) for path in iter_paths("*_master_table.csv", "master_table.csv"))
    changed_data = sum(cleanup_data_json(path) for path in iter_paths("*_data.json", "data.json"))
    print(
        "cleaned",
        f"master_full={changed_full}",
        f"master_csv={changed_csv}",
        f"data_json={changed_data}",
    )


if __name__ == "__main__":
    main()
