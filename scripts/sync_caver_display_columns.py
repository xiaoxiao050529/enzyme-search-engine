#!/usr/bin/env python3

import csv
import json
import re
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "backend" / "data"
DISPLAY_COLUMNS = [
    "苯乙酮_CAVER通道判定",
    "苯乙酮_CAVER口袋深度",
    "苯乙酮_CAVER口袋深度判定列",
    "苯乙酮_CAVER口袋深度判定值",
]


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload):
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def row_id(row):
    return str((row or {}).get("Representative") or (row or {}).get("_id") or (row or {}).get("id") or "").strip().upper()


def parse_depth_basis(text):
    raw = str(text or "").strip()
    if not raw:
        return "", ""
    match = re.match(r"^([A-Za-z0-9_]+)\s*=\s*(.+)$", raw)
    if not match:
        return "", raw
    return match.group(1).strip(), match.group(2).strip()


def parse_caver_note_fields(note):
    text = str(note or "").strip()
    best_match = re.search(r"best_bottleneck=([^;；。]+)", text)
    channel_match = re.search(r"通道判定=([^;；。]+)", text)
    depth_match = re.search(r"口袋深度=([^（(；;。]+)[（(]([^）)]+)[）)]", text)
    return {
        "best_bottleneck": best_match.group(1).replace("Å", "").strip() if best_match else "",
        "channel_judge": channel_match.group(1).strip() if channel_match else "",
        "depth_label": depth_match.group(1).strip() if depth_match else "",
        "depth_basis": depth_match.group(2).strip() if depth_match else "",
    }


def extract_original_note(note):
    text = str(note or "").strip()
    if "原备注：" in text:
        return text.split("原备注：", 1)[1].strip()
    return ""


def build_note(row, original_note):
    seed_mode = str(row.get("CAVER_SeedMode", "") or "").strip()
    channel_count = str(row.get("CAVER_ChannelCount", "") or "").strip()
    if not (seed_mode or channel_count or original_note):
        return str(row.get("苯乙酮_判断说明", "") or "").strip()
    note = f"CAVER pocket-seed 重算：{seed_mode}; channels={channel_count}。详细判定见专用列。".strip()
    if original_note:
        note += f" 原备注：{original_note}"
    return note


def apply_display_fields(row):
    row = dict(row or {})
    parsed = parse_caver_note_fields(row.get("苯乙酮_判断说明", ""))
    best_bottleneck = str(row.get("苯乙酮_通道瓶颈半径_A", "") or "").replace("Å", "").strip()
    if not best_bottleneck:
        best_bottleneck = str(row.get("CAVER_BestBottleneckRadius_A", "") or "").replace("Å", "").strip() or parsed["best_bottleneck"]
        if best_bottleneck:
            row["苯乙酮_通道瓶颈半径_A"] = best_bottleneck
    channel_judge = str(row.get("苯乙酮_CAVER通道判定", "") or "").strip() or str(row.get("CAVER_AcetophenonePass", "") or "").strip() or parsed["channel_judge"]
    depth_label = str(row.get("苯乙酮_CAVER口袋深度", "") or "").strip() or str(row.get("CAVER_PocketDepthLabel", "") or "").strip() or parsed["depth_label"]
    depth_basis = str(row.get("CAVER_PocketDepthBasis", "") or "").strip() or parsed["depth_basis"]
    basis_col, basis_value = parse_depth_basis(depth_basis)
    if channel_judge:
        row["苯乙酮_CAVER通道判定"] = channel_judge
    if depth_label:
        row["苯乙酮_CAVER口袋深度"] = depth_label
    if basis_col:
        row["苯乙酮_CAVER口袋深度判定列"] = basis_col
    if basis_value:
        row["苯乙酮_CAVER口袋深度判定值"] = basis_value
    original_note = extract_original_note(row.get("苯乙酮_判断说明", ""))
    if any(str(row.get(col, "") or "").strip() for col in ["CAVER_SeedMode", "CAVER_ChannelCount", "CAVER_AcetophenonePass", "CAVER_PocketDepthLabel", "CAVER_PocketDepthBasis"]):
        row["苯乙酮_判断说明"] = build_note(row, original_note)
    return row


def should_include_display_columns(rows):
    for row in rows:
        if any(str((row or {}).get(col, "") or "").strip() for col in DISPLAY_COLUMNS):
            return True
        if any(str((row or {}).get(col, "") or "").strip() for col in ["CAVER_BestBottleneckRadius_A", "CAVER_AcetophenonePass", "CAVER_PocketDepthLabel", "CAVER_PocketDepthBasis"]):
            return True
    return False


def sync_bundle(full_path: Path):
    prefix = full_path.name[:-len("_master_full.json")] if full_path.name.endswith("_master_full.json") else ""
    data_name = f"{prefix}_data.json" if prefix else "data.json"
    csv_name = f"{prefix}_master_table.csv" if prefix else "master_table.csv"
    data_path = DATA_DIR / data_name
    csv_path = DATA_DIR / csv_name

    full_obj = read_json(full_path)
    header = list(full_obj.get("header") or [])
    old_rows = list(full_obj.get("rows") or [])
    rows = [apply_display_fields(row) for row in old_rows]
    changed = rows != old_rows

    if should_include_display_columns(rows):
      for col in DISPLAY_COLUMNS:
          if col not in header:
              header.append(col)
              changed = True

    if changed:
        full_obj["header"] = header
        full_obj["rows"] = [{col: row.get(col, "") for col in header} for row in rows]
        write_json(full_path, full_obj)

        if csv_path.exists():
            with csv_path.open("w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(header)
                for row in rows:
                    writer.writerow([row.get(col, "") for col in header])

        if data_path.exists():
            data_obj = read_json(data_path)
            items = list(data_obj.get("items") or [])
            item_map = {row_id(item): dict(item) for item in items if row_id(item)}
            out_items = []
            for row in rows:
                rid = row_id(row)
                if not rid:
                    continue
                item = item_map.get(rid, {"id": rid})
                for col in ["苯乙酮_通道瓶颈半径_A", *DISPLAY_COLUMNS, "苯乙酮_判断说明"]:
                    if str(row.get(col, "") or "").strip():
                        item[col] = row.get(col, "")
                out_items.append(item)
            data_obj["items"] = out_items
            write_json(data_path, data_obj)
    return changed


def main():
    changed = 0
    for full_path in sorted(DATA_DIR.glob("*_master_full.json")):
        changed += 1 if sync_bundle(full_path) else 0
    for singleton in [DATA_DIR / "master_full.json"]:
        if singleton.exists():
            changed += 1 if sync_bundle(singleton) else 0
    print("synced", changed)


if __name__ == "__main__":
    main()
