#!/usr/bin/env python3
import csv
import json
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "backend" / "data"
REPORT_CSV = PROJECT_DIR / "backend" / "runtime" / "reports" / "workflow_table_2_caver_pocket_seed_v1_acetophenone_summary.csv"
TARGET_TABLES = ["workflow_table_3", "workflow_table_5"]


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload):
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def build_summary():
    with REPORT_CSV.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        rows = list(csv.DictReader(f))
    groups = {"yes": [], "borderline": [], "no": []}
    for row in rows:
        key = str(row.get("Acetophenone_size_pass_by_caver", "") or "").strip()
        if key not in groups:
            continue
        groups[key].append(
            {
                "pdb_id": row.get("PDB_ID", ""),
                "best_bottleneck_radius_a": row.get("Best_bottleneck_radius_A", ""),
            }
        )
    summary = {
        "method": {
            "seed_mode": "best_pocket_center",
            "probe_radius": 0.6,
            "shell_radius": 5,
            "shell_depth": 4,
            "pass_radius_a": 2.2,
            "borderline_radius_a": 1.8,
        },
        "counts": {k: len(v) for k, v in groups.items()},
        "groups": groups,
        "text": (
            "CAVER pocket-seed 重算结论："
            f"yes={len(groups['yes'])}, "
            f"borderline={len(groups['borderline'])}, "
            f"no={len(groups['no'])}."
        ),
    }
    return summary


def attach_to_table(table_id: str, summary: dict):
    full_path = DATA_DIR / f"{table_id}_master_full.json"
    data_path = DATA_DIR / f"{table_id}_data.json"
    full_obj = read_json(full_path)
    data_obj = read_json(data_path)
    full_obj["analysis_summary"] = summary
    data_obj["analysis_summary"] = summary
    write_json(full_path, full_obj)
    write_json(data_path, data_obj)


def update_registry(summary: dict):
    path = DATA_DIR / "table_registry.json"
    obj = read_json(path)
    tables = obj.get("tables", [])
    text = summary.get("text", "")
    for t in tables:
        if t.get("id") in TARGET_TABLES:
            base = "Master Table 手工保存结果"
            source = str(t.get("description", "") or "")
            if "来源：" in source:
                base = source.split("；")[0].split("。")[0]
            t["description"] = f"{base}；{text}"
    write_json(path, obj)


def main():
    summary = build_summary()
    for table_id in TARGET_TABLES:
        attach_to_table(table_id, summary)
    update_registry(summary)
    print("attached", TARGET_TABLES)


if __name__ == "__main__":
    main()
