#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "backend" / "data"
REPORT_CSV = PROJECT_DIR / "backend" / "runtime" / "reports" / "workflow_table_2_caver_pocket_seed_v1_acetophenone_summary.csv"

TARGET_TABLES = [
    "workflow_table_3",
    "workflow_table_5",
]
REMOVED_COLUMNS = [
    "苯乙酮_活性口袋命中Zn",
    "Zn第一配位壳最近空隙_A",
    "ZN_弯曲通道半径_A",
    "ZN_弯曲通道长度_A",
]
DISPLAY_COLUMNS = [
    "苯乙酮_CAVER通道判定",
    "苯乙酮_CAVER口袋深度",
    "苯乙酮_CAVER口袋深度判定列",
    "苯乙酮_CAVER口袋深度判定值",
]


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload):
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def load_report():
    with REPORT_CSV.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        rows = []
        for raw in csv.DictReader(f):
            rows.append({str(k or "").strip(): v for k, v in raw.items()})
    return {r["PDB_ID"].strip().upper(): r for r in rows if str(r.get("PDB_ID", "")).strip()}


def row_id(row):
    return str(row.get("Representative") or row.get("PDB_ID") or row.get("_id") or "").strip().upper()


def parse_depth_basis(text):
    raw = str(text or "").strip()
    if not raw:
        return "", ""
    match = re.match(r"^([A-Za-z0-9_]+)\s*=\s*(.+)$", raw)
    if not match:
        return "", raw
    return match.group(1).strip(), match.group(2).strip()


def apply_caver_display_fields(row, depth_label, depth_basis):
    row["苯乙酮_CAVER通道判定"] = row.get("CAVER_AcetophenonePass", "")
    row["苯乙酮_CAVER口袋深度"] = depth_label
    basis_col, basis_value = parse_depth_basis(depth_basis)
    row["苯乙酮_CAVER口袋深度判定列"] = basis_col
    row["苯乙酮_CAVER口袋深度判定值"] = basis_value
    return row


def update_row(row, rep):
    depth_label, depth_basis = classify_depth(row)
    row["苯乙酮_综合判断"] = (
        "可以" if rep.get("Acetophenone_size_pass_by_caver") == "yes"
        else ("边缘" if rep.get("Acetophenone_size_pass_by_caver") == "borderline" else "不支持")
    )
    row["苯乙酮_尺寸可进入"] = rep.get("Acetophenone_size_pass_by_caver", "")
    row["苯乙酮_通道瓶颈半径_A"] = rep.get("Best_bottleneck_radius_A", "")
    if str(row.get("苯乙酮_口袋最小半径_A", "") or "").strip() == "":
        row["苯乙酮_口袋最小半径_A"] = rep.get("Best_bottleneck_radius_A", "")
    row["CAVER_SeedMode"] = rep.get("Caver_seed_mode", "")
    row["CAVER_ChannelCount"] = rep.get("Channel_count", "")
    row["CAVER_BestBottleneckRadius_A"] = rep.get("Best_bottleneck_radius_A", "")
    row["CAVER_BestThroughput"] = rep.get("Best_throughput", "")
    row["CAVER_AcetophenonePass"] = rep.get("Acetophenone_size_pass_by_caver", "")
    row["CAVER_SeedXYZ"] = rep.get("Seed_xyz", "")
    row["CAVER_JobDir"] = rep.get("Job_dir", "")
    row["CAVER_PocketDepthLabel"] = depth_label
    row["CAVER_PocketDepthBasis"] = depth_basis
    apply_caver_display_fields(row, depth_label, depth_basis)
    row["苯乙酮_判断说明"] = (
        f"CAVER pocket-seed 重算：{rep.get('Caver_seed_mode','')}; "
        f"channels={rep.get('Channel_count','')}。详细判定见专用列。"
        + (f" 原备注：{rep.get('Existing_note','')}" if str(rep.get('Existing_note','')).strip() else "")
    ).strip()
    for col in REMOVED_COLUMNS:
        row.pop(col, None)
    return row


def data_item_from_row(row):
    return {
        "id": row_id(row),
        "name": str(row.get("Protein_Name", "") or ""),
        "cluster": str(row.get("Cluster_ID", "") or ""),
        "receptor_pdb": str(row.get("Receptor_PDB", "") or ""),
        "receptor_rel": str(row.get("Receptor_PDB", "") or ""),
        "best_sdf": str(row.get("Best_SDF", "") or ""),
        "species": str(row.get("Species", "") or ""),
        "monomer_seq": str(row.get("MonomerSeq", "") or ""),
        "residue_length": row.get("Zn_Binding_MonomerResidueCount", ""),
        "zn_coord_site": str(row.get("Zn_CoordSite", row.get("Zn_CoordSite", "")) or ""),
        "best_pocket_id": str(row.get("BestPocket_ID", "") or ""),
        "best_pocket_min_dist_to_zn": str(row.get("BestPocket_MinDistToZn", "") or ""),
        "acetophenone_judgement": str(row.get("苯乙酮_综合判断", "") or ""),
        "caver_best_bottleneck_radius_a": str(row.get("CAVER_BestBottleneckRadius_A", "") or ""),
        "caver_pocket_depth_label": str(row.get("CAVER_PocketDepthLabel", "") or ""),
    }


def classify_depth(row):
    depth = to_float(row.get("ZN_Depth", ""))
    surface = to_float(row.get("ZN_Surface_Distance", ""))
    if surface is not None and surface <= 0.5:
        return "表面暴露", f"ZN_Surface_Distance={surface:.3f} Å"
    if depth is not None and depth <= 2.0:
        return "表面暴露", f"ZN_Depth={depth:.3f}"
    if surface is not None and surface <= 2.5:
        return "半埋藏", f"ZN_Surface_Distance={surface:.3f} Å"
    if depth is not None and depth <= 4.5:
        return "半埋藏", f"ZN_Depth={depth:.3f}"
    if depth is not None:
        return "深埋藏", f"ZN_Depth={depth:.3f}"
    if surface is not None:
        return "深埋藏", f"ZN_Surface_Distance={surface:.3f} Å"
    return "未知", "缺少 ZN_Depth / ZN_Surface_Distance"


def to_float(value):
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except Exception:
        return None


def backfill_table(table_id: str, report_by_pid):
    full_path = DATA_DIR / f"{table_id}_master_full.json"
    data_path = DATA_DIR / f"{table_id}_data.json"
    csv_path = DATA_DIR / f"{table_id}_master_table.csv"
    ids_path = DATA_DIR / f"{table_id}_ids.json"

    full_obj = read_json(full_path)
    rows = list(full_obj.get("rows") or [])
    header = [col for col in list(full_obj.get("header") or []) if col not in REMOVED_COLUMNS]
    extra_cols = [
        "CAVER_SeedMode",
        "CAVER_ChannelCount",
        "CAVER_BestBottleneckRadius_A",
        "CAVER_BestThroughput",
        "CAVER_AcetophenonePass",
        "CAVER_PocketDepthLabel",
        "CAVER_PocketDepthBasis",
        "CAVER_SeedXYZ",
        "CAVER_JobDir",
        *DISPLAY_COLUMNS,
    ]
    for col in extra_cols:
        if col not in header:
            header.append(col)

    updated = 0
    out_rows = []
    for row in rows:
        row = dict(row)
        for col in REMOVED_COLUMNS:
            row.pop(col, None)
        pid = row_id(row)
        rep = report_by_pid.get(pid)
        if rep:
            row = update_row(row, rep)
            updated += 1
        out_rows.append(row)

    write_json(full_path, {"header": header, "rows": out_rows})

    with csv_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in out_rows:
            writer.writerow([row.get(col, "") for col in header])

    ids = [row_id(r) for r in out_rows if row_id(r)]
    write_json(ids_path, {"table_id": table_id, "ids": ids, "count": len(ids), "label": read_json(ids_path).get("label", table_id) if ids_path.exists() else table_id})
    write_json(data_path, {"items": [data_item_from_row(r) for r in out_rows if row_id(r)]})
    return updated, len(out_rows)


def main():
    report_by_pid = load_report()
    for table_id in TARGET_TABLES:
        updated, total = backfill_table(table_id, report_by_pid)
        print(table_id, "updated", updated, "total", total)


if __name__ == "__main__":
    main()
