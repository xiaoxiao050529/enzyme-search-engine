#!/usr/bin/env python3
"""
Annotate a master-table bundle with acetophenone passage heuristics.

This version intentionally avoids fpocket alpha-sphere minima as the primary
"radius" metric because those values are strongly compressed by fpocket's own
radius floor and can be misleading when interpreted as entrance radii.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
BACKEND_DIR = PROJECT_DIR / "backend"
DATA_DIR = BACKEND_DIR / "data"
sys.path.insert(0, str(BACKEND_DIR))

from diffdock_api_server import _pdbzn_structure_atom_rows  # noqa: E402


REQUIRED_RADIUS = 2.94
REMOVED_COLUMNS = [
    "苯乙酮_活性口袋命中Zn",
    "Zn第一配位壳最近空隙_A",
    "ZN_弯曲通道半径_A",
    "ZN_弯曲通道长度_A",
]
NEW_COLUMNS = [
    "苯乙酮_综合判断",
    "苯乙酮_尺寸可进入",
    "苯乙酮_O到Zn最短距离_A",
    "苯乙酮_任意原子到Zn最短距离_A",
    "苯乙酮_判断说明",
]
FRONT_COLUMNS = [
    "Representative",
    "苯乙酮_综合判断",
    "苯乙酮_尺寸可进入",
    "苯乙酮_O到Zn最短距离_A",
    "苯乙酮_任意原子到Zn最短距离_A",
    "苯乙酮_判断说明",
    "Zn_CoordResidueCount",
    "Zn_CoordHisCount",
    "Zn_CoordNonHisCount",
    "Zn_CoordSite",
    "Zn_CoordResidues",
    "ZN_Depth",
    "ZN_Depth_Rounded",
    "ZN_SASA",
    "ZN_Surface_Distance",
]

ZH_LABELS = {
    "Representative": "PDB编号",
    "苯乙酮_综合判断": "苯乙酮综合判断",
    "苯乙酮_尺寸可进入": "苯乙酮尺寸可进入",
    "苯乙酮_O到Zn最短距离_A": "苯乙酮O到Zn最短距离(Å)",
    "苯乙酮_任意原子到Zn最短距离_A": "苯乙酮任意原子到Zn最短距离(Å)",
    "苯乙酮_判断说明": "苯乙酮判断说明",
    "Zn_CoordResidueCount": "Zn配位残基总数",
    "Zn_CoordHisCount": "Zn配位His数",
    "Zn_CoordNonHisCount": "Zn配位非His数",
    "Zn_CoordSite": "Zn配位位点",
    "Zn_CoordResidues": "Zn配位残基明细",
    "ZN_Depth": "Zn埋藏深度",
    "ZN_Depth_Rounded": "Zn埋藏深度(分箱)",
    "ZN_SASA": "Zn溶剂可及面积",
    "ZN_Surface_Distance": "Zn到表面最短直线逃逸距离",
}

VDW = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "SE": 1.90,
    "ZN": 1.39,
}
WATER_RESIDUES = {"HOH", "WAT", "DOD"}
LOCAL_MARGIN = 18.0
PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
    "GLY", "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS",
    "LYN", "MET", "MSE", "PHE", "PRO", "SEC", "SER", "THR", "TRP", "TYR", "VAL",
}


def as_text(value) -> str:
    return "" if value is None else str(value)


def fmt_num(value) -> str:
    if value is None:
        return ""
    return f"{float(value):.3f}"


def dist(a, b) -> float:
    dx = float(a[0]) - float(b[0])
    dy = float(a[1]) - float(b[1])
    dz = float(a[2]) - float(b[2])
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def reorder_keys(existing):
    keys = [as_text(x).strip() for x in existing if as_text(x).strip() and as_text(x).strip() not in REMOVED_COLUMNS]
    out = []
    seen = set()
    for key in FRONT_COLUMNS + keys:
        if not key or key in seen or key not in keys:
            continue
        seen.add(key)
        out.append(key)
    return out


def parse_sdf_atoms(path: Path):
    lines = path.read_text(errors="ignore").splitlines()
    if len(lines) < 4:
        return []
    try:
        atom_count = int(lines[3][:3])
    except Exception:
        return []
    atoms = []
    for line in lines[4 : 4 + atom_count]:
        if len(line) < 34:
            continue
        try:
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
            element = line[31:34].strip().upper()
        except Exception:
            continue
        atoms.append((element, (x, y, z)))
    return atoms


def first_local_structure(pid: str):
    upper = as_text(pid).strip().upper()
    for cand in [
        DATA_DIR / "structures" / f"{upper}.pdb",
        DATA_DIR / "structures" / f"{upper}.cif",
        DATA_DIR / "structures" / f"{upper.lower()}.pdb",
        DATA_DIR / "structures" / f"{upper.lower()}.cif",
    ]:
        if cand.exists():
            return cand
    return None


def find_target_zn(row, structure_path: Path):
    atoms = _pdbzn_structure_atom_rows(structure_path)
    site = as_text((row or {}).get("Zn_CoordSite")).strip()
    want_chain = ""
    want_seq = ""
    if ":" in site:
        want_chain, want_seq = site.split(":", 1)
    for atom in atoms:
        if as_text(atom.get("element")).upper() != "ZN":
            continue
        if want_chain and as_text(atom.get("chain")).strip() != want_chain:
            continue
        if want_seq and as_text(atom.get("seq")).strip() != want_seq:
            continue
        return atom, atoms
    for atom in atoms:
        if as_text(atom.get("element")).upper() == "ZN":
            return atom, atoms
    return None, atoms


def best_pose_distances(pid: str, zn_point):
    pose_dir = DATA_DIR / "diffdock" / as_text(pid).strip().upper()
    if zn_point is None or not pose_dir.exists():
        return None, None
    best_any = None
    best_o = None
    for pose_path in sorted(pose_dir.glob("rank*.sdf")):
        atoms = parse_sdf_atoms(pose_path)
        if not atoms:
            continue
        any_min = min(dist(point, zn_point) for _, point in atoms)
        oxygen_dists = [dist(point, zn_point) for element, point in atoms if element == "O"]
        o_min = min(oxygen_dists) if oxygen_dists else None
        if best_any is None or any_min < best_any:
            best_any = any_min
        if o_min is not None and (best_o is None or o_min < best_o):
            best_o = o_min
    return best_any, best_o


def protein_atoms_same_chain(atoms, chain_id: str):
    out = []
    for atom in atoms:
        if chain_id and as_text(atom.get("chain")).strip() != chain_id:
            continue
        comp = as_text(atom.get("comp")).strip().upper()
        if comp in WATER_RESIDUES or comp not in PROTEIN_RESIDUES:
            continue
        element = as_text(atom.get("element")).strip().upper()
        if element == "H":
            continue
        out.append(atom)
    return out


def local_surface_clearances(atoms, zn_atom):
    zn_point = (float(zn_atom["x"]), float(zn_atom["y"]), float(zn_atom["z"]))
    out = []
    for atom in atoms:
        if atom is zn_atom:
            continue
        center = (float(atom["x"]), float(atom["y"]), float(atom["z"]))
        if max(abs(center[i] - zn_point[i]) for i in range(3)) > LOCAL_MARGIN:
            continue
        radius = VDW.get(as_text(atom.get("element")).upper(), 1.70)
        out.append(max(0.0, dist(center, zn_point) - radius))
    out.sort()
    return out


def pocket_match(row) -> bool:
    return as_text((row or {}).get("BestPocket_ZnMatch")).strip() == "1"


def pocket_min_radius(row) -> float | None:
    # The old fpocket-derived "minimum radius" is intentionally not reused.
    # Return None for now rather than a misleading 3.4A floor-truncated value.
    return None


def size_pass(entrance_radius: float | None, pocket_has_zn: bool) -> str:
    if not pocket_has_zn:
        return "unknown"
    if entrance_radius is None:
        return ""
    return "yes" if entrance_radius >= REQUIRED_RADIUS else "no"


def assessment(size_flag: str, pocket_has_zn: bool, best_any, best_o):
    if not pocket_has_zn:
        if best_o is not None and best_o <= 4.0:
            return "待复核", "best_pocket 未命中 Zn，但现有 pose 已把羰基 O 送到 Zn 近邻"
        if best_any is not None and best_any <= 4.0:
            return "待复核", "best_pocket 未命中 Zn，但现有 pose 已进入 Zn 邻域，建议单独复核真实活性口袋"
        return "偏不支持", "当前 best_pocket 未命中 Zn 活性位点，且现有苯乙酮 pose 未显示稳定靠近 Zn"

    if size_flag == "no":
        return "不顺利", "按 Zn 邻域最近原子表面空隙估计，活性位点入口第一层空间偏紧"

    if best_o is not None and best_o <= 3.70:
        return "可以", "活性位点附近存在可接近 Zn 的 pose，且羰基 O 已靠近 Zn（<= 3.7 Å）"
    if best_o is not None and best_o <= 4.00:
        return "大概率可以", "活性位点附近存在可接近 Zn 的 pose，且羰基 O 已进入可接受 Zn 邻域（<= 4.0 Å）"
    if best_any is not None and best_any <= 4.00:
        return "边缘", "底物整体已能进入 Zn 邻域，但羰基 O 朝向仍未对准 Zn"
    return "不顺利", "现有 pose 仍未把羰基 O 或主体稳定送到 Zn 位点"


def build_updates(row):
    pid = as_text((row or {}).get("Representative") or (row or {}).get("_id")).strip().upper()
    if not pid:
        return {}

    structure_path = first_local_structure(pid)
    if structure_path is None:
        return {
            "苯乙酮_综合判断": "",
            "苯乙酮_尺寸可进入": "",
            "苯乙酮_O到Zn最短距离_A": "",
            "苯乙酮_任意原子到Zn最短距离_A": "",
            "苯乙酮_判断说明": "缺少本地结构文件，未计算",
        }

    zn_atom, atoms = find_target_zn(row, structure_path)
    zn_point = None if zn_atom is None else (float(zn_atom["x"]), float(zn_atom["y"]), float(zn_atom["z"]))
    best_any, best_o = best_pose_distances(pid, zn_point)
    chain_atoms = protein_atoms_same_chain(atoms, as_text(zn_atom.get("chain")).strip() if zn_atom else "")
    local_clearances = local_surface_clearances(chain_atoms, zn_atom) if zn_atom else []
    entrance_radius = local_clearances[0] if local_clearances else None
    p_has_zn = pocket_match(row)
    p_min = None
    sflag = size_pass(entrance_radius, p_has_zn)
    grade, note = assessment(sflag, p_has_zn, best_any, best_o)
    return {
        "苯乙酮_综合判断": grade,
        "苯乙酮_尺寸可进入": sflag,
        "苯乙酮_O到Zn最短距离_A": fmt_num(best_o),
        "苯乙酮_任意原子到Zn最短距离_A": fmt_num(best_any),
        "苯乙酮_判断说明": note,
    }


def update_master_full(path: Path, updates_by_pid):
    obj = json.loads(path.read_text(encoding="utf-8"))
    rows = list(obj.get("rows") or [])
    for row in rows:
        pid = as_text((row or {}).get("Representative") or (row or {}).get("_id")).strip().upper()
        for key, value in updates_by_pid.get(pid, {}).items():
            row[key] = value
        row.pop("中文名", None)
        for key in REMOVED_COLUMNS:
            row.pop(key, None)
    header = list(obj.get("header") or [])
    header = [h for h in header if h != "中文名" and h not in REMOVED_COLUMNS]
    for key in NEW_COLUMNS:
        if key not in header:
            header.append(key)
    header = reorder_keys(header)
    obj["header"] = header
    obj["rows"] = [{k: row.get(k, "") for k in header} for row in rows]
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def update_csv_table(path: Path, updates_by_pid):
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
        fieldnames = list(rows[0].keys()) if rows else []
    for row in rows:
        pid = as_text((row or {}).get("Representative") or (row or {}).get("_id")).strip().upper()
        for key, value in updates_by_pid.get(pid, {}).items():
            row[key] = value
        row.pop("中文名", None)
        for key in REMOVED_COLUMNS:
            row.pop(key, None)
    fieldnames = [f for f in fieldnames if f != "中文名" and f not in REMOVED_COLUMNS]
    for key in NEW_COLUMNS:
        if key not in fieldnames:
            fieldnames.append(key)
    fieldnames = reorder_keys(fieldnames)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def update_data_json(path: Path, updates_by_pid):
    obj = json.loads(path.read_text(encoding="utf-8"))
    for item in obj.get("items", []):
        pid = as_text((item or {}).get("id")).strip().upper()
        for key, value in updates_by_pid.get(pid, {}).items():
            item[key] = value
        item.pop("中文名", None)
        for key in REMOVED_COLUMNS:
            item.pop(key, None)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table-id", default="workflow_table", help="Table id/prefix to update")
    parser.add_argument("--ids", default="", help="Optional comma-separated Representatives to update")
    parser.add_argument("--channel-json", default="", help="Optional JSON results file from analyze_zn_surface_channels.py")
    args = parser.parse_args()

    prefix = as_text(args.table_id).strip()
    if not prefix:
        raise SystemExit("table id is required")

    full_path = DATA_DIR / f"{prefix}_master_full.json"
    csv_path = DATA_DIR / f"{prefix}_master_table.csv"
    data_path = DATA_DIR / f"{prefix}_data.json"
    if not full_path.exists() or not csv_path.exists() or not data_path.exists():
        raise SystemExit(f"missing one of: {full_path.name}, {csv_path.name}, {data_path.name}")

    obj = json.loads(full_path.read_text(encoding="utf-8"))
    rows = list(obj.get("rows") or [])
    wanted = {x.strip().upper() for x in as_text(args.ids).split(",") if x.strip()}
    channel_map = {}
    if args.channel_json:
        raw = json.loads(Path(args.channel_json).read_text(encoding="utf-8"))
        for item in raw:
            pid = as_text(item.get("pid")).strip().upper()
            if pid:
                channel_map[pid] = item
    updates_by_pid = {}
    for row in rows:
        pid = as_text((row or {}).get("Representative") or (row or {}).get("_id")).strip().upper()
        if not pid or pid in updates_by_pid:
            continue
        if wanted and pid not in wanted:
            continue
        updates_by_pid[pid] = build_updates(row)

    update_master_full(full_path, updates_by_pid)
    update_csv_table(csv_path, updates_by_pid)
    update_data_json(data_path, updates_by_pid)
    print(f"updated {prefix}: {len(updates_by_pid)} rows")


if __name__ == "__main__":
    main()
