import csv
import json
import os
from pathlib import Path

from diffdock_api_server import _pdbzn_dist, _pdbzn_find_tmalign_pdb_path, _pdbzn_structure_atom_rows

BASE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = BASE_DIR.parent
DATA_DIR = BASE_DIR / "data"

HIS_RESIDUES = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}
PROTEIN_DONOR_ATOMS = {
    "ALA": set(),
    "ARG": {"NE", "NH1", "NH2"},
    "ASN": {"OD1", "ND2"},
    "ASP": {"OD1", "OD2"},
    "ASH": {"OD1", "OD2"},
    "CYS": {"SG"},
    "CYM": {"SG"},
    "CYX": {"SG"},
    "GLN": {"OE1", "NE2"},
    "GLU": {"OE1", "OE2"},
    "GLH": {"OE1", "OE2"},
    "GLY": set(),
    "HIS": {"ND1", "NE2"},
    "HID": {"ND1", "NE2"},
    "HIE": {"ND1", "NE2"},
    "HIP": {"ND1", "NE2"},
    "HSD": {"ND1", "NE2"},
    "HSE": {"ND1", "NE2"},
    "HSP": {"ND1", "NE2"},
    "ILE": set(),
    "LEU": set(),
    "LYS": {"NZ"},
    "LYN": {"NZ"},
    "MET": {"SD"},
    "MSE": {"SE"},
    "PHE": set(),
    "PRO": set(),
    "SEC": {"SE"},
    "SER": {"OG"},
    "THR": {"OG1"},
    "TRP": {"NE1"},
    "TYR": {"OH"},
    "VAL": set(),
}

COORD_COLUMN_RESIDUE = "Zn_CoordResidueCount"
COORD_COLUMN_HIS = "Zn_CoordHisCount"
COORD_COLUMN_NON_HIS = "Zn_CoordNonHisCount"
COORD_COLUMN_RESIDUES = "Zn_CoordResidues"
COORD_COLUMN_SITE = "Zn_CoordSite"
DROP_PREFIXES = ("STEP5_",)

def load_master():
    env = str(os.environ.get("PDBZN_MASTER_TABLE", "") or "").strip()
    paths = []
    if env:
        paths.append(Path(env).expanduser())
    paths.extend(
        [
            DATA_DIR / "master_table.csv",
            PROJECT_DIR / "PDB_ZN" / "zn_his_master_table7.csv",
            Path("/root/PDB_ZN/zn_his_master_table7.csv"),
            Path("/root/zn_his_master_table.csv"),
        ]
    )
    for p in paths:
        if p.exists():
            with p.open(newline="", encoding="utf-8", errors="ignore") as f:
                rows = list(csv.reader(f))
            return str(p), rows
    raise SystemExit("master table not found")

def parse_lengths(s):
    if not s:
        return None, None, None
    parts = s.split(";")
    if len(parts) != 3:
        return None, None, None
    try:
        return int(parts[0]), int(parts[1]), int(parts[2])
    except:
        return parts[0], parts[1], parts[2]


def residue_key(atom):
    return (
        str((atom or {}).get("chain", "") or "").strip(),
        str((atom or {}).get("seq", "") or "").strip(),
        str((atom or {}).get("icode", "") or "").strip(),
        str((atom or {}).get("comp", "") or "").strip().upper(),
    )


def residue_label(res_key):
    chain, seq, icode, comp = res_key
    c = chain or "?"
    s = seq or "?"
    ic = "" if (not icode or icode in {"?", "."}) else icode
    return f"{comp}:{c}:{s}{ic}"


def donor_atom_allowed(atom):
    comp = str((atom or {}).get("comp", "") or "").strip().upper()
    atom_name = str((atom or {}).get("atom", "") or "").strip().upper()
    return atom_name in PROTEIN_DONOR_ATOMS.get(comp, set())


def select_best_zn_site(atoms):
    zn_atoms = [a for a in atoms if str(a.get("comp") or "").upper() == "ZN" or str(a.get("element") or "").upper() == "ZN"]
    donor_atoms = [a for a in atoms if donor_atom_allowed(a)]
    best = None
    for zn in zn_atoms:
        nearest = {}
        for atom in donor_atoms:
            d = _pdbzn_dist(zn, atom)
            if d > 3.2:
                continue
            rk = residue_key(atom)
            old = nearest.get(rk)
            if old is None or d < old["distance"]:
                nearest[rk] = {"atom": atom, "distance": float(d)}
        his_count = sum(1 for rk in nearest if rk[3] in HIS_RESIDUES)
        total_count = len(nearest)
        dist_sum = sum(sorted(x["distance"] for x in nearest.values())[:3]) if nearest else 9999.0
        candidate = {
            "zn": zn,
            "nearest": nearest,
            "his_count": his_count,
            "total_count": total_count,
            "dist_sum": dist_sum,
        }
        if best is None:
            best = candidate
            continue
        better = (
            candidate["his_count"] > best["his_count"]
            or (
                candidate["his_count"] == best["his_count"]
                and (
                    candidate["total_count"] > best["total_count"]
                    or (
                        candidate["total_count"] == best["total_count"]
                        and candidate["dist_sum"] < best["dist_sum"]
                    )
                )
            )
        )
        if better:
            best = candidate
    return best


def coordination_metrics(rep):
    fp = _pdbzn_find_tmalign_pdb_path(rep)
    if not fp:
        return {
            COORD_COLUMN_RESIDUE: 0,
            COORD_COLUMN_HIS: 0,
            COORD_COLUMN_NON_HIS: 0,
            COORD_COLUMN_RESIDUES: "",
            COORD_COLUMN_SITE: "",
        }
    atoms = _pdbzn_structure_atom_rows(fp)
    if not atoms:
        return {
            COORD_COLUMN_RESIDUE: 0,
            COORD_COLUMN_HIS: 0,
            COORD_COLUMN_NON_HIS: 0,
            COORD_COLUMN_RESIDUES: "",
            COORD_COLUMN_SITE: "",
        }
    best = select_best_zn_site(atoms)
    if not best:
        return {
            COORD_COLUMN_RESIDUE: 0,
            COORD_COLUMN_HIS: 0,
            COORD_COLUMN_NON_HIS: 0,
            COORD_COLUMN_RESIDUES: "",
            COORD_COLUMN_SITE: "",
        }
    nearest = best["nearest"]
    ordered = sorted(nearest.items(), key=lambda kv: kv[1]["distance"])
    residue_count = len(ordered)
    his_count = sum(1 for rk, _ in ordered if rk[3] in HIS_RESIDUES)
    residue_labels = [f"{residue_label(rk)}:{round(item['distance'], 3)}" for rk, item in ordered]
    zn = best["zn"]
    return {
        COORD_COLUMN_RESIDUE: residue_count,
        COORD_COLUMN_HIS: his_count,
        COORD_COLUMN_NON_HIS: max(0, residue_count - his_count),
        COORD_COLUMN_RESIDUES: "; ".join(residue_labels),
        COORD_COLUMN_SITE: f"{str(zn.get('chain') or '?')}:{str(zn.get('seq') or '?')}",
    }


def write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)


def write_csv(path: Path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in header})


def load_table1_ids():
    p = DATA_DIR / "table1_ids.json"
    if not p.exists():
        return []
    try:
        obj = json.load(p.open("r", encoding="utf-8"))
    except Exception:
        return []
    ids = []
    for x in obj.get("ids", []):
        s = str(x or "").strip().upper()
        if s and s not in ids:
            ids.append(s)
    return ids


def load_registry():
    p = DATA_DIR / "table_registry.json"
    if not p.exists():
        return []
    try:
        obj = json.load(p.open("r", encoding="utf-8"))
    except Exception:
        return []
    tables = obj.get("tables", [])
    return [t for t in tables if isinstance(t, dict)]


def save_registry(tables):
    def sort_key(entry):
        table_id = str((entry or {}).get("id", "") or "")
        digits = "".join(ch for ch in table_id if ch.isdigit())
        order = int(digits) if digits else 10**9
        return (order, table_id)
    write_json(DATA_DIR / "table_registry.json", {"tables": sorted(tables, key=sort_key)})


def keep_column(name):
    raw = str(name or "").strip()
    if not raw or raw == "_id":
        return False
    upper = raw.upper()
    for prefix in DROP_PREFIXES:
        if upper.startswith(prefix):
            return False
    return True

def main():
    master_path, rows = load_master()
    raw_header = rows[0]
    header = [name for name in raw_header if keep_column(name)]
    for extra in [COORD_COLUMN_RESIDUE, COORD_COLUMN_HIS, COORD_COLUMN_NON_HIS, COORD_COLUMN_RESIDUES, COORD_COLUMN_SITE]:
        if extra not in header:
            header.append(extra)
    idx = {name: i for i, name in enumerate(raw_header)}
    out_all = []
    full_all = []
    for r in rows[1:]:
        def get(name):
            i = idx.get(name)
            return r[i] if i is not None and i < len(r) else ""
        rep = get("Representative") or os.path.splitext(os.path.basename(get("Receptor_PDB")))[0]
        coord = coordination_metrics(rep)
        protein_name = get("Protein_Name")
        cluster = get("Cluster_ID")
        receptor_pdb = get("Receptor_PDB")
        best_sdf = get("Best_SDF")
        species = get("Species")
        monomer_seq = get("MonomerSeq")
        length_field = get("Length;Oligomer;Monomer")
        residue_len, oligomer, monomer_len = parse_lengths(length_field)
        # build relative receptor path under server root if file exists
        ext = ".pdb"
        if receptor_pdb and receptor_pdb.lower().endswith(".cif"):
            ext = ".cif"
        receptor_rel = f"../backend/data/structures/{rep}{ext}"
        item = {
            "id": rep,
            "name": protein_name,
            "cluster": cluster,
            "receptor_pdb": receptor_pdb,
            "receptor_rel": receptor_rel,
            "best_sdf": best_sdf,
            "species": species,
            "monomer_seq": monomer_seq,
            "residue_length": residue_len,
            "oligomer": oligomer,
            "monomer_length": monomer_len,
            "zn_coord_residue_count": coord[COORD_COLUMN_RESIDUE],
            "zn_coord_his_count": coord[COORD_COLUMN_HIS],
            "zn_coord_non_his_count": coord[COORD_COLUMN_NON_HIS],
        }
        out_all.append(item)
        obj = {}
        rep_val = ""
        rp_val = ""
        rep_idx = idx.get("Representative")
        if rep_idx is not None and rep_idx < len(r):
            rep_val = r[rep_idx]
        rp_idx = idx.get("Receptor_PDB")
        if rp_idx is not None and rp_idx < len(r):
            rp_val = r[rp_idx]
        for k in header:
            i = idx.get(k)
            if k == "BestPocket_ID" and "_id" not in obj:
                if rep_val:
                    obj["_id"] = rep_val
                elif rp_val:
                    obj["_id"] = os.path.splitext(os.path.basename(rp_val))[0]
            if k in coord:
                obj[k] = coord[k]
            else:
                obj[k] = r[i] if i is not None and i < len(r) else ""
        if "_id" not in obj:
            if rep_val:
                obj["_id"] = rep_val
            elif rp_val:
                obj["_id"] = os.path.splitext(os.path.basename(rp_val))[0]
        full_all.append(obj)
    table1_ids = load_table1_ids()
    if table1_ids:
        want = set(table1_ids)
        item_map = {}
        row_map = {}
        for item in out_all:
            pid = str(item.get("id", "")).strip().upper()
            if pid and pid not in item_map:
                item_map[pid] = item
        for row in full_all:
            pid = str(row.get("_id", "") or row.get("Representative", "")).strip().upper()
            if pid and pid not in row_map:
                row_map[pid] = row
        out = [item_map[pid] for pid in table1_ids if pid in item_map]
        full = [row_map[pid] for pid in table1_ids if pid in row_map]
    else:
        out = [item for item in out_all if str(item.get("id", "")).strip()]
        full = [row for row in full_all if str(row.get("_id", "")).strip()]
    data_payload = {"items": out}
    full_payload = {"rows": full, "header": header}
    write_json(DATA_DIR / "data.json", data_payload)
    print("OK", DATA_DIR / "data.json", len(out))
    write_json(DATA_DIR / "master_full.json", full_payload)
    print("OK", DATA_DIR / "master_full.json", len(full))
    write_json(DATA_DIR / "table1_data.json", data_payload)
    print("OK", DATA_DIR / "table1_data.json", len(out))
    write_json(DATA_DIR / "table1_master_full.json", full_payload)
    print("OK", DATA_DIR / "table1_master_full.json", len(full))
    write_csv(DATA_DIR / "table1_master_table.csv", ["_id"] + header, full)
    print("OK", DATA_DIR / "table1_master_table.csv", len(full))
    existing = [t for t in load_registry() if str((t or {}).get("id", "") or "") != "table1"]
    existing.append(
        {
            "id": "table1",
            "label": f"Table 1 ({len(full)} proteins)",
            "full": "table1_master_full.json",
            "data": "table1_data.json",
            "csv": "table1_master_table.csv",
            "ids": "table1_ids.json",
            "row_count": len(full),
            "description": "主表快照（Table 1）",
        }
    )
    save_registry(existing)
    print("OK", DATA_DIR / "table_registry.json", len(existing))

if __name__ == "__main__":
    main()
