import csv
import json
import os
import re
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
COORD_COLUMN_NOTE = "Zn_CoordNote"
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


def site_sort_key(site):
    return (
        -int(site.get("his_count", 0)),
        int(site.get("non_his_count", 10**9)),
        int(site.get("residue_count", 10**9)),
        float(site.get("dist_sum", 9999.0)),
        str(site.get("site", "")),
        int(site.get("index", 10**9)),
    )


def summarize_site(site):
    if not site:
        return ""
    return (
        f"Zn{site.get('index', '?')}@{site.get('site', '?')}: "
        f"total={site.get('residue_count', 0)}, "
        f"his={site.get('his_count', 0)}, "
        f"non_his={site.get('non_his_count', 0)}, "
        f"residues={site.get('residues_text', '') or '-'}"
    )


def collect_zn_sites(atoms):
    zn_atoms = [a for a in atoms if str(a.get("comp") or "").upper() == "ZN" or str(a.get("element") or "").upper() == "ZN"]
    donor_atoms = [a for a in atoms if donor_atom_allowed(a)]
    sites = []
    for idx, zn in enumerate(zn_atoms, start=1):
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
        ordered = sorted(nearest.items(), key=lambda kv: (kv[1]["distance"], residue_label(kv[0])))
        residue_labels = [f"{residue_label(rk)}:{round(item['distance'], 3)}" for rk, item in ordered]
        site = {
            "index": idx,
            "zn": zn,
            "nearest": nearest,
            "ordered": ordered,
            "his_count": his_count,
            "residue_count": total_count,
            "non_his_count": max(0, total_count - his_count),
            "dist_sum": dist_sum,
            "site": f"{str(zn.get('chain') or '?')}:{str(zn.get('seq') or '?')}",
            "residues_text": "; ".join(residue_labels),
        }
        sites.append(site)
    return sorted(sites, key=site_sort_key)


def coordination_metrics(rep):
    fp = _pdbzn_find_tmalign_pdb_path(rep)
    if not fp:
        return {
            COORD_COLUMN_RESIDUE: 0,
            COORD_COLUMN_HIS: 0,
            COORD_COLUMN_NON_HIS: 0,
            COORD_COLUMN_RESIDUES: "",
            COORD_COLUMN_SITE: "",
            COORD_COLUMN_NOTE: "",
        }
    atoms = _pdbzn_structure_atom_rows(fp)
    if not atoms:
        return {
            COORD_COLUMN_RESIDUE: 0,
            COORD_COLUMN_HIS: 0,
            COORD_COLUMN_NON_HIS: 0,
            COORD_COLUMN_RESIDUES: "",
            COORD_COLUMN_SITE: "",
            COORD_COLUMN_NOTE: "",
        }
    sites = collect_zn_sites(atoms)
    if not sites:
        return {
            COORD_COLUMN_RESIDUE: 0,
            COORD_COLUMN_HIS: 0,
            COORD_COLUMN_NON_HIS: 0,
            COORD_COLUMN_RESIDUES: "",
            COORD_COLUMN_SITE: "",
            COORD_COLUMN_NOTE: "",
        }
    is_multi_zn = len(sites) > 1
    if is_multi_zn:
        qualified = [s for s in sites if int(s.get("his_count", 0)) >= 3]
        if qualified:
            best = sorted(
                qualified,
                key=lambda s: (
                    int(s.get("non_his_count", 10**9)),
                    -int(s.get("his_count", 0)),
                    int(s.get("residue_count", 10**9)),
                    float(s.get("dist_sum", 9999.0)),
                    str(s.get("site", "")),
                    int(s.get("index", 10**9)),
                ),
            )[0]
            non_his_value = best["non_his_count"]
            note_sites = sorted(
                qualified,
                key=lambda s: (
                    int(s.get("non_his_count", 10**9)),
                    -int(s.get("his_count", 0)),
                    str(s.get("site", "")),
                    int(s.get("index", 10**9)),
                ),
            )
        else:
            best = sorted(sites, key=site_sort_key)[0]
            non_his_value = ""
            note_sites = sorted(sites, key=site_sort_key)
    else:
        best = sorted(sites, key=site_sort_key)[0]
        non_his_value = best["non_his_count"]
        note_sites = sites
    others = [s for s in note_sites if s is not best]
    return {
        COORD_COLUMN_RESIDUE: best["residue_count"],
        COORD_COLUMN_HIS: best["his_count"],
        COORD_COLUMN_NON_HIS: non_his_value,
        COORD_COLUMN_RESIDUES: best["residues_text"],
        COORD_COLUMN_SITE: best["site"],
        COORD_COLUMN_NOTE: " | ".join(summarize_site(site) for site in others),
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


def parse_pose_meta(path: Path):
    name = path.name
    rank_match = re.search(r"rank(\d+)", name, flags=re.IGNORECASE)
    conf_match = re.search(r"confidence-([-+]?[0-9]*\.?[0-9]+)", name, flags=re.IGNORECASE)
    rank = int(rank_match.group(1)) if rank_match else None
    confidence = float(conf_match.group(1)) if conf_match else None
    return {"rank": rank, "confidence": confidence}


def better_pose(a, b):
    if a is None:
        return True
    ac = a.get("confidence")
    bc = b.get("confidence")
    if bc is not None and ac is None:
        return True
    if bc is None and ac is not None:
        return False
    if bc is not None and ac is not None and bc != ac:
        return bc > ac
    ar = a.get("rank")
    br = b.get("rank")
    if br is not None and ar is None:
        return True
    if br is None and ar is not None:
        return False
    if br is not None and ar is not None and br != ar:
        return br < ar
    return str(b.get("filename") or "") < str(a.get("filename") or "")


def build_diffdock_index(valid_ids):
    diffdock_root = DATA_DIR / "diffdock"
    items = []
    if not diffdock_root.exists():
        return {"items": [], "total_proteins": 0}
    valid = {str(x or "").strip().upper() for x in (valid_ids or []) if str(x or "").strip()}
    for d in sorted([x for x in diffdock_root.iterdir() if x.is_dir()], key=lambda p: p.name.upper()):
        pid = d.name.strip().upper()
        if valid and pid not in valid:
            continue
        structure_exists = any((DATA_DIR / "structures" / f"{pid}{ext}").exists() for ext in [".pdb", ".ent", ".cif"])
        if not structure_exists:
            continue
        grouped = {}
        fallback = []
        for sdf in sorted(d.glob("*.sdf"), key=lambda p: p.name):
            meta = parse_pose_meta(sdf)
            rec = {
                "file": f"../backend/data/diffdock/{pid}/{sdf.name}",
                "filename": sdf.name,
                "rank": meta["rank"],
                "confidence": meta["confidence"],
            }
            if meta["rank"] is None:
                continue
            if meta["confidence"] is None:
                fallback.append(rec)
                continue
            old = grouped.get(meta["rank"])
            if better_pose(old, rec):
                grouped[meta["rank"]] = rec
        for rec in fallback:
            rank = rec.get("rank")
            if rank is None or rank in grouped:
                continue
            grouped[rank] = rec
        poses = [grouped[k] for k in sorted(grouped.keys())]
        if not poses:
            continue
        items.append({
            "id": pid,
            "pose_count": len(poses),
            "poses": poses,
        })
    return {"items": items, "total_proteins": len(items)}

def main():
    master_path, rows = load_master()
    raw_header = rows[0]
    header = [name for name in raw_header if keep_column(name)]
    for extra in [COORD_COLUMN_RESIDUE, COORD_COLUMN_HIS, COORD_COLUMN_NON_HIS, COORD_COLUMN_RESIDUES, COORD_COLUMN_SITE, COORD_COLUMN_NOTE]:
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
    diffdock_payload = build_diffdock_index([item.get("id", "") for item in out])
    write_json(DATA_DIR / "diffdock_index.json", diffdock_payload)
    print("OK", DATA_DIR / "diffdock_index.json", diffdock_payload.get("total_proteins", 0))
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
