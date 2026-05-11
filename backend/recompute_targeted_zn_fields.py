#!/usr/bin/env python3

import csv
import json
import math
import urllib.error
import urllib.request
from pathlib import Path

from diffdock_api_server import (
    DATA_DIR,
    _pdbzn_residue_key,
    _pdbzn_residue_label,
    _pdbzn_structure_atom_rows,
)

HIS_RESIDUES = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}
WATER_RESIDUES = {"HOH", "WAT", "DOD"}
PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
    "GLY", "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS",
    "LYN", "MET", "MSE", "PHE", "PRO", "SEC", "SER", "THR", "TRP", "TYR", "VAL",
}
PROTEIN_DONOR_ATOMS = {
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
    "HIS": {"ND1", "NE2"},
    "HID": {"ND1", "NE2"},
    "HIE": {"ND1", "NE2"},
    "HIP": {"ND1", "NE2"},
    "HSD": {"ND1", "NE2"},
    "HSE": {"ND1", "NE2"},
    "HSP": {"ND1", "NE2"},
    "LYS": {"NZ"},
    "LYN": {"NZ"},
    "MET": {"SD"},
    "MSE": {"SE"},
    "SEC": {"SE"},
    "SER": {"OG"},
    "THR": {"OG1"},
    "TRP": {"NE1"},
    "TYR": {"OH"},
}
VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "SE": 1.90,
    "ZN": 1.39,
}
PROBE_RADIUS = 1.40
ZN_RADIUS = 1.39
SITE_DONOR_CUTOFF = 3.2
OCCLUSION_CUTOFF = 4.0
RAY_SAMPLES = 960
DOWNLOAD_TIMEOUT = 45
DOWNLOAD_BASE = "https://files.rcsb.org/download"
MASTER_JSONS = [
    DATA_DIR / "table1_master_full.json",
    DATA_DIR / "master_full.json",
]
MASTER_CSVS = [
    DATA_DIR / "table1_master_table.csv",
    DATA_DIR / "master_table.csv",
]
ITEM_JSONS = [
    DATA_DIR / "table1_data.json",
    DATA_DIR / "data.json",
]
JSON_FRONT_COLUMNS = [
    "Representative",
    "Zn_CoordResidueCount",
    "Zn_CoordHisCount",
    "Zn_CoordNonHisCount",
    "Zn_CoordSite",
    "Zn_CoordResidues",
    "ZN_Depth",
    "ZN_Depth_Rounded",
    "ZN_SASA",
    "ZN_Surface_Distance",
    "ZN_5A_ResidueCount",
    "ZN_5A_Residues",
    "TriHis_SameChain",
    "Zn_Binding_MonomerResidueCount",
]
CSV_FRONT_COLUMNS = [
    "_id",
    "Representative",
    "Zn_CoordResidueCount",
    "Zn_CoordHisCount",
    "Zn_CoordNonHisCount",
    "Zn_CoordSite",
    "Zn_CoordResidues",
    "ZN_Depth",
    "ZN_Depth_Rounded",
    "ZN_SASA",
    "ZN_Surface_Distance",
    "ZN_5A_ResidueCount",
    "ZN_5A_Residues",
    "TriHis_SameChain",
    "Zn_Binding_MonomerResidueCount",
]


def safe_exists(path_like):
    try:
        path = Path(path_like)
        return path.exists() and path.is_file()
    except Exception:
        return False


def dist(a, b):
    dx = float(a["x"]) - float(b["x"])
    dy = float(a["y"]) - float(b["y"])
    dz = float(a["z"]) - float(b["z"])
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def donor_atom_allowed(atom):
    comp = str((atom or {}).get("comp", "") or "").strip().upper()
    atom_name = str((atom or {}).get("atom", "") or "").strip().upper()
    return atom_name in PROTEIN_DONOR_ATOMS.get(comp, set())


def parse_num(value):
    try:
        return float(str(value).strip())
    except Exception:
        return None


def parse_int(value):
    try:
        return int(round(float(str(value).strip())))
    except Exception:
        return None


def golden_spiral_points(n):
    out = []
    ga = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        z = 1.0 - (2 * i + 1) / float(n)
        r = math.sqrt(max(0.0, 1.0 - z * z))
        t = i * ga
        out.append((math.cos(t) * r, math.sin(t) * r, z))
    return out


SPHERE_POINTS = golden_spiral_points(RAY_SAMPLES)
ZN_SPHERE_AREA = 4.0 * math.pi * ((ZN_RADIUS + PROBE_RADIUS) ** 2)


def structure_candidates(rep):
    pid = str(rep or "").strip().upper()
    return [
        DATA_DIR / "structures" / f"{pid}.pdb",
        DATA_DIR / "structures" / f"{pid}.cif",
        DATA_DIR / "structures" / f"{pid}.cif.gz",
        DATA_DIR / "structures" / f"{pid.lower()}.pdb",
        DATA_DIR / "structures" / f"{pid.lower()}.cif",
        DATA_DIR / "structures" / f"{pid.lower()}.cif.gz",
    ]


def first_local_structure(rep):
    for cand in structure_candidates(rep):
        if safe_exists(cand):
            return cand
    return None


def download_structure(rep):
    pid = str(rep or "").strip().upper()
    if not pid:
        return None
    struct_dir = DATA_DIR / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdb", "cif"):
        url = f"{DOWNLOAD_BASE}/{pid}.{ext}"
        target = struct_dir / f"{pid}.{ext}"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        try:
            with urllib.request.urlopen(req, timeout=DOWNLOAD_TIMEOUT) as resp:
                data = resp.read()
        except urllib.error.HTTPError:
            continue
        except Exception:
            continue
        if not data or len(data) < 100:
            continue
        target.write_bytes(data)
        return target
    return None


def ensure_local_structure(rep):
    hit = first_local_structure(rep)
    if hit is not None:
        return hit
    return download_structure(rep)


def old_site_label(row):
    return str((row or {}).get("Zn_CoordSite", "") or "").strip()


def parse_length_oligomer_monomer(value):
    parts = [x.strip() for x in str(value or "").split(";")]
    nums = [parse_int(x) for x in parts[:3]]
    while len(nums) < 3:
        nums.append(None)
    return nums[0], nums[1], nums[2]


def local_best_sdf_path(rep, row):
    raw = str((row or {}).get("Best_SDF", "") or "").strip()
    fname = Path(raw).name
    if not fname:
        return None
    path = DATA_DIR / "diffdock" / rep / fname
    return path if path.exists() else None


def read_sdf_points(path):
    pts = []
    if path is None or not path.exists():
        return pts
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 4:
        return pts
    try:
        atom_count = int(lines[3][:3])
    except Exception:
        return pts
    for line in lines[4 : 4 + atom_count]:
        try:
            pts.append(
                {
                    "x": float(line[0:10]),
                    "y": float(line[10:20]),
                    "z": float(line[20:30]),
                }
            )
        except Exception:
            continue
    return pts


def tri_his_sites(atoms):
    zn_atoms = [
        atom
        for atom in atoms
        if str(atom.get("comp") or "").upper() == "ZN" or str(atom.get("element") or "").upper() == "ZN"
    ]
    donor_atoms = [atom for atom in atoms if donor_atom_allowed(atom)]
    sites = []
    max_his = 0
    for zn in zn_atoms:
        zn_chain = str(zn.get("chain") or "").strip()
        nearest = {}
        for atom in donor_atoms:
            if zn_chain and str(atom.get("chain") or "").strip() != zn_chain:
                continue
            d = dist(zn, atom)
            if d > SITE_DONOR_CUTOFF:
                continue
            rk = _pdbzn_residue_key(atom)
            old = nearest.get(rk)
            if old is None or d < old["distance"]:
                nearest[rk] = {"atom": atom, "distance": float(d)}
        ordered = sorted(nearest.items(), key=lambda kv: kv[1]["distance"])
        his_count = sum(1 for rk, _ in ordered if rk[3] in HIS_RESIDUES)
        max_his = max(max_his, his_count)
        if his_count < 3:
            continue
        sites.append(
            {
                "label": f"{str(zn.get('chain') or '?')}:{str(zn.get('seq') or '?')}",
                "zn": zn,
                "nearest": nearest,
                "ordered": ordered,
                "his_count": his_count,
            }
        )
    return sites, max_his


def pick_site_by_points(sites, points):
    if not sites:
        return None
    if not points:
        return sites[0]
    best = None
    for site in sites:
        dmin = min(dist(site["zn"], pt) for pt in points)
        item = (float(dmin), site["label"])
        if best is None or item < best[0]:
            best = (item, site)
    return best[1]


def mean_site_distance(site):
    if not site or not site.get("ordered"):
        return 999.0
    return sum(item["distance"] for _, item in site["ordered"]) / float(len(site["ordered"]))


def pick_site(rep, row, sites):
    if not sites:
        return None
    label = old_site_label(row)
    if label:
        for site in sites:
            if site["label"] == label:
                return site
    sdf_points = read_sdf_points(local_best_sdf_path(rep, row))
    if sdf_points:
        return pick_site_by_points(sites, sdf_points)
    return sorted(
        sites,
        key=lambda site: (
            -int(site.get("his_count") or 0),
            -len(site.get("ordered") or []),
            mean_site_distance(site),
            site.get("label") or "",
        ),
    )[0]


def is_protein_residue_name(comp):
    return str(comp or "").strip().upper() in PROTEIN_RESIDUES


def residue_label_from_parts(comp, chain, seq):
    return f"{str(comp or '').upper()}:{str(chain or '').strip()}:{str(seq or '').strip()}"


def site_chain(site):
    if site is None:
        return ""
    zn = site.get("zn") or {}
    return str(zn.get("chain") or "").strip()


def filter_atoms_to_site_chain(atoms, site):
    chain = site_chain(site)
    if not chain:
        return list(atoms or [])
    out = []
    for atom in atoms or []:
        atom_chain = str((atom or {}).get("chain") or "").strip()
        if atom_chain == chain:
            out.append(atom)
    return out


def build_occ_spheres(atoms, zn_atom):
    occupiers = []
    for atom in atoms:
        if atom is zn_atom:
            continue
        comp = str(atom.get("comp") or "").upper()
        if comp in WATER_RESIDUES:
            continue
        element = str(atom.get("element") or "").upper()
        if element == "H":
            continue
        center_dist = dist(zn_atom, atom)
        if center_dist > OCCLUSION_CUTOFF:
            continue
        occupiers.append(
            (
                float(atom["x"]) - float(zn_atom["x"]),
                float(atom["y"]) - float(zn_atom["y"]),
                float(atom["z"]) - float(zn_atom["z"]),
                VDW_RADII.get(element, 1.70) + PROBE_RADIUS,
            )
        )
    return occupiers


def build_occ_spheres_protein_only(atoms, zn_atom, site):
    occupiers = []
    allowed_keys = set()
    zn_chain = site_chain(site)
    if site is not None:
        for residue_key, _ in site.get("ordered") or []:
            allowed_keys.add(residue_key)
    for atom in atoms:
        if atom is zn_atom:
            continue
        atom_chain = str(atom.get("chain") or "").strip()
        if zn_chain and atom_chain != zn_chain:
            continue
        comp = str(atom.get("comp") or "").upper()
        if comp in WATER_RESIDUES:
            continue
        record = str(atom.get("record") or "").upper()
        rk = _pdbzn_residue_key(atom)
        if record == "HETATM" and rk not in allowed_keys:
            continue
        element = str(atom.get("element") or "").upper()
        if element == "H":
            continue
        center_dist = dist(zn_atom, atom)
        if center_dist > OCCLUSION_CUTOFF:
            continue
        occupiers.append(
            (
                float(atom["x"]) - float(zn_atom["x"]),
                float(atom["y"]) - float(zn_atom["y"]),
                float(atom["z"]) - float(zn_atom["z"]),
                VDW_RADII.get(element, 1.70) + PROBE_RADIUS,
            )
        )
    return occupiers


def build_occ_spheres_coordination_only(site):
    occupiers = []
    if site is None:
        return occupiers
    zn_atom = site["zn"]
    for _, item in site.get("ordered") or []:
        atom = item.get("atom")
        if not atom:
            continue
        element = str(atom.get("element") or "").upper()
        if element == "H":
            continue
        occupiers.append(
            (
                float(atom["x"]) - float(zn_atom["x"]),
                float(atom["y"]) - float(zn_atom["y"]),
                float(atom["z"]) - float(zn_atom["z"]),
                VDW_RADII.get(element, 1.70) + PROBE_RADIUS,
            )
        )
    return occupiers


def ray_escape_distance(occupiers, direction):
    ux, uy, uz = direction
    intervals = []
    for cx, cy, cz, radius in occupiers:
        proj = cx * ux + cy * uy + cz * uz
        perp2 = (cx * cx + cy * cy + cz * cz) - proj * proj
        r2 = radius * radius
        if perp2 >= r2:
            continue
        half = math.sqrt(max(0.0, r2 - perp2))
        start = proj - half
        end = proj + half
        if end <= 0.0:
            continue
        if start < 0.0:
            start = 0.0
        intervals.append((start, end))
    if not intervals:
        return 0.0
    intervals.sort()
    if intervals[0][0] > 1e-6:
        return 0.0
    current_end = intervals[0][1]
    for start, end in intervals[1:]:
        if start > current_end + 1e-6:
            break
        if end > current_end:
            current_end = end
    return current_end


def zn_occlusion_metrics(atoms, site):
    if site is None:
        return {"ZN_SASA": "", "ZN_Depth": "", "ZN_Depth_Rounded": "", "ZN_Surface_Distance": ""}
    shell_r = ZN_RADIUS + PROBE_RADIUS
    monomer_atoms = filter_atoms_to_site_chain(atoms, site)
    depth_occupiers = build_occ_spheres_protein_only(monomer_atoms, site["zn"], site)
    sasa_occupiers = build_occ_spheres_coordination_only(site)
    if not depth_occupiers:
        return {
            "ZN_SASA": f"{ZN_SPHERE_AREA:.3f}",
            "ZN_Depth": "0.000",
            "ZN_Depth_Rounded": "0.0",
            "ZN_Surface_Distance": "0.000",
        }
    exposed = 0
    directional_depths = []
    for direction in SPHERE_POINTS:
        escape = ray_escape_distance(depth_occupiers, direction)
        directional_depths.append(max(0.0, escape - ZN_RADIUS))
        if not sasa_occupiers:
            exposed += 1
            continue
        coord_escape = ray_escape_distance(sasa_occupiers, direction)
        if coord_escape < shell_r - 1e-6:
            exposed += 1
    sasa = ZN_SPHERE_AREA * exposed / float(len(SPHERE_POINTS))
    mean_depth = sum(directional_depths) / float(len(directional_depths))
    surface_distance = min(directional_depths) if directional_depths else 0.0
    return {
        "ZN_SASA": f"{sasa:.3f}",
        "ZN_Depth": f"{mean_depth:.3f}",
        "ZN_Depth_Rounded": f"{mean_depth:.1f}",
        "ZN_Surface_Distance": f"{surface_distance:.3f}",
    }


def coord_metrics(site, max_his):
    if site is None:
        return {
            "TriHisSatisfied": "false",
            "TriHisCountMax": str(int(max_his or 0)),
            "Zn_CoordResidueCount": "0",
            "Zn_CoordHisCount": "0",
            "Zn_CoordNonHisCount": "0",
            "Zn_CoordResidues": "",
            "Zn_CoordSite": "",
            "TriHis_SameChain": "false",
        }
    ordered = list(site["ordered"])
    residue_count = len(ordered)
    his_count = sum(1 for rk, _ in ordered if rk[3] in HIS_RESIDUES)
    his_chains = sorted({str(rk[0] or "").strip() for rk, _ in ordered if rk[3] in HIS_RESIDUES})
    residue_labels = [f"{_pdbzn_residue_label(rk)}:{item['distance']:.3f}" for rk, item in ordered]
    return {
        "TriHisSatisfied": "true" if his_count >= 3 else "false",
        "TriHisCountMax": str(int(max_his or his_count)),
        "Zn_CoordResidueCount": str(residue_count),
        "Zn_CoordHisCount": str(his_count),
        "Zn_CoordNonHisCount": str(max(0, residue_count - his_count)),
        "Zn_CoordResidues": "; ".join(residue_labels),
        "Zn_CoordSite": site["label"],
        "TriHis_SameChain": "true" if his_count >= 3 and len(his_chains) == 1 else "false",
    }


def zn_5a_metrics(atoms, site):
    if site is None:
        return {"ZN_5A_ResidueCount": "0", "ZN_5A_Residues": ""}
    zn = site["zn"]
    chain = site_chain(site)
    nearest = {}
    for atom in atoms:
        if chain and str(atom.get("chain") or "").strip() != chain:
            continue
        comp = str(atom.get("comp") or "").upper()
        if not is_protein_residue_name(comp):
            continue
        key = _pdbzn_residue_key(atom)
        d = dist(zn, atom)
        old = nearest.get(key)
        if old is None or d < old["distance"]:
            nearest[key] = {"distance": float(d), "atom": atom}
    ordered = sorted(
        [(rk, item) for rk, item in nearest.items() if item["distance"] <= 5.0],
        key=lambda kv: (kv[1]["distance"], _pdbzn_residue_label(kv[0])),
    )
    labels = [f"{_pdbzn_residue_label(rk)}:{item['distance']:.3f}" for rk, item in ordered]
    return {
        "ZN_5A_ResidueCount": str(len(ordered)),
        "ZN_5A_Residues": "; ".join(labels),
    }


def zn_binding_monomer_residue_count(row, atoms, site):
    _, oligomer, monomer = parse_length_oligomer_monomer((row or {}).get("Length;Oligomer;Monomer"))
    if monomer is not None and monomer > 0:
        return str(monomer)
    total_len, oligomer, _ = parse_length_oligomer_monomer((row or {}).get("Length;Oligomer;Monomer"))
    if total_len is not None and oligomer is not None and oligomer > 0:
        return str(int(round(float(total_len) / float(oligomer))))
    if site is None:
        return ""
    zn_chain = str(site["zn"].get("chain") or "").strip()
    residues = set()
    for atom in atoms:
        if str(atom.get("chain") or "").strip() != zn_chain:
            continue
        comp = str(atom.get("comp") or "").upper()
        if not is_protein_residue_name(comp):
            continue
        residues.add(_pdbzn_residue_key(atom))
    return str(len(residues)) if residues else ""


def recompute_one(rep, row):
    structure = ensure_local_structure(rep)
    if structure is None:
        return {}
    atoms = _pdbzn_structure_atom_rows(structure)
    sites, max_his = tri_his_sites(atoms)
    site = pick_site(rep, row, sites)
    updates = {}
    updates.update(coord_metrics(site, max_his))
    updates.update(zn_occlusion_metrics(atoms, site))
    updates.update(zn_5a_metrics(atoms, site))
    updates["Zn_Binding_MonomerResidueCount"] = zn_binding_monomer_residue_count(row, atoms, site)
    return updates


def reorder_columns(existing, preferred):
    existing_list = [str(x or "").strip() for x in existing if str(x or "").strip()]
    seen = set()
    out = []
    for name in list(preferred) + existing_list:
        key = str(name or "").strip()
        if not key or key in seen or key not in existing_list:
            continue
        seen.add(key)
        out.append(key)
    return out


def reordered_row(row, ordered_keys):
    out = {}
    for key in ordered_keys:
        if key in row:
            out[key] = row.get(key, "")
    for key, value in row.items():
        if key not in out:
            out[key] = value
    return out


def update_csv(path, updates_by_rep):
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        rows = list(csv.DictReader(f))
        fieldnames = list(rows[0].keys()) if rows else []
    changed = 0
    for row in rows:
        rep = str((row.get("Representative") or row.get("_id") or "")).strip().upper()
        updates = updates_by_rep.get(rep)
        if not updates:
            continue
        for key, value in updates.items():
            row[key] = value
        changed += 1
    for updates in updates_by_rep.values():
        for key in updates.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    fieldnames = reorder_columns(fieldnames, CSV_FRONT_COLUMNS)
    rows = [reordered_row(row, fieldnames) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    return changed


def update_master_full(path, updates_by_rep):
    obj = json.loads(path.read_text(encoding="utf-8"))
    rows = obj.get("rows", [])
    changed = 0
    for row in rows:
        rep = str((row.get("Representative") or row.get("_id") or "")).strip().upper()
        updates = updates_by_rep.get(rep)
        if not updates:
            continue
        for key, value in updates.items():
            row[key] = value
        changed += 1
    header = list(obj.get("header") or list(rows[0].keys() if rows else []))
    for updates in updates_by_rep.values():
        for key in updates.keys():
            if key not in header:
                header.append(key)
    header = reorder_columns(header, JSON_FRONT_COLUMNS)
    obj["header"] = header
    obj["rows"] = [reordered_row(row, header) for row in rows]
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")
    return changed


def update_data_json(path, updates_by_rep):
    obj = json.loads(path.read_text(encoding="utf-8"))
    changed = 0
    for item in obj.get("items", []):
        rep = str((item.get("id") or "")).strip().upper()
        updates = updates_by_rep.get(rep)
        if not updates:
            continue
        for src, dst in [
            ("Zn_CoordResidueCount", "zn_coord_residue_count"),
            ("Zn_CoordHisCount", "zn_coord_his_count"),
            ("Zn_CoordNonHisCount", "zn_coord_non_his_count"),
        ]:
            if src in updates and dst in item:
                try:
                    item[dst] = int(float(updates[src]))
                except Exception:
                    item[dst] = updates[src]
        text_map = [
            ("ZN_Surface_Distance", "zn_surface_distance"),
            ("ZN_5A_ResidueCount", "zn_5a_residue_count"),
            ("ZN_5A_Residues", "zn_5a_residues"),
            ("TriHis_SameChain", "tri_his_same_chain"),
            ("Zn_Binding_MonomerResidueCount", "zn_binding_monomer_residue_count"),
        ]
        for src, dst in text_map:
            if src in updates:
                item[dst] = updates[src]
        changed += 1
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")
    return changed


def load_primary_rows():
    with (DATA_DIR / "table1_master_table.csv").open(newline="", encoding="utf-8", errors="ignore") as f:
        return list(csv.DictReader(f))


def main():
    rows = load_primary_rows()
    reps = []
    by_rep = {}
    for row in rows:
        rep = str((row.get("Representative") or row.get("_id") or "")).strip().upper()
        if not rep or rep in by_rep:
            continue
        reps.append(rep)
        by_rep[rep] = row

    updates_by_rep = {}
    structure_hits = 0
    for rep in reps:
        updates = recompute_one(rep, by_rep[rep])
        if updates:
            structure_hits += 1
            updates_by_rep[rep] = updates

    for path in MASTER_CSVS:
        if path.exists():
            print("updated csv", path, update_csv(path, updates_by_rep))
    for path in MASTER_JSONS:
        if path.exists():
            print("updated master_full", path, update_master_full(path, updates_by_rep))
    for path in ITEM_JSONS:
        if path.exists():
            print("updated data", path, update_data_json(path, updates_by_rep))

    print("structure_hits", structure_hits, "of", len(reps))
    for rep in reps[:10]:
        print(rep, json.dumps(updates_by_rep.get(rep, {}), ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
