#!/usr/bin/env python3
"""
Approximate curved solvent channels from Zn to protein exterior on a 3D grid.

Outputs two channel-like quantities:
1. shortest_path_probe_radius_A
   Along the geometrically shortest escape path from Zn neighborhood to the
   local exterior, what is the narrowest clearance radius?
2. best_path_probe_radius_A
   Among all escape paths to the local exterior, what is the largest possible
   bottleneck clearance radius?

This is still an approximation, but it is materially closer to a true curved
channel search than fpocket alpha-sphere minima or straight-line ray tests.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
import sys
from pathlib import Path

import numpy as np


PROJECT_DIR = Path(__file__).resolve().parents[1]
BACKEND_DIR = PROJECT_DIR / "backend"
DATA_DIR = BACKEND_DIR / "data"
sys.path.insert(0, str(BACKEND_DIR))

from diffdock_api_server import _pdbzn_structure_atom_rows  # noqa: E402


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
ZN_RADIUS = 1.39
WATER_PROBE = 1.40
LOCAL_MARGIN = 16.0
BOUNDARY_PADDING = 5.0


def as_text(v) -> str:
    return "" if v is None else str(v)


def parse_ids_arg(raw: str):
    return [x.strip().upper() for x in as_text(raw).split(",") if x.strip()]


def dist(a, b) -> float:
    dx = float(a[0]) - float(b[0])
    dy = float(a[1]) - float(b[1])
    dz = float(a[2]) - float(b[2])
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def load_rows(table_id: str):
    path = DATA_DIR / f"{table_id}_master_full.json"
    obj = json.loads(path.read_text(encoding="utf-8"))
    return {as_text(r.get("Representative") or r.get("_id")).strip().upper(): r for r in obj.get("rows", [])}


def structure_path_for(pid: str):
    for cand in [
        DATA_DIR / "structures" / f"{pid}.pdb",
        DATA_DIR / "structures" / f"{pid}.cif",
        DATA_DIR / "structures" / f"{pid.lower()}.pdb",
        DATA_DIR / "structures" / f"{pid.lower()}.cif",
    ]:
        if cand.exists():
            return cand
    return None


def find_target_zn(atoms, row):
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
        return atom
    for atom in atoms:
        if as_text(atom.get("element")).upper() == "ZN":
            return atom
    return None


def protein_atoms_same_chain(atoms, chain_id: str):
    out = []
    for atom in atoms:
        if chain_id and as_text(atom.get("chain")).strip() != chain_id:
            continue
        comp = as_text(atom.get("comp")).strip().upper()
        if comp in WATER_RESIDUES:
            continue
        element = as_text(atom.get("element")).strip().upper()
        if element == "H":
            continue
        out.append(atom)
    return out


def build_local_atom_cloud(atom_rows, zn_atom):
    zn = np.array([float(zn_atom["x"]), float(zn_atom["y"]), float(zn_atom["z"])], dtype=np.float32)
    coords = []
    radii = []
    for atom in atom_rows:
        if atom is zn_atom:
            continue
        center = np.array([float(atom["x"]), float(atom["y"]), float(atom["z"])], dtype=np.float32)
        if np.max(np.abs(center - zn)) > LOCAL_MARGIN:
            continue
        coords.append(center)
        radii.append(float(VDW.get(as_text(atom.get("element")).upper(), 1.70)))
    if not coords:
        return np.zeros((0, 3), dtype=np.float32), np.zeros((0,), dtype=np.float32), zn
    return np.vstack(coords), np.array(radii, dtype=np.float32), zn


def make_clearance_grid(coords, radii, zn, spacing: float):
    if coords.shape[0] == 0:
        mins = zn - (LOCAL_MARGIN + BOUNDARY_PADDING)
        maxs = zn + (LOCAL_MARGIN + BOUNDARY_PADDING)
    else:
        mins = np.min(coords - radii[:, None], axis=0)
        maxs = np.max(coords + radii[:, None], axis=0)
        mins = np.minimum(mins, zn - LOCAL_MARGIN)
        maxs = np.maximum(maxs, zn + LOCAL_MARGIN)
        mins = mins - BOUNDARY_PADDING
        maxs = maxs + BOUNDARY_PADDING

    shape = np.ceil((maxs - mins) / spacing).astype(int) + 1
    xs = mins[0] + np.arange(shape[0], dtype=np.float32) * spacing
    ys = mins[1] + np.arange(shape[1], dtype=np.float32) * spacing
    zs = mins[2] + np.arange(shape[2], dtype=np.float32) * spacing
    clearance = np.full((shape[0], shape[1], shape[2]), np.inf, dtype=np.float32)

    for ix0 in range(0, shape[0], 8):
        ix1 = min(shape[0], ix0 + 8)
        X, Y, Z = np.meshgrid(xs[ix0:ix1], ys, zs, indexing="ij")
        chunk = np.full(X.shape, np.inf, dtype=np.float32)
        for center, radius in zip(coords, radii):
            dx = X - center[0]
            dy = Y - center[1]
            dz = Z - center[2]
            d = np.sqrt(dx * dx + dy * dy + dz * dz) - radius
            chunk = np.minimum(chunk, d.astype(np.float32))
        clearance[ix0:ix1, :, :] = chunk

    start_idx = np.round((zn - mins) / spacing).astype(int)
    start_idx = np.clip(start_idx, 0, shape - 1)
    return clearance, tuple(start_idx.tolist())


def neighbor_offsets():
    out = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if dx == dy == dz == 0:
                    continue
                out.append((dx, dy, dz, math.sqrt(dx * dx + dy * dy + dz * dz)))
    return out


NEIGHBORS = neighbor_offsets()


def is_boundary(idx, shape):
    x, y, z = idx
    sx, sy, sz = shape
    return x == 0 or y == 0 or z == 0 or x == sx - 1 or y == sy - 1 or z == sz - 1


def raw_accessible(clearance_value: float) -> float:
    return float(clearance_value)


def seed_accessible_voxels(clearance, start_idx, max_shell=6):
    sx, sy, sz = start_idx
    shape = clearance.shape
    seeds = []
    for shell in range(0, max_shell + 1):
        for i in range(max(0, sx - shell), min(shape[0], sx + shell + 1)):
            for j in range(max(0, sy - shell), min(shape[1], sy + shell + 1)):
                for k in range(max(0, sz - shell), min(shape[2], sz + shell + 1)):
                    if max(abs(i - sx), abs(j - sy), abs(k - sz)) != shell:
                        continue
                    local_probe = max(0.0, raw_accessible(clearance[i, j, k]))
                    if local_probe > 0:
                        seeds.append(((i, j, k), local_probe))
        if seeds:
            return seeds
    return seeds


def shortest_surface_distance(clearance, start_idx):
    shape = clearance.shape
    seeds = seed_accessible_voxels(clearance, start_idx)
    if not seeds:
        return None
    pq = []
    seen = {}
    for idx, _probe in seeds:
        seed_cost = dist(idx, start_idx)
        pq.append((seed_cost, idx))
        seen[idx] = seed_cost
    heapq.heapify(pq)
    while pq:
        path_len, idx = heapq.heappop(pq)
        if path_len > seen.get(idx, float("inf")) + 1e-9:
            continue
        if is_boundary(idx, shape):
            return path_len
        x, y, z = idx
        for dx, dy, dz, step_len in NEIGHBORS:
            nx, ny, nz = x + dx, y + dy, z + dz
            if nx < 0 or ny < 0 or nz < 0 or nx >= shape[0] or ny >= shape[1] or nz >= shape[2]:
                continue
            if raw_accessible(clearance[nx, ny, nz]) <= 0:
                continue
            nidx = (nx, ny, nz)
            nd = path_len + step_len
            if nd + 1e-9 < seen.get(nidx, float("inf")):
                seen[nidx] = nd
                heapq.heappush(pq, (nd, nidx))
    return None


def shortest_path_probe_radius(clearance, start_idx):
    shape = clearance.shape
    seeds = seed_accessible_voxels(clearance, start_idx)
    if not seeds:
        return None
    pq = []
    best_len = {}
    best_bottle = {}
    for idx, start_probe in seeds:
        seed_cost = dist(idx, start_idx)
        pq.append((seed_cost, idx, start_probe))
        best_len[idx] = seed_cost
        best_bottle[idx] = start_probe
    heapq.heapify(pq)
    while pq:
        path_len, idx, bottle = heapq.heappop(pq)
        if path_len > best_len.get(idx, float("inf")) + 1e-9:
            continue
        if is_boundary(idx, shape):
            return bottle
        x, y, z = idx
        for dx, dy, dz, step_len in NEIGHBORS:
            nx, ny, nz = x + dx, y + dy, z + dz
            if nx < 0 or ny < 0 or nz < 0 or nx >= shape[0] or ny >= shape[1] or nz >= shape[2]:
                continue
            local_probe = max(0.0, raw_accessible(clearance[nx, ny, nz]))
            if local_probe <= 0:
                continue
            nidx = (nx, ny, nz)
            nd = path_len + step_len
            nb = min(bottle, local_probe)
            if nd + 1e-9 < best_len.get(nidx, float("inf")):
                best_len[nidx] = nd
                best_bottle[nidx] = nb
                heapq.heappush(pq, (nd, nidx, nb))
    return None


def best_path_probe_radius(clearance, start_idx):
    shape = clearance.shape
    seeds = seed_accessible_voxels(clearance, start_idx)
    if not seeds:
        return None
    best = {}
    pq = []
    for idx, start_probe in seeds:
        best[idx] = start_probe
        pq.append((-start_probe, idx))
    heapq.heapify(pq)
    while pq:
        neg_score, idx = heapq.heappop(pq)
        score = -neg_score
        if score + 1e-9 < best.get(idx, -1.0):
            continue
        if is_boundary(idx, shape):
            return score
        x, y, z = idx
        for dx, dy, dz, _step_len in NEIGHBORS:
            nx, ny, nz = x + dx, y + dy, z + dz
            if nx < 0 or ny < 0 or nz < 0 or nx >= shape[0] or ny >= shape[1] or nz >= shape[2]:
                continue
            local_probe = max(0.0, raw_accessible(clearance[nx, ny, nz]))
            if local_probe <= 0:
                continue
            nidx = (nx, ny, nz)
            cand = min(score, local_probe)
            if cand > best.get(nidx, -1.0) + 1e-9:
                best[nidx] = cand
                heapq.heappush(pq, (-cand, nidx))
    return None


def analyze_one(pid: str, row, spacing: float):
    spath = structure_path_for(pid)
    if spath is None:
        return {"pid": pid, "error": "structure missing"}
    atoms = _pdbzn_structure_atom_rows(spath)
    zn_atom = find_target_zn(atoms, row)
    if zn_atom is None:
        return {"pid": pid, "error": "zn missing"}
    chain_id = as_text(zn_atom.get("chain")).strip()
    scoped_atoms = protein_atoms_same_chain(atoms, chain_id)
    coords, radii, zn = build_local_atom_cloud(scoped_atoms, zn_atom)
    clearance, start_idx = make_clearance_grid(coords, radii, zn, spacing)
    return {
        "pid": pid,
        "spacing_A": spacing,
        "shortest_surface_path_length_A": None if shortest_surface_distance(clearance, start_idx) is None else round(shortest_surface_distance(clearance, start_idx) * spacing, 3),
        "shortest_path_probe_radius_A": None if shortest_path_probe_radius(clearance, start_idx) is None else round(shortest_path_probe_radius(clearance, start_idx), 3),
        "best_path_probe_radius_A": None if best_path_probe_radius(clearance, start_idx) is None else round(best_path_probe_radius(clearance, start_idx), 3),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table-id", default="workflow_table")
    parser.add_argument("--ids", required=True, help="Comma-separated PDB ids")
    parser.add_argument("--spacing", type=float, default=1.0)
    parser.add_argument("--output-json", default="", help="Optional JSON file to write results")
    args = parser.parse_args()

    rows = load_rows(args.table_id)
    results = []
    for pid in parse_ids_arg(args.ids):
        row = rows.get(pid)
        if row is None:
            item = {"pid": pid, "error": "row missing"}
            print(json.dumps(item, ensure_ascii=False))
            results.append(item)
            continue
        item = analyze_one(pid, row, args.spacing)
        print(json.dumps(item, ensure_ascii=False))
        results.append(item)
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(results, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
