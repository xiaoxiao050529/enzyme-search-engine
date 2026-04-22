import json
import gzip
import math
import mimetypes
import os
import re
import shlex
import shutil
import sqlite3
import subprocess
import threading
import time
import traceback
import uuid
import csv
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse


BASE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = BASE_DIR.parent
DATA_DIR = BASE_DIR / "data"
RUNTIME_DIR = BASE_DIR / "runtime"
TABLE_REGISTRY_PATH = DATA_DIR / "table_registry.json"
IMPORT_SELECTION_PATH = RUNTIME_DIR / "import_selection.json"
IMPORT_STRUCT_DIR = RUNTIME_DIR / "import_structures"

_pdbzn_tmalign_cache = {}
_pdbzn_tmalign_cache_lock = threading.Lock()
_pdbzn_table_registry_lock = threading.Lock()
_pdbzn_pocket_table_cache = None
_pdbzn_pocket_table_cache_lock = threading.Lock()

IMPORT_STRUCT_DIR.mkdir(parents=True, exist_ok=True)


def _env_dir(name):
    raw = str(os.environ.get(name, "") or "").strip()
    if not raw:
        return None
    p = Path(raw).expanduser()
    if p.exists() and p.is_dir():
        return p.resolve()
    return None


def _first_existing_dir(candidates):
    for p in candidates:
        try:
            if p and p.exists() and p.is_dir():
                return p.resolve()
        except Exception:
            continue
    return None


def _safe_path_exists(path_like):
    try:
        return Path(path_like).exists()
    except Exception:
        return False


def _safe_existing_file(path_like):
    if not path_like:
        return None
    try:
        p = Path(path_like)
        if p.exists() and p.is_file():
            return p
    except Exception:
        return None
    return None


def _safe_existing_dir(path_like):
    if not path_like:
        return None
    try:
        p = Path(path_like)
        if p.exists() and p.is_dir():
            return p
    except Exception:
        return None
    return None


def discover_pdbzn_base_dir():
    env = _env_dir("PDBZN_BASE_DIR")
    if env:
        return env
    found = _first_existing_dir(
        [
            PROJECT_DIR / "PDB_ZN",
            DATA_DIR,
            Path("/root/PDB_ZN"),
        ]
    )
    return found or DATA_DIR


def discover_pdbzn_pdb_dir(base_dir: Path):
    env = _env_dir("PDBZN_PDB_DIR")
    if env:
        return env
    found = _first_existing_dir(
        [
            base_dir / "pdb_structures",
            base_dir / "structures",
            DATA_DIR / "structures",
            base_dir,
        ]
    )
    return found or (DATA_DIR / "structures")


WORK_DIR = RUNTIME_DIR / "diffdock_jobs"
WORK_DIR.mkdir(parents=True, exist_ok=True)
FPOCKET_WORK_DIR = RUNTIME_DIR / "fpocket_jobs"
FPOCKET_WORK_DIR.mkdir(parents=True, exist_ok=True)
TMALIGN_WORK_DIR = RUNTIME_DIR / "tmalign_jobs"
TMALIGN_WORK_DIR.mkdir(parents=True, exist_ok=True)
PDBZN_BASE_DIR = discover_pdbzn_base_dir()
PDBZN_PDB_DIR = discover_pdbzn_pdb_dir(PDBZN_BASE_DIR)
PDBZN_DB_PATH = RUNTIME_DIR / "pdbzn.sqlite"
PDBZN_EXPORT_DIR = RUNTIME_DIR / "exports"
PDBZN_EXPORT_DIR.mkdir(parents=True, exist_ok=True)
HOST = os.environ.get("DIFFDOCK_API_HOST", "127.0.0.1")
PORT = int(os.environ.get("DIFFDOCK_API_PORT", "8015"))


def discover_infer_py():
    env = os.environ.get("DIFFDOCK_INFER_PY", "").strip()
    if env and _safe_path_exists(env):
        return env
    candidates = [
        str(PROJECT_DIR / ".tools" / "src" / "DiffDock" / "inference.py"),
        "/root/DiffDock/inference.py",
        "/root/diffdock/inference.py",
        "/opt/DiffDock/inference.py",
    ]
    for p in candidates:
        if _safe_path_exists(p):
            return p
    return ""


INFER_PY = discover_infer_py()
PYTHON_BIN = os.environ.get(
    "DIFFDOCK_PYTHON_BIN",
    str(PROJECT_DIR / ".tools" / "diffdock-venv-sys" / "bin" / "python") if _safe_path_exists(PROJECT_DIR / ".tools" / "diffdock-venv-sys" / "bin" / "python") else (shutil.which("python3") or "python3"),
)
FPOCKET_BIN = os.environ.get("FPOCKET_BIN", "fpocket")
TMALIGN_BIN = os.environ.get("TMALIGN_BIN", "TMalign")

jobs = {}
jobs_lock = threading.Lock()
fpocket_jobs = {}
fpocket_jobs_lock = threading.Lock()
tmalign_jobs = {}
tmalign_jobs_lock = threading.Lock()
pdbzn_db_lock = threading.Lock()
PDBZN_LARGE_LIGANDS = {
    "NAD", "NAI", "NDP", "FAD", "FMN",
    "ATP", "ADP", "AMP", "GTP", "GDP", "SAM", "SAH",
    "COA", "ACO",
    "HEM", "HEC", "HOM",
    "PLP", "TPP", "BTN", "THF", "CNC",
    "NAG", "NDG", "MAN", "GLC",
}
PDBZN_METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA",
    "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG", "AL",
    "GA", "IN", "TL", "CR", "MO", "W", "V", "TI", "Y", "ZR",
}
PDBZN_HIS_RESIDUES = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}
PDBZN_PROTEIN_DONOR_ATOMS = {
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


def now_ts():
    return int(time.time())


def safe_name(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name or "file")


def list_pose_files(out_dir: Path):
    files = sorted([p for p in out_dir.glob("**/*.sdf") if p.is_file()], key=lambda x: x.name)
    out = []
    for p in files:
        out.append({
            "name": p.name,
            "size": p.stat().st_size,
            "url": f"/api/diffdock/file/{p.parent.name}/{p.name}",
        })
    return out


def discover_fpocket_bin():
    env = os.environ.get("FPOCKET_BIN", "").strip()
    if env:
        if _safe_path_exists(env):
            return env
        found = shutil.which(env)
        if found:
            return found
    local_candidates = [
        PROJECT_DIR / ".tools" / "bin" / "fpocket",
        PROJECT_DIR / ".tools" / "src" / "fpocket-src" / "bin" / "fpocket",
    ]
    for p in local_candidates:
        if _safe_path_exists(p):
            return str(p)
    found_default = shutil.which("fpocket")
    return found_default or ""


def list_fpocket_outputs(job_id: str):
    job_dir = FPOCKET_WORK_DIR / job_id
    roots = []
    out_dir = job_dir / "output"
    if out_dir.exists():
        roots.append(out_dir)
    in_dir = job_dir / "input"
    if in_dir.exists():
        for d in sorted(in_dir.glob("*_out")):
            if d.is_dir():
                roots.append(d)
    if not roots:
        return []
    items = []
    seen = set()
    for root in roots:
        for p in sorted(root.rglob("*")):
            if not p.is_file():
                continue
            rel = p.relative_to(job_dir).as_posix()
            if rel in seen:
                continue
            seen.add(rel)
            items.append({
                "name": p.name,
                "path": rel,
                "size": p.stat().st_size,
                "url": f"/api/fpocket/file/{job_id}/{rel}",
            })
    items.sort(key=lambda x: x.get("path", ""))
    return items


def parse_fpocket_info_file(info_path: Path):
    if not info_path.exists():
        return {"columns": [], "rows": []}
    text = info_path.read_text(encoding="utf-8", errors="ignore")
    rows = []
    current = None
    col_order = []
    for raw in text.splitlines():
        line = raw.strip()
        m = re.match(r"^Pocket\s+(\d+)\s*:", line, flags=re.IGNORECASE)
        if m:
            if current:
                rows.append(current)
            current = {"pocket_id": int(m.group(1)), "params": {}}
            continue
        if not current:
            continue
        if not line:
            continue
        if ":" not in line:
            continue
        key, val = line.split(":", 1)
        key = key.strip()
        val = val.strip()
        parsed_val = val
        try:
            parsed_val = float(val)
            if parsed_val.is_integer():
                parsed_val = int(parsed_val)
        except Exception:
            pass
        current["params"][key] = parsed_val
        if key not in col_order:
            col_order.append(key)
    if current:
        rows.append(current)
    return {"columns": col_order, "rows": rows}


def parse_fpocket_vert_file(vert_path: Path):
    points = []
    if not vert_path.exists():
        return points
    for raw in vert_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.rstrip("\n")
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        x = y = z = None
        r = None
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except Exception:
            pass
        toks = line.split()
        if x is None or y is None or z is None:
            if len(toks) >= 8:
                try:
                    x = float(toks[5])
                    y = float(toks[6])
                    z = float(toks[7])
                except Exception:
                    continue
        if len(toks) >= 2:
            try:
                r = float(toks[-1])
            except Exception:
                r = None
        points.append({"x": x, "y": y, "z": z, "r": r})
    return points


def build_fpocket_pocket_payload(job_id: str, job: dict):
    job_dir = FPOCKET_WORK_DIR / job_id
    in_dir = job_dir / "input"
    pdb_name = str((job or {}).get("pdb_name", "") or "")
    if not pdb_name and in_dir.exists():
        pdb_files = sorted([p for p in in_dir.iterdir() if p.is_file() and p.suffix.lower() in {".pdb", ".ent"}], key=lambda p: p.name)
        if pdb_files:
            pdb_name = pdb_files[0].name
    stem = Path(pdb_name).stem if pdb_name else ""
    root_candidates = []
    if stem:
        d = in_dir / f"{stem}_out"
        if d.exists() and d.is_dir():
            root_candidates.append(d)
    if in_dir.exists():
        for d in sorted(in_dir.glob("*_out")):
            if d.is_dir() and d not in root_candidates:
                root_candidates.append(d)
    if not root_candidates:
        return {"ok": False, "error": "未找到fpocket输出目录", "pockets": [], "columns": [], "rows": [], "structure_pdb_text": ""}
    out_root = root_candidates[0]
    info_path = out_root / f"{out_root.name.replace('_out', '')}_info.txt"
    if not info_path.exists():
        info_files = sorted(out_root.glob("*_info.txt"))
        if info_files:
            info_path = info_files[0]
    info = parse_fpocket_info_file(info_path)
    rows = info.get("rows", [])
    by_id = {int(r.get("pocket_id", -1)): r for r in rows if isinstance(r.get("pocket_id"), int)}
    pockets_dir = out_root / "pockets"
    vert_files = []
    if pockets_dir.exists() and pockets_dir.is_dir():
        vert_files = sorted(pockets_dir.glob("pocket*_vert.pqr"), key=lambda p: int(re.findall(r"\d+", p.name)[0]) if re.findall(r"\d+", p.name) else 10**9)
    if not vert_files:
        one = out_root / f"{out_root.name.replace('_out', '')}_pockets.pqr"
        if one.exists():
            vert_files = [one]
    pockets = []
    for vf in vert_files:
        nums = re.findall(r"\d+", vf.name)
        pocket_id = int(nums[0]) if nums else 1
        pts = parse_fpocket_vert_file(vf)
        row = by_id.get(pocket_id, {"pocket_id": pocket_id, "params": {}})
        pockets.append({
            "pocket_id": pocket_id,
            "point_count": len(pts),
            "points": pts,
            "params": row.get("params", {}),
            "vert_file": vf.relative_to(job_dir).as_posix(),
        })
    rows_sorted = sorted(rows, key=lambda x: int(x.get("pocket_id", 10**9)))
    structure_path = in_dir / pdb_name if pdb_name else None
    structure_text = ""
    if structure_path and structure_path.exists():
        structure_text = structure_path.read_text(encoding="utf-8", errors="ignore")
    if not structure_text:
        struct_candidates = []
        if stem:
            struct_candidates.append(out_root / f"{stem}_out.pdb")
        struct_candidates.extend(sorted(out_root.glob("*_out.pdb"), key=lambda p: p.name))
        for cand in struct_candidates:
            if cand.exists() and cand.is_file():
                structure_text = cand.read_text(encoding="utf-8", errors="ignore")
                break
    return {
        "ok": True,
        "output_root": out_root.relative_to(job_dir).as_posix(),
        "columns": info.get("columns", []),
        "rows": rows_sorted,
        "pockets": pockets,
        "structure_pdb_text": structure_text,
    }


def discover_tmalign_bin():
    env = os.environ.get("TMALIGN_BIN", "").strip()
    if env:
        if _safe_path_exists(env):
            return env
        found = shutil.which(env)
        if found:
            return found
    candidates = [
        str(PROJECT_DIR / ".tools" / "bin" / "TMalign"),
        "/root/tools/TMalign",
        "/usr/local/bin/TMalign",
        "/usr/bin/TMalign",
    ]
    for p in candidates:
        if _safe_path_exists(p):
            return p
    found_default = shutil.which("TMalign")
    return found_default or ""


def parse_tmalign_log(text: str):
    log = text or ""
    aligned_len = None
    rmsd = None
    seq_id = None
    tm_scores = []
    for m in re.finditer(r"TM-score=\s*([0-9]*\.?[0-9]+)", log):
        try:
            tm_scores.append(float(m.group(1)))
        except Exception:
            pass
    m = re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9]*\.?[0-9]+),\s*Seq_ID=n_identical/n_aligned=\s*([0-9]*\.?[0-9]+)", log)
    if m:
        try:
            aligned_len = int(m.group(1))
        except Exception:
            aligned_len = None
        try:
            rmsd = float(m.group(2))
        except Exception:
            rmsd = None
        try:
            seq_id = float(m.group(3))
        except Exception:
            seq_id = None
    return {
        "tm_score_1": tm_scores[0] if len(tm_scores) > 0 else None,
        "tm_score_2": tm_scores[1] if len(tm_scores) > 1 else None,
        "tm_score_max": max(tm_scores) if tm_scores else None,
        "aligned_length": aligned_len,
        "rmsd": rmsd,
        "seq_id": seq_id,
    }


def _pdbzn_data_aliases(name):
    alias_map = {
        "zn_his_master_table7.csv": [DATA_DIR / "master_table.csv"],
        "zn_his_master_table.csv": [DATA_DIR / "master_table.csv"],
        "zn_his_similarity_ranking2.csv": [DATA_DIR / "master_table.csv"],
        "zn_his_similarity_ranking.csv": [DATA_DIR / "master_table.csv"],
    }
    return list(alias_map.get(name, []))


def _pdbzn_candidate_files(name):
    cands = []
    cands.extend(_pdbzn_data_aliases(name))
    cands.extend(
        [
            PDBZN_BASE_DIR / name,
            DATA_DIR / name,
            BASE_DIR / name,
            PROJECT_DIR / name,
        ]
    )
    out = []
    seen = set()
    for p in cands:
        key = str(p)
        if key in seen:
            continue
        seen.add(key)
        out.append(p)
    return out


def _pdbzn_parse_uploaded_ids(raw):
    ids = []
    seen = set()
    text = str(raw or "")
    if not text.strip():
        return ids
    for raw_line in text.splitlines():
        line = str(raw_line or "").strip()
        if not line or line.startswith("#"):
            continue
        for token in re.split(r"[\s,;]+", line):
            seg = str(token or "").strip().strip("\"'`")
            if not seg or seg.startswith("#"):
                continue
            name = Path(seg).name
            lower = name.lower()
            for suffix in [".cif.gz", ".pdb.gz", ".ent.gz", ".cif", ".pdb", ".ent"]:
                if lower.endswith(suffix):
                    name = name[: -len(suffix)]
                    break
            pid = _pdbzn_guess_pdb_id(name)
            if not pid or pid in seen:
                continue
            seen.add(pid)
            ids.append(pid)
    return ids


def _pdbzn_active_import_selection():
    hit = _safe_existing_file(IMPORT_SELECTION_PATH)
    if hit is None:
        return None
    obj = _pdbzn_read_json(hit)
    if not isinstance(obj, dict):
        return None
    ids = []
    seen = set()
    for item in obj.get("ids", []) or []:
        pid = _pdbzn_guess_pdb_id(str(item or ""))
        if not pid or pid in seen:
            continue
        seen.add(pid)
        ids.append(pid)
    if not ids:
        return None
    obj["ids"] = ids
    return obj


def _pdbzn_write_import_selection(ids, source_name=""):
    payload = {
        "mode": "upload_list",
        "source_name": str(source_name or "").strip(),
        "ids": list(ids or []),
        "count": len(list(ids or [])),
        "updated_at": now_ts(),
        "updated_at_text": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
    }
    _pdbzn_write_json(IMPORT_SELECTION_PATH, payload)
    return payload


def _pdbzn_clear_import_selection():
    try:
        IMPORT_SELECTION_PATH.unlink()
    except FileNotFoundError:
        pass
    except Exception:
        pass


def _pdbzn_structure_roots():
    roots = [IMPORT_STRUCT_DIR, PDBZN_PDB_DIR, DATA_DIR / "structures", PDBZN_BASE_DIR]
    out = []
    seen = set()
    for root in roots:
        if not root:
            continue
        try:
            path = Path(root)
        except Exception:
            continue
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        if path.exists() and path.is_dir():
            out.append(path)
    return out


def _pdbzn_local_structure_path(pdb_id, prefer_pdb=False):
    pid = str(pdb_id or "").strip().upper()
    if not pid:
        return None
    suffixes = [".pdb", ".ent", ".cif.gz", ".cif"] if prefer_pdb else [".cif.gz", ".cif", ".pdb", ".ent"]
    candidates = []
    for root in _pdbzn_structure_roots():
        for stem in [pid, pid.lower()]:
            for suffix in suffixes:
                candidates.append(root / f"{stem}{suffix}")
    seen = set()
    for cand in candidates:
        key = str(cand)
        if key in seen:
            continue
        seen.add(key)
        hit = _safe_existing_file(cand)
        if hit is not None:
            return hit
    return None


def _pdbzn_download_structure_file(pdb_id):
    pid = str(pdb_id or "").strip().upper()
    if not pid:
        return None
    IMPORT_STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    attempts = [
        (f"{pid}.pdb", f"https://files.rcsb.org/download/{pid}.pdb"),
        (f"{pid}.cif.gz", f"https://files.rcsb.org/download/{pid}.cif.gz"),
        (f"{pid}.cif", f"https://files.rcsb.org/download/{pid}.cif"),
    ]
    for fname, url in attempts:
        out_path = IMPORT_STRUCT_DIR / fname
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Codex-PDBZN/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
            if not data:
                continue
            out_path.write_bytes(data)
            return out_path
        except Exception:
            try:
                out_path.unlink()
            except Exception:
                pass
            continue
    return None


def _pdbzn_resolve_import_file_map(ids, allow_download=False):
    file_map = {}
    missing = []
    for pid in ids or []:
        hit = _pdbzn_local_structure_path(pid, prefer_pdb=True)
        if hit is None and allow_download:
            hit = _pdbzn_download_structure_file(pid)
        if hit is None:
            missing.append(pid)
            continue
        file_map[pid] = hit
    return file_map, missing


def _pdbzn_source_structure_files():
    selection = _pdbzn_active_import_selection()
    if selection and selection.get("ids"):
        file_map, _ = _pdbzn_resolve_import_file_map(selection.get("ids") or [], allow_download=False)
        out = list(file_map.values())
        out.sort(key=lambda p: p.name.lower())
        return out
    source_dirs = []
    for d in _pdbzn_structure_roots():
        if d.exists() and d.is_dir():
            source_dirs.append(d)
    preferred = {}
    ext_rank = {".cif.gz": 0, ".cif": 1, ".pdb": 2}
    for root in source_dirs:
        for path in sorted(root.glob("*")):
            if not path.is_file():
                continue
            name = path.name.lower()
            if not name.endswith((".cif", ".cif.gz", ".pdb")):
                continue
            rep = _pdbzn_guess_pdb_id(path.name)
            if not rep:
                continue
            rank = 9
            for suffix, suffix_rank in ext_rank.items():
                if name.endswith(suffix):
                    rank = suffix_rank
                    break
            old = preferred.get(rep)
            if old is None or rank < old[0]:
                preferred[rep] = (rank, path)
    out = [item[1] for item in preferred.values()]
    out.sort(key=lambda p: p.name.lower())
    return out


def _pdbzn_find_structure_path(pdb_id):
    return _pdbzn_local_structure_path(pdb_id, prefer_pdb=False)


def _pdbzn_existing_file(candidates):
    for name in candidates:
        for p in _pdbzn_candidate_files(name):
            if p.exists() and p.is_file():
                return p
    return None


def _pdbzn_file_map():
    return {
        "step1_clusters": _pdbzn_existing_file(["optimized_tmalign_clusters1.csv", "optimized_tmalign_clusters.csv"]),
        "step2_similarity": _pdbzn_existing_file(["zn_his_similarity_ranking2.csv", "zn_his_similarity_ranking.csv"]),
        "step3_ligand_filter": _pdbzn_existing_file(["zn_his_similarity_ranking_filtered23.csv", "zn_his_similarity_ranking_filtered.csv"]),
        "step4_depth_sort": _pdbzn_existing_file(["zn_his_similarity_ranking_sorted_by_depth4.csv", "zn_his_similarity_ranking_sorted_by_depth.csv"]),
        "step5_diffdock": _pdbzn_existing_file(["diffdock_screen5.csv", "diffdock_screen.csv"]),
        "step6_master": _pdbzn_existing_file(["zn_his_master_table7.csv", "zn_his_master_table.csv"]),
    }


def _pdbzn_to_num(v):
    if v is None:
        return None
    s = str(v).strip()
    if s == "":
        return None
    try:
        return float(s)
    except Exception:
        return None


def _pdbzn_header_unique(header):
    out = []
    used = {}
    for i, h in enumerate(header):
        k = (h or "").strip() or f"col_{i + 1}"
        if k not in used:
            used[k] = 1
            out.append(k)
            continue
        used[k] += 1
        out.append(f"{k}_{used[k]}")
    return out


def _pdbzn_read_table(path: Path):
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return {"columns": [], "rows": []}
    columns = _pdbzn_header_unique(rows[0])
    body = []
    n = len(columns)
    for raw in rows[1:]:
        raw = list(raw)
        if len(raw) < n:
            raw.extend([""] * (n - len(raw)))
        row = {columns[i]: raw[i] for i in range(n)}
        if len(raw) > n:
            row["__extra"] = " | ".join(raw[n:])
        body.append(row)
    return {"columns": columns, "rows": body}


def _pdbzn_rep(row):
    rep = str((row or {}).get("Representative", "") or "").strip()
    return rep.upper() if rep else ""


def _pdbzn_step_payload(step_id, title, columns, rows, before_count, row_limit):
    limited = rows[:row_limit] if row_limit > 0 else rows
    return {
        "id": step_id,
        "title": title,
        "columns": columns,
        "rows": limited,
        "total_before_filter": before_count,
        "total_after_filter": len(rows),
        "rows_returned": len(limited),
        "truncated": len(limited) < len(rows),
    }


def _pdbzn_workflow_defaults():
    return {
        "row_limit": 500,
        "step1": {
            "metal_ion": "ZN",
            "residue_requirements": "HIS:3",
            "protein_name_contains": "",
            "pdb_contains": "",
        },
        "step2": {
            "cluster_threshold": 0.7,
            "keep_only_cluster_members": True,
            "mode": "recompute_tmalign",
        },
        "step3": {
            "ligand_distance": 5.0,
            "metal_distance": 8.0,
            "max_residue_per_oligomer": 400,
        },
        "step4": {
            "cluster_threshold": 0.7,
            "require_his3": True,
        },
        "step5": {
            "run_fpocket": True,
            "fpocket_max_runs": 20,
            "fpocket_timeout_sec": 120,
            "write_to_master": False,
            "table_label": "",
            "output_file": "zn_his_master_table_step5.csv",
        },
    }


def _pdbzn_guess_pdb_id(file_name: str):
    stem = Path(file_name or "").stem
    if stem.lower().endswith(".cif"):
        stem = Path(stem).stem
    m = re.search(r"([0-9A-Za-z]{4})", stem or "")
    if m:
        return m.group(1).upper()
    return (stem or "").strip().upper()


def _pdbzn_parse_neighbor_counts(neighbor_list):
    text = str(neighbor_list or "").strip()
    out = {}
    if not text:
        return out
    for item in text.split(";"):
        seg = item.split(":")
        if len(seg) < 2:
            continue
        res = str(seg[1] or "").strip().upper()
        if not res:
            continue
        out[res] = int(out.get(res, 0)) + 1
    return out


def _pdbzn_parse_residue_requirements(raw):
    out = {}
    if isinstance(raw, dict):
        for k, v in raw.items():
            key = str(k or "").strip().upper()
            n = _pdbzn_to_num(v)
            if key and n is not None and n > 0:
                out[key] = int(n)
        return out
    txt = str(raw or "").replace("，", ",").replace("；", ";").strip()
    if not txt:
        return out
    for seg in re.split(r"[,\n;]+", txt):
        token = seg.strip()
        if not token:
            continue
        pair = re.split(r"[:= ]+", token)
        if len(pair) < 2:
            continue
        key = str(pair[0] or "").strip().upper()
        n = _pdbzn_to_num(pair[1])
        if key and n is not None and n > 0:
            out[key] = int(n)
    return out


def _pdbzn_neighbor_similarity(a, b):
    ak = set((a or {}).keys())
    bk = set((b or {}).keys())
    keys = sorted(ak | bk)
    if not keys:
        return 0.0
    num = 0.0
    den = 0.0
    for k in keys:
        av = max(0.0, float((a or {}).get(k, 0) or 0))
        bv = max(0.0, float((b or {}).get(k, 0) or 0))
        num += min(av, bv)
        den += max(av, bv)
    if den <= 0:
        return 0.0
    return max(0.0, min(1.0, num / den))


def _pdbzn_cluster_rows(rows, threshold):
    arr = list(rows or [])
    n = len(arr)
    if n < 1:
        return []
    t = _pdbzn_to_num(threshold)
    if t is None:
        t = 0.7
    t = max(0.0, min(1.0, float(t)))
    parent = list(range(n))
    rank = [0] * n

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra = find(a)
        rb = find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
            return
        if rank[ra] > rank[rb]:
            parent[rb] = ra
            return
        parent[rb] = ra
        rank[ra] += 1

    for i in range(n):
        ci = arr[i].get("_neighbor_counts") or {}
        for j in range(i + 1, n):
            cj = arr[j].get("_neighbor_counts") or {}
            sim = _pdbzn_neighbor_similarity(ci, cj)
            if sim >= t:
                union(i, j)

    groups = {}
    for i in range(n):
        r = find(i)
        groups.setdefault(r, []).append(i)
    group_list = sorted(groups.values(), key=lambda idxs: (-len(idxs), min(idxs)))
    out = []
    cid = 0
    for idxs in group_list:
        cid += 1
        members = [str(arr[k].get("PDB_ID", "")) for k in idxs if str(arr[k].get("PDB_ID", ""))]
        if not members:
            continue
        leader = sorted([(arr[k].get("Similarity_Score"), arr[k].get("PDB_ID", "")) for k in idxs], key=lambda x: (-(x[0] if x[0] is not None else -999.0), str(x[1] or "")))[0]
        sims = [arr[k].get("Similarity_Score") for k in idxs if arr[k].get("Similarity_Score") is not None]
        mean_sim = None if not sims else (sum(float(x) for x in sims) / len(sims))
        out.append({
            "Cluster_ID": cid,
            "Size": len(members),
            "Representative": str(leader[1] or ""),
            "Members": ",".join(members),
            "Mean_Similarity_Score": None if mean_sim is None else round(float(mean_sim), 6),
        })
    return out


def _pdbzn_cluster_rows_from_step1_clusters(rows, keep_only_members=True):
    arr = list(rows or [])
    if not arr:
        return {"clusters": [], "outside_members": []}
    f_clusters = _pdbzn_existing_file(["optimized_tmalign_clusters1.csv", "optimized_tmalign_clusters.csv"])
    if not f_clusters or (not f_clusters.exists()):
        return {"clusters": [], "outside_members": []}
    table = _pdbzn_read_table(f_clusters)
    member_to_group = {}
    for idx, r in enumerate(table.get("rows", [])):
        cid_raw = str((r or {}).get("Cluster_ID", "") or "").strip() or str(idx + 1)
        rep = _pdbzn_guess_pdb_id(str((r or {}).get("Representative", "") or ""))
        members_raw = str((r or {}).get("Members", "") or "")
        tokens = [x for x in re.split(r"[\s,;|]+", members_raw) if str(x or "").strip()]
        members = []
        if rep:
            members.append(rep)
        for t in tokens:
            m = _pdbzn_guess_pdb_id(t)
            if m:
                members.append(m)
        dedup = []
        seen = set()
        for m in members:
            if m in seen:
                continue
            seen.add(m)
            dedup.append(m)
        for m in dedup:
            member_to_group[m] = {
                "source_order": idx,
                "source_cluster_id": cid_raw,
                "source_representative": rep,
            }
    if not member_to_group:
        return {"clusters": [], "outside_members": []}
    row_map = {str(r.get("PDB_ID", "")).upper(): r for r in arr if str(r.get("PDB_ID", ""))}
    grouped = {}
    grouped_meta = {}
    outside_members = []
    for pid in row_map.keys():
        g = member_to_group.get(pid)
        if g:
            key = f"M:{g['source_order']}:{g['source_cluster_id']}"
            grouped_meta[key] = g
        else:
            outside_members.append(pid)
            if keep_only_members:
                continue
            key = f"S:{pid}"
            grouped_meta[key] = {
                "source_order": 10**9,
                "source_cluster_id": "",
                "source_representative": "",
            }
        grouped.setdefault(key, []).append(pid)
    ordered_keys = sorted(
        grouped.keys(),
        key=lambda k: (
            -len(grouped.get(k, [])),
            int(grouped_meta.get(k, {}).get("source_order", 10**9)),
            k,
        ),
    )
    out = []
    cid = 0
    for k in ordered_keys:
        members = sorted(grouped.get(k, []))
        if not members:
            continue
        cid += 1
        g = grouped_meta.get(k, {})
        rep = str(g.get("source_representative", "") or "").upper()
        if rep not in members:
            best = sorted(
                [
                    (
                        row_map[m].get("Similarity_Score"),
                        m,
                    )
                    for m in members
                ],
                key=lambda x: (-(x[0] if x[0] is not None else -999.0), str(x[1] or "")),
            )[0]
            rep = str(best[1] or "")
        sims = [row_map[m].get("Similarity_Score") for m in members if row_map[m].get("Similarity_Score") is not None]
        mean_sim = None if not sims else (sum(float(x) for x in sims) / len(sims))
        out.append(
            {
                "Cluster_ID": cid,
                "Size": len(members),
                "Representative": rep,
                "Members": ",".join(members),
                "Mean_Similarity_Score": None if mean_sim is None else round(float(mean_sim), 6),
            }
        )
    return {"clusters": out, "outside_members": sorted(outside_members)}


def _pdbzn_parse_len_oligomer_monomer(raw):
    txt = str(raw or "").strip()
    if not txt:
        return (None, None, None)
    parts = [p.strip() for p in txt.split(";")]
    vals = []
    for p in parts[:3]:
        n = _pdbzn_to_num(p)
        vals.append(None if n is None else float(n))
    while len(vals) < 3:
        vals.append(None)
    return (vals[0], vals[1], vals[2])


def _pdbzn_find_cif_path(pdb_id):
    return _pdbzn_find_structure_path(pdb_id)


def _pdbzn_mmcif_atom_rows(cif_path: Path):
    rows = []
    open_fn = gzip.open if str(cif_path).lower().endswith(".gz") else open
    headers = []
    in_atom_loop = False
    with open_fn(cif_path, "rt", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = (raw or "").strip()
            if not line:
                continue
            if line == "loop_":
                in_atom_loop = False
                headers = []
                continue
            if line.startswith("_atom_site."):
                in_atom_loop = True
                headers.append(line)
                continue
            if in_atom_loop and line.startswith("_"):
                in_atom_loop = False
                headers = []
                continue
            if in_atom_loop and line.startswith("#"):
                in_atom_loop = False
                headers = []
                continue
            if not in_atom_loop:
                continue
            if not headers:
                continue
            try:
                vals = shlex.split(line)
            except Exception:
                vals = line.split()
            if len(vals) < len(headers):
                continue
            m = {headers[i]: vals[i] for i in range(len(headers))}
            comp = str(m.get("_atom_site.auth_comp_id", "") or m.get("_atom_site.label_comp_id", "") or "").strip().upper()
            atom_name = str(m.get("_atom_site.auth_atom_id", "") or m.get("_atom_site.label_atom_id", "") or "").strip().upper()
            element = str(m.get("_atom_site.type_symbol", "") or "").strip().upper()
            if not element:
                letters = re.sub(r"[^A-Za-z]+", "", atom_name or "")
                if letters:
                    element = letters[:2].upper()
            x = _pdbzn_to_num(m.get("_atom_site.Cartn_x"))
            y = _pdbzn_to_num(m.get("_atom_site.Cartn_y"))
            z = _pdbzn_to_num(m.get("_atom_site.Cartn_z"))
            if x is None or y is None or z is None:
                continue
            chain = str(m.get("_atom_site.auth_asym_id", "") or m.get("_atom_site.label_asym_id", "") or "").strip()
            seq = str(m.get("_atom_site.auth_seq_id", "") or m.get("_atom_site.label_seq_id", "") or "").strip()
            icode = str(m.get("_atom_site.pdbx_PDB_ins_code", "") or "").strip()
            rows.append(
                {
                    "comp": comp,
                    "atom": atom_name,
                    "element": element,
                    "x": float(x),
                    "y": float(y),
                    "z": float(z),
                    "chain": chain,
                    "seq": seq,
                    "icode": icode,
                }
            )
    return rows


def _pdbzn_pdb_atom_rows(pdb_path: Path):
    rows = []
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = (raw or "").rstrip("\n")
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = str(line[12:16] or "").strip().upper()
            comp = str(line[17:20] or "").strip().upper()
            chain = str(line[21:22] or "").strip()
            seq = str(line[22:26] or "").strip()
            icode = str(line[26:27] or "").strip()
            try:
                x = float(str(line[30:38] or "").strip())
                y = float(str(line[38:46] or "").strip())
                z = float(str(line[46:54] or "").strip())
            except Exception:
                continue
            element = str(line[76:78] or "").strip().upper()
            if not element:
                letters = re.sub(r"[^A-Za-z]+", "", atom_name or "")
                if letters:
                    element = letters[:2].upper()
            rows.append(
                {
                    "comp": comp,
                    "atom": atom_name,
                    "element": element,
                    "x": float(x),
                    "y": float(y),
                    "z": float(z),
                    "chain": chain,
                    "seq": seq,
                    "icode": icode,
                }
            )
    return rows


def _pdbzn_structure_atom_rows(structure_path: Path):
    suffixes = [s.lower() for s in structure_path.suffixes]
    if suffixes[-2:] == [".cif", ".gz"] or (suffixes and suffixes[-1] == ".cif"):
        return _pdbzn_mmcif_atom_rows(structure_path)
    if suffixes and suffixes[-1] == ".pdb":
        return _pdbzn_pdb_atom_rows(structure_path)
    return _pdbzn_mmcif_atom_rows(structure_path)


def _pdbzn_dist(a, b):
    dx = float(a["x"]) - float(b["x"])
    dy = float(a["y"]) - float(b["y"])
    dz = float(a["z"]) - float(b["z"])
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _pdbzn_residue_key(atom):
    chain = str((atom or {}).get("chain", "") or "").strip()
    seq = str((atom or {}).get("seq", "") or "").strip()
    icode = str((atom or {}).get("icode", "") or "").strip()
    comp = str((atom or {}).get("comp", "") or "").strip().upper()
    return (chain, seq, icode, comp)


def _pdbzn_residue_label(res_key):
    chain, seq, icode, comp = res_key
    c = chain or "?"
    s = seq or "?"
    ic = "" if (not icode or icode in {"?", "."}) else icode
    return f"{comp}:{c}:{s}{ic}"


def _pdbzn_is_truthy(value):
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    return text in {"1", "true", "yes", "y", "ok"}


def _pdbzn_coord_donor_atom_allowed(atom):
    comp = str((atom or {}).get("comp", "") or "").strip().upper()
    atom_name = str((atom or {}).get("atom", "") or "").strip().upper()
    return atom_name in PDBZN_PROTEIN_DONOR_ATOMS.get(comp, set())


def _pdbzn_coord_site_sort_key(site):
    return (
        -int(site.get("his_count", 0)),
        int(site.get("non_his_count", 10**9)),
        int(site.get("residue_count", 10**9)),
        float(site.get("dist_sum", 9999.0)),
        str(site.get("site", "")),
        int(site.get("index", 10**9)),
    )


def _pdbzn_coord_site_summary(site):
    if not site:
        return ""
    return (
        f"Zn{site.get('index', '?')}@{site.get('site', '?')}: "
        f"total={site.get('residue_count', 0)}, "
        f"his={site.get('his_count', 0)}, "
        f"non_his={site.get('non_his_count', 0)}, "
        f"residues={site.get('residues_text', '') or '-'}"
    )


def _pdbzn_collect_coord_sites(atoms):
    zn_atoms = [
        a for a in atoms
        if str(a.get("comp") or "").upper() == "ZN" or str(a.get("element") or "").upper() == "ZN"
    ]
    donor_atoms = [a for a in atoms if _pdbzn_coord_donor_atom_allowed(a)]
    sites = []
    for idx, zn in enumerate(zn_atoms, start=1):
        nearest = {}
        for atom in donor_atoms:
            d = _pdbzn_dist(zn, atom)
            if d > 3.2:
                continue
            rk = _pdbzn_residue_key(atom)
            old = nearest.get(rk)
            if old is None or d < old["distance"]:
                nearest[rk] = {"atom": atom, "distance": float(d)}
        his_count = sum(1 for rk in nearest if rk[3] in PDBZN_HIS_RESIDUES)
        total_count = len(nearest)
        dist_sum = sum(sorted(x["distance"] for x in nearest.values())[:3]) if nearest else 9999.0
        ordered = sorted(nearest.items(), key=lambda kv: (kv[1]["distance"], _pdbzn_residue_label(kv[0])))
        residue_labels = [f"{_pdbzn_residue_label(rk)}:{round(item['distance'], 3)}" for rk, item in ordered]
        sites.append(
            {
                "index": idx,
                "his_count": his_count,
                "residue_count": total_count,
                "non_his_count": max(0, total_count - his_count),
                "dist_sum": dist_sum,
                "site": f"{str(zn.get('chain') or '?')}:{str(zn.get('seq') or '?')}",
                "residues_text": "; ".join(residue_labels),
            }
        )
    return sorted(sites, key=_pdbzn_coord_site_sort_key)


def _pdbzn_coordination_metrics(rep, receptor_file=None):
    fp = _safe_existing_file(receptor_file) if receptor_file else None
    if fp is None:
        fp = _pdbzn_find_tmalign_pdb_path(rep)
    if fp is None:
        return None
    try:
        atoms = _pdbzn_structure_atom_rows(fp)
    except Exception:
        return None
    if not atoms:
        return None
    sites = _pdbzn_collect_coord_sites(atoms)
    if not sites:
        return None
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
            non_his_value = best.get("non_his_count", "")
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
            best = sorted(sites, key=_pdbzn_coord_site_sort_key)[0]
            non_his_value = ""
            note_sites = sorted(sites, key=_pdbzn_coord_site_sort_key)
    else:
        best = sorted(sites, key=_pdbzn_coord_site_sort_key)[0]
        non_his_value = best.get("non_his_count", "")
        note_sites = sites
    others = [s for s in note_sites if s is not best]
    return {
        "Zn_CoordResidueCount": best.get("residue_count", 0),
        "Zn_CoordHisCount": best.get("his_count", 0),
        "Zn_CoordNonHisCount": non_his_value,
        "Zn_CoordResidues": best.get("residues_text", ""),
        "Zn_CoordSite": best.get("site", ""),
        "Zn_CoordNote": " | ".join(_pdbzn_coord_site_summary(site) for site in others),
    }


def _pdbzn_angle_deg(v1, v2):
    n1 = math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])
    n2 = math.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])
    if n1 <= 0.0 or n2 <= 0.0:
        return None
    dot = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2)
    dot = max(-1.0, min(1.0, float(dot)))
    return math.degrees(math.acos(dot))


def _pdbzn_step4_geometry_report(pdb_id, radius=5.0):
    fp = _pdbzn_find_cif_path(pdb_id)
    if not fp:
        return {
            "ok": False,
            "reasons": [f"缺少结构文件 {pdb_id}"],
            "tri_his_count": 0,
        }
    try:
        atoms = _pdbzn_structure_atom_rows(fp)
    except Exception:
        return {
            "ok": False,
            "reasons": [f"结构解析失败 {pdb_id}"],
            "tri_his_count": 0,
        }
    if not atoms:
        return {
            "ok": False,
            "reasons": [f"结构无坐标原子 {pdb_id}"],
            "tri_his_count": 0,
        }
    zn_atoms = [a for a in atoms if (a.get("comp") == "ZN" or a.get("element") == "ZN")]
    if not zn_atoms:
        return {
            "ok": False,
            "reasons": ["未检测到ZN原子"],
            "tri_his_count": 0,
        }
    his_atoms = [
        a
        for a in atoms
        if str(a.get("comp") or "").upper() in {"HIS", "HID", "HIE", "HIP"}
        and str(a.get("atom") or "").upper() in {"ND1", "NE2"}
    ]
    best_site = None
    for z in zn_atoms:
        nearest_by_res = {}
        for ha in his_atoms:
            d = _pdbzn_dist(z, ha)
            if d > 3.2:
                continue
            rk = _pdbzn_residue_key(ha)
            old = nearest_by_res.get(rk)
            if (old is None) or (d < old["distance"]):
                nearest_by_res[rk] = {"atom": ha, "distance": float(d)}
        picked = sorted(nearest_by_res.items(), key=lambda kv: kv[1]["distance"])
        tri = picked[:3]
        score_count = len(picked)
        score_sum3 = sum(x[1]["distance"] for x in tri) if len(tri) >= 3 else 9999.0
        candidate = {
            "zn": z,
            "picked": picked,
            "tri": tri,
            "score_count": score_count,
            "score_sum3": score_sum3,
        }
        if best_site is None:
            best_site = candidate
        else:
            better = (candidate["score_count"] > best_site["score_count"]) or (
                candidate["score_count"] == best_site["score_count"] and candidate["score_sum3"] < best_site["score_sum3"]
            )
            if better:
                best_site = candidate
    if best_site is None:
        return {
            "ok": False,
            "reasons": ["未找到可用ZN位点"],
            "tri_his_count": 0,
        }
    tri = list(best_site["tri"])
    tri_his_count = len(tri)
    reasons = []
    if tri_his_count < 3:
        reasons.append(f"ZN周围可配位组氨酸不足3个（检测到{tri_his_count}）")
    dist_items = []
    dist_vals = []
    labels = []
    coords = []
    for rk, item in tri:
        a = item["atom"]
        d = float(item["distance"])
        lab = _pdbzn_residue_label(rk)
        labels.append(lab)
        dist_vals.append(d)
        coords.append((float(a["x"]), float(a["y"]), float(a["z"])))
        dist_items.append(f"{lab}-{str(a.get('atom') or '').upper()}:{round(d, 3)}")
    if len(set((round(x[0], 3), round(x[1], 3), round(x[2], 3)) for x in coords)) < len(coords):
        reasons.append("疑似对称构象重复原子（坐标重合）")
    for d in dist_vals:
        if d < 1.7 or d > 3.1:
            reasons.append(f"HIS-ZN键长异常（{round(d, 3)}A）")
            break
    his_pair_min = None
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            pi = {"x": coords[i][0], "y": coords[i][1], "z": coords[i][2]}
            pj = {"x": coords[j][0], "y": coords[j][1], "z": coords[j][2]}
            dd = _pdbzn_dist(pi, pj)
            if his_pair_min is None or dd < his_pair_min:
                his_pair_min = float(dd)
    if his_pair_min is not None and his_pair_min < 0.8:
        reasons.append(f"疑似对称构象重叠（HIS-HIS最小距离 {round(his_pair_min, 3)}A）")
    angle_items = []
    z = best_site["zn"]
    v = []
    for rk, item in tri:
        a = item["atom"]
        v.append((float(a["x"]) - float(z["x"]), float(a["y"]) - float(z["y"]), float(a["z"]) - float(z["z"])))
    angle_vals = []
    for i in range(len(v)):
        for j in range(i + 1, len(v)):
            ang = _pdbzn_angle_deg(v[i], v[j])
            if ang is None:
                continue
            angle_vals.append(float(ang))
            li = labels[i] if i < len(labels) else f"HIS{i+1}"
            lj = labels[j] if j < len(labels) else f"HIS{j+1}"
            angle_items.append(f"{li}~{lj}:{round(float(ang), 2)}")
    for a in angle_vals:
        if a < 70.0 or a > 170.0:
            reasons.append(f"配位角异常（{round(a, 2)}°）")
            break
    close_res = {}
    for a in atoms:
        d = _pdbzn_dist(z, a)
        if d > float(radius):
            continue
        rk = _pdbzn_residue_key(a)
        old = close_res.get(rk)
        if old is None or d < old:
            close_res[rk] = float(d)
    close_items = sorted(close_res.items(), key=lambda kv: kv[1])
    close_labels = [f"{_pdbzn_residue_label(k)}:{round(vv, 3)}" for k, vv in close_items]
    zn_site = f"{str(z.get('chain') or '?')}:{str(z.get('seq') or '?')}"
    return {
        "ok": len(reasons) == 0 and tri_his_count >= 3,
        "reasons": reasons,
        "tri_his_count": tri_his_count,
        "tri_his_residues": "; ".join(labels),
        "his_zn_bond_lengths": "; ".join(dist_items),
        "his_zn_bond_angles": "; ".join(angle_items),
        "residues_within_5a": "; ".join(close_labels),
        "residues_within_5a_count": len(close_labels),
        "zn_site": zn_site,
        "symmetry_spatial_filter_pass": (len(reasons) == 0),
    }


def _pdbzn_step3_structure_reasons(pdb_id, ligand_distance, metal_distance):
    fp = _pdbzn_find_cif_path(pdb_id)
    if not fp:
        return [f"缺少结构文件 {pdb_id}"], {"max_his_coord": None}
    try:
        atoms = _pdbzn_structure_atom_rows(fp)
    except Exception:
        return [f"结构解析失败 {pdb_id}"], {"max_his_coord": None}
    if not atoms:
        return [f"结构无坐标原子 {pdb_id}"], {"max_his_coord": None}
    zn_atoms = [a for a in atoms if (a.get("comp") == "ZN" or a.get("element") == "ZN")]
    ligand_atoms = [a for a in atoms if str(a.get("comp") or "").upper() in PDBZN_LARGE_LIGANDS]
    other_metal_atoms = [a for a in atoms if str(a.get("element") or "").upper() in PDBZN_METAL_ELEMENTS and str(a.get("comp") or "").upper() != "ZN"]
    his_atoms = [a for a in atoms if str(a.get("comp") or "").upper() in {"HIS", "HID", "HIE", "HIP"} and str(a.get("atom") or "").upper() in {"ND1", "NE2"}]
    reasons = []
    for z in zn_atoms:
        for la in ligand_atoms:
            d = _pdbzn_dist(z, la)
            if d <= float(ligand_distance):
                reasons.append(f"ZN附近{ligand_distance}A内存在大配体 {la.get('comp')} ({round(d, 2)}A)")
                break
        if reasons:
            break
    for z in zn_atoms:
        for ma in other_metal_atoms:
            d = _pdbzn_dist(z, ma)
            if d <= float(metal_distance):
                reasons.append(f"ZN附近{metal_distance}A内存在其他金属 {ma.get('comp') or ma.get('element')} ({round(d, 2)}A)")
                break
        if any("其他金属" in x for x in reasons):
            break
    max_his = 0
    for z in zn_atoms:
        coord_set = set()
        for ha in his_atoms:
            if _pdbzn_dist(z, ha) <= 3.0:
                coord_set.add((ha.get("chain"), ha.get("seq"), ha.get("icode")))
        max_his = max(max_his, len(coord_set))
    if zn_atoms and max_his < 3:
        reasons.append(f"ZN周围配位组氨酸不足3个（检测到{max_his}）")
    return reasons, {"max_his_coord": max_his}


def _pdbzn_find_tmalign_pdb_path(pdb_id):
    return _pdbzn_local_structure_path(pdb_id, prefer_pdb=True)


def _pdbzn_tmalign_score(p1: Path, p2: Path):
    exe = discover_tmalign_bin()
    if not exe:
        return None
    try:
        a = str(Path(p1).resolve())
        b = str(Path(p2).resolve())
        if a > b:
            a, b = b, a
        key = (
            a,
            int(Path(a).stat().st_mtime),
            b,
            int(Path(b).stat().st_mtime),
            str(exe),
        )
    except Exception:
        key = None
    if key is not None:
        with _pdbzn_tmalign_cache_lock:
            cached = _pdbzn_tmalign_cache.get(key)
        if cached is not None:
            return cached
    try:
        proc = subprocess.run(
            [exe, str(p1), str(p2)],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=20,
        )
        out = str(proc.stdout or "")
        scores = []
        for m in re.finditer(r"TM-score=\s*([0-9]*\.?[0-9]+)", out):
            try:
                scores.append(float(m.group(1)))
            except Exception:
                pass
        if not scores:
            score = 0.0
        else:
            score = float(max(scores))
        if key is not None:
            with _pdbzn_tmalign_cache_lock:
                _pdbzn_tmalign_cache[key] = score
        return score
    except Exception:
        return 0.0


def _pdbzn_cluster_rows_tmalign(rows, threshold):
    arr = list(rows or [])
    n = len(arr)
    if n < 1:
        return []
    t = _pdbzn_to_num(threshold)
    if t is None:
        t = 0.7
    t = max(0.0, min(1.0, float(t)))
    parent = list(range(n))
    rank = [0] * n

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra = find(a)
        rb = find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
            return
        if rank[ra] > rank[rb]:
            parent[rb] = ra
            return
        parent[rb] = ra
        rank[ra] += 1

    valid = []
    for i, r in enumerate(arr):
        pid = str(r.get("Representative", "") or r.get("PDB_ID", "") or "").strip().upper()
        fp = _pdbzn_find_tmalign_pdb_path(pid)
        if fp:
            valid.append((i, pid, fp))
    for a in range(len(valid)):
        ia, _, f1 = valid[a]
        for b in range(a + 1, len(valid)):
            ib, _, f2 = valid[b]
            s = _pdbzn_tmalign_score(f1, f2)
            if s is not None and float(s) > t:
                union(ia, ib)

    groups = {}
    for i in range(n):
        r = find(i)
        groups.setdefault(r, []).append(i)
    group_list = sorted(groups.values(), key=lambda idxs: (-len(idxs), min(idxs)))
    out = []
    cid = 0
    for idxs in group_list:
        cid += 1
        members = []
        for k in idxs:
            pid = str(arr[k].get("Representative", "") or arr[k].get("PDB_ID", "") or "").strip().upper()
            if pid:
                members.append(pid)
        if not members:
            continue
        leader = sorted(
            [(arr[k].get("Similarity_Score"), str(arr[k].get("Representative", "") or arr[k].get("PDB_ID", "") or "").strip().upper()) for k in idxs],
            key=lambda x: (-(x[0] if x[0] is not None else -999.0), x[1]),
        )[0]
        rep = str(leader[1] or members[0] or "")
        out.append(
            {
                "Cluster_ID": cid,
                "Size": len(members),
                "Representative": rep,
                "Members": ",".join(members),
                "TMAlign_Threshold": t,
            }
        )
    return out


def _pdbzn_db_connect():
    conn = sqlite3.connect(str(PDBZN_DB_PATH))
    conn.row_factory = sqlite3.Row
    return conn


def _pdbzn_db_ensure(conn):
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS proteins (
            pdb_id TEXT PRIMARY KEY,
            file_name TEXT NOT NULL,
            file_path TEXT NOT NULL,
            has_zn_his_cluster INTEGER NOT NULL DEFAULT 0,
            metal_ions TEXT NOT NULL DEFAULT '',
            his_count_max INTEGER NOT NULL DEFAULT 0,
            neighbor_residue_counts TEXT NOT NULL DEFAULT '{}',
            similarity_score REAL,
            protein_name TEXT,
            protein_category TEXT,
            details TEXT,
            imported_at INTEGER NOT NULL DEFAULT 0
        )
        """
    )
    conn.commit()


def _pdbzn_db_stats():
    if not PDBZN_DB_PATH.exists():
        return {"db_exists": False, "total_proteins": 0, "with_master": 0, "last_imported_at": 0}
    conn = _pdbzn_db_connect()
    try:
        _pdbzn_db_ensure(conn)
        total = int((conn.execute("SELECT COUNT(*) AS c FROM proteins").fetchone() or {"c": 0})["c"])
        with_master = int((conn.execute("SELECT COUNT(*) AS c FROM proteins WHERE protein_name IS NOT NULL AND TRIM(protein_name) != ''").fetchone() or {"c": 0})["c"])
        imported_at = int((conn.execute("SELECT COALESCE(MAX(imported_at), 0) AS t FROM proteins").fetchone() or {"t": 0})["t"])
        return {"db_exists": True, "total_proteins": total, "with_master": with_master, "last_imported_at": imported_at}
    finally:
        conn.close()


def _pdbzn_source_total_count():
    try:
        return len(_pdbzn_source_structure_files())
    except Exception:
        return 0


def _pdbzn_import_database(import_options=None):
    options = import_options if isinstance(import_options, dict) else {}
    upload_ids = []
    if isinstance(options.get("pdb_ids"), list):
        for item in options.get("pdb_ids") or []:
            pid = _pdbzn_guess_pdb_id(str(item or ""))
            if pid and pid not in upload_ids:
                upload_ids.append(pid)
    if not upload_ids:
        upload_ids = _pdbzn_parse_uploaded_ids(options.get("pdb_list_text", ""))
    source_name = str(options.get("source_name", "") or "").strip()
    clear_selection = bool(options.get("clear_selection", False))
    if clear_selection:
        _pdbzn_clear_import_selection()
    if upload_ids:
        all_map, missing = _pdbzn_resolve_import_file_map(upload_ids, allow_download=True)
        if missing:
            return {
                "ok": False,
                "error": "部分 PDB 结构无法从本地库或 RCSB 获取",
                "missing_ids": missing,
                "resolved_count": len(all_map),
            }
        selection = _pdbzn_write_import_selection(upload_ids, source_name)
        source_mode = "upload_list"
    else:
        active_selection = _pdbzn_active_import_selection()
        if active_selection and active_selection.get("ids"):
            all_map, missing = _pdbzn_resolve_import_file_map(active_selection.get("ids") or [], allow_download=True)
            if missing:
                return {
                    "ok": False,
                    "error": "当前上传列表中存在无法解析的结构文件",
                    "missing_ids": missing,
                    "resolved_count": len(all_map),
                }
            selection = active_selection
            source_mode = "saved_upload_list"
        else:
            files_all = list(_pdbzn_source_structure_files())
            if not files_all:
                return {"ok": False, "error": f"未找到可导入的结构文件（.pdb/.cif）: {PDBZN_PDB_DIR}"}
            all_map = {(_pdbzn_guess_pdb_id(p.name)): p for p in files_all}
            selection = None
            source_mode = "all_structure_files"
    if not all_map:
        return {"ok": False, "error": "入库范围为空：未解析到可用结构文件"}
    cluster_flags = {}
    f_opt = _pdbzn_existing_file(["optimized-main1.csv"])
    if f_opt and f_opt.exists():
        t_opt = _pdbzn_read_table(f_opt)
        for r in t_opt["rows"]:
            fn = str(r.get("Filename", "") or "").strip()
            if not fn:
                continue
            rep = _pdbzn_guess_pdb_id(fn)
            flag = str(r.get("Has_ZN_HIS_Cluster", "") or "").strip().upper()
            cluster_flags[rep] = 1 if flag in {"YES", "Y", "TRUE", "1"} else 0
    file_map = dict(all_map)
    master_rows = {}
    f_master = _pdbzn_existing_file(["zn_his_master_table7.csv", "zn_his_master_table.csv"])
    if f_master and f_master.exists():
        t_master = _pdbzn_read_table(f_master)
        for r in t_master["rows"]:
            rep = _pdbzn_rep(r)
            if rep:
                master_rows[rep] = r
    import_at = now_ts()
    inserted = 0
    updated = 0
    deleted = 0
    with pdbzn_db_lock:
        conn = _pdbzn_db_connect()
        try:
            _pdbzn_db_ensure(conn)
            keep_ids = sorted(file_map.keys())
            if keep_ids:
                q = ",".join(["?"] * len(keep_ids))
                cur = conn.execute(f"DELETE FROM proteins WHERE pdb_id NOT IN ({q})", tuple(keep_ids))
                deleted = int(cur.rowcount or 0)
            for rep, fp in file_map.items():
                mr = master_rows.get(rep, {})
                neighbor_counts = _pdbzn_parse_neighbor_counts((mr or {}).get("Neighbor_List", ""))
                his_neighbor = int(neighbor_counts.get("HIS", 0))
                his_tri = _pdbzn_to_num((mr or {}).get("TriHisCountMax"))
                his_max = his_neighbor if his_tri is None else max(his_neighbor, int(his_tri))
                if (his_max < 3) and int(cluster_flags.get(rep, 0)) == 1:
                    his_max = 3
                sim = _pdbzn_to_num((mr or {}).get("Similarity_Score"))
                exists = conn.execute("SELECT 1 FROM proteins WHERE pdb_id = ?", (rep,)).fetchone()
                conn.execute(
                    """
                    INSERT INTO proteins (
                        pdb_id, file_name, file_path, has_zn_his_cluster, metal_ions, his_count_max,
                        neighbor_residue_counts, similarity_score, protein_name, protein_category, details, imported_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(pdb_id) DO UPDATE SET
                        file_name=excluded.file_name,
                        file_path=excluded.file_path,
                        has_zn_his_cluster=excluded.has_zn_his_cluster,
                        metal_ions=excluded.metal_ions,
                        his_count_max=excluded.his_count_max,
                        neighbor_residue_counts=excluded.neighbor_residue_counts,
                        similarity_score=excluded.similarity_score,
                        protein_name=excluded.protein_name,
                        protein_category=excluded.protein_category,
                        details=excluded.details,
                        imported_at=excluded.imported_at
                    """,
                    (
                        rep,
                        fp.name,
                        str(fp),
                        int(cluster_flags.get(rep, 0)),
                        "ZN",
                        int(his_max),
                        json.dumps(neighbor_counts, ensure_ascii=False),
                        sim,
                        str((mr or {}).get("Protein_Name", "") or ""),
                        str((mr or {}).get("Protein_Category", "") or ""),
                        str((mr or {}).get("Details", "") or ""),
                        import_at,
                    ),
                )
                if exists:
                    updated += 1
                else:
                    inserted += 1
            conn.commit()
        finally:
            conn.close()
    stats = _pdbzn_db_stats()
    return {
        "ok": True,
        "db_path": str(PDBZN_DB_PATH),
        "imported_files": len(file_map),
        "inserted": inserted,
        "updated": updated,
        "deleted": deleted,
        "source_mode": source_mode,
        "source_dir": str(PDBZN_BASE_DIR),
        "selection": selection,
        "selection_count": len(selection.get("ids") or []) if isinstance(selection, dict) else 0,
        "missing_ids": [],
        "stats": stats,
    }


def _pdbzn_run_workflow(filters):
    f = filters or {}
    defaults = _pdbzn_workflow_defaults()
    stats_now = _pdbzn_db_stats()
    source_total = _pdbzn_source_total_count()
    if (not PDBZN_DB_PATH.exists()) or (stats_now.get("total_proteins", 0) < 1) or (source_total > 0 and stats_now.get("total_proteins", 0) < source_total):
        imported = _pdbzn_import_database()
        if not imported.get("ok"):
            return imported
    try:
        row_limit = int(f.get("row_limit", defaults["row_limit"]) or defaults["row_limit"])
    except Exception:
        row_limit = defaults["row_limit"]
    if row_limit < 1:
        row_limit = defaults["row_limit"]
    step1f = f.get("step1", {}) or {}
    filter_used = {**defaults["step1"], **step1f}
    metal_ion = str(filter_used.get("metal_ion", "ZN") or "").strip().upper()
    residue_requirements = _pdbzn_parse_residue_requirements(filter_used.get("residue_requirements", ""))
    pname_q = str(filter_used.get("protein_name_contains", "") or "").strip().lower()
    pdb_q = str(filter_used.get("pdb_contains", "") or "").strip().lower()
    conn = _pdbzn_db_connect()
    try:
        _pdbzn_db_ensure(conn)
        rows_source = conn.execute(
            """
            SELECT
                pdb_id, file_name, has_zn_his_cluster, metal_ions, his_count_max, neighbor_residue_counts,
                similarity_score, protein_name, protein_category, details
            FROM proteins
            ORDER BY COALESCE(similarity_score, -999) DESC, pdb_id ASC
            """
        ).fetchall()
    finally:
        conn.close()
    rows1 = []
    for r in rows_source:
        metal_text = str(r["metal_ions"] or "").upper()
        if metal_ion and metal_ion not in metal_text:
            continue
        if pname_q and pname_q not in str(r["protein_name"] or "").lower():
            continue
        if pdb_q and pdb_q not in str(r["pdb_id"] or "").lower():
            continue
        try:
            neighbor_counts = json.loads(str(r["neighbor_residue_counts"] or "{}"))
            if not isinstance(neighbor_counts, dict):
                neighbor_counts = {}
        except Exception:
            neighbor_counts = {}
        his_count = int(r["his_count_max"] or 0)
        has_cluster = bool(r["has_zn_his_cluster"])
        hit = True
        for res, need in residue_requirements.items():
            have = int(neighbor_counts.get(res, 0))
            if res == "HIS":
                have = max(have, his_count)
                if has_cluster and int(need) <= 3:
                    have = max(have, 3)
            if have < int(need):
                hit = False
                break
        if not hit:
            continue
        rows1.append({
            "PDB_ID": r["pdb_id"],
            "Metal_Ions": r["metal_ions"],
            "HIS_Count": his_count,
            "Neighbor_Residue_Counts": json.dumps(neighbor_counts, ensure_ascii=False),
            "Similarity_Score": r["similarity_score"],
            "Protein_Name": r["protein_name"] or "",
            "Protein_Category": r["protein_category"] or "",
            "Has_ZN_HIS_Cluster": has_cluster,
            "Matched_By_ClusterFlag": bool(has_cluster and his_count < 3 and residue_requirements.get("HIS", 0) <= 3),
            "Details": r["details"] or "",
            "File_Name": r["file_name"],
            "_neighbor_counts": neighbor_counts,
        })
    columns = [
        "PDB_ID",
        "Metal_Ions",
        "HIS_Count",
        "Neighbor_Residue_Counts",
        "Similarity_Score",
        "Protein_Name",
        "Protein_Category",
        "Has_ZN_HIS_Cluster",
        "Matched_By_ClusterFlag",
        "Details",
        "File_Name",
    ]
    reps = [str(x.get("PDB_ID", "")) for x in rows1 if str(x.get("PDB_ID", ""))]
    rows1_public = []
    for x in rows1:
        row = dict(x)
        if "_neighbor_counts" in row:
            del row["_neighbor_counts"]
        rows1_public.append(row)
    steps = []
    step = _pdbzn_step_payload("step1_metal_residue", "Step 1: 金属离子 + 邻域残基筛选", columns, rows1_public, len(rows_source), row_limit)
    steps.append(step)
    return {
        "ok": True,
        "db_path": str(PDBZN_DB_PATH),
        "filters_used": {
            "row_limit": row_limit,
            "step1": filter_used,
            "residue_requirements_parsed": residue_requirements,
        },
        "steps": steps,
        "summary": {
            "final_count": len(rows1),
            "step1_count": len(rows1),
            "cluster_count": None,
            "final_representatives": reps[:200],
        },
    }


def _pdbzn_cluster_workflow(filters):
    f = filters or {}
    defaults = _pdbzn_workflow_defaults()
    stats_now = _pdbzn_db_stats()
    source_total = _pdbzn_source_total_count()
    if (not PDBZN_DB_PATH.exists()) or (stats_now.get("total_proteins", 0) < 1) or (source_total > 0 and stats_now.get("total_proteins", 0) < source_total):
        imported = _pdbzn_import_database()
        if not imported.get("ok"):
            return imported
    try:
        row_limit = int(f.get("row_limit", defaults["row_limit"]) or defaults["row_limit"])
    except Exception:
        row_limit = defaults["row_limit"]
    if row_limit < 1:
        row_limit = defaults["row_limit"]
    step1f = f.get("step1", {}) or {}
    filter_used = {**defaults["step1"], **step1f}
    step2f = f.get("step2", {}) or {}
    step2_used = {**defaults.get("step2", {}), **step2f}
    cluster_threshold = _pdbzn_to_num(step2_used.get("cluster_threshold", 0.7))
    keep_only_cluster_members = bool(step2_used.get("keep_only_cluster_members", True))
    cluster_mode = str(step2_used.get("mode", "recompute_tmalign") or "recompute_tmalign").strip().lower()
    if cluster_mode not in {"preset_clusters", "recompute_tmalign", "recompute_neighbor"}:
        cluster_mode = "recompute_tmalign"
    if cluster_threshold is None:
        cluster_threshold = 0.7
    cluster_threshold = max(0.0, min(1.0, float(cluster_threshold)))
    step1_result = _pdbzn_run_workflow({"row_limit": max(row_limit, 1000000), "step1": filter_used})
    if not step1_result.get("ok"):
        return step1_result
    base_rows = []
    for r in (((step1_result.get("steps") or [{}])[0].get("rows")) or []):
        try:
            neighbor_counts = json.loads(str(r.get("Neighbor_Residue_Counts") or "{}"))
            if not isinstance(neighbor_counts, dict):
                neighbor_counts = {}
        except Exception:
            neighbor_counts = {}
        row = dict(r)
        row["_neighbor_counts"] = neighbor_counts
        base_rows.append(row)
    step1_count = len(base_rows)
    rows_with_pdb = []
    no_pdb_members = []
    for r in base_rows:
        pid = str(r.get("Representative", "") or r.get("PDB_ID", "") or "").strip().upper()
        if not pid:
            continue
        fp = _pdbzn_find_tmalign_pdb_path(pid)
        if fp:
            rows_with_pdb.append(r)
        else:
            no_pdb_members.append(pid)
    dedup_no_pdb = []
    seen_no_pdb = set()
    for pid in no_pdb_members:
        if pid in seen_no_pdb:
            continue
        seen_no_pdb.add(pid)
        dedup_no_pdb.append(pid)
    outside_members = []
    clusters = []
    if cluster_mode == "recompute_tmalign":
        clusters = _pdbzn_cluster_rows_tmalign(rows_with_pdb, cluster_threshold)
    elif cluster_mode == "recompute_neighbor":
        clusters = _pdbzn_cluster_rows(rows_with_pdb, cluster_threshold)
    else:
        preset_result = _pdbzn_cluster_rows_from_step1_clusters(rows_with_pdb, keep_only_cluster_members)
        clusters = list((preset_result or {}).get("clusters") or [])
        outside_members = list((preset_result or {}).get("outside_members") or [])
        if not clusters:
            clusters = _pdbzn_cluster_rows(rows_with_pdb, cluster_threshold)
            outside_members = []
    ccols = ["Cluster_ID", "Size", "Representative", "Members", "Mean_Similarity_Score"]
    step2 = _pdbzn_step_payload("step2_clustering", "Step 2: 蛋白聚类", ccols, clusters, len(rows_with_pdb), row_limit)
    return {
        "ok": True,
        "db_path": str(PDBZN_DB_PATH),
        "filters_used": {
            "row_limit": row_limit,
            "step1": filter_used,
            "step2": {
                **step2_used,
                "keep_only_cluster_members": keep_only_cluster_members,
                "mode": cluster_mode,
            },
        },
        "steps": [step2],
        "summary": {
            "final_count": len(clusters),
            "step1_count": step1_count,
            "step2_input_count": len(rows_with_pdb),
            "cluster_count": len(clusters),
            "cluster_mode": cluster_mode,
            "dropped_outside_cluster_count": len(outside_members),
            "dropped_outside_cluster_members": outside_members[:200],
            "dropped_no_pdb_count": len(dedup_no_pdb),
            "dropped_no_pdb_members": dedup_no_pdb[:200],
            "final_representatives": [str(x.get("Representative", "")) for x in clusters if str(x.get("Representative", ""))][:200],
        },
    }


def _pdbzn_step3_filter_workflow(filters):
    f = filters or {}
    defaults = _pdbzn_workflow_defaults()
    stats_now = _pdbzn_db_stats()
    source_total = _pdbzn_source_total_count()
    if (not PDBZN_DB_PATH.exists()) or (stats_now.get("total_proteins", 0) < 1) or (source_total > 0 and stats_now.get("total_proteins", 0) < source_total):
        imported = _pdbzn_import_database()
        if not imported.get("ok"):
            return imported
    try:
        row_limit = int(f.get("row_limit", defaults["row_limit"]) or defaults["row_limit"])
    except Exception:
        row_limit = defaults["row_limit"]
    if row_limit < 1:
        row_limit = defaults["row_limit"]
    step1_used = {**defaults["step1"], **(f.get("step1", {}) or {})}
    step2_used = {**defaults["step2"], **(f.get("step2", {}) or {})}
    step3_used = {**defaults.get("step3", {}), **(f.get("step3", {}) or {})}
    ligand_distance = _pdbzn_to_num(step3_used.get("ligand_distance", 5.0))
    metal_distance = _pdbzn_to_num(step3_used.get("metal_distance", 8.0))
    max_residue_per_oligomer = _pdbzn_to_num(step3_used.get("max_residue_per_oligomer", 400))
    if ligand_distance is None:
        ligand_distance = 5.0
    if metal_distance is None:
        metal_distance = 8.0
    if max_residue_per_oligomer is None:
        max_residue_per_oligomer = 400.0
    ligand_distance = float(max(0.0, ligand_distance))
    metal_distance = float(max(0.0, metal_distance))
    max_residue_per_oligomer = float(max(1.0, max_residue_per_oligomer))
    step2_result = _pdbzn_cluster_workflow({"row_limit": max(row_limit, 1000000), "step1": step1_used, "step2": step2_used})
    if not step2_result.get("ok"):
        return step2_result
    cluster_rows = list((((step2_result.get("steps") or [{}])[0].get("rows")) or []))
    master_map = {}
    f_master = _pdbzn_existing_file(["zn_his_master_table7.csv", "zn_his_master_table.csv"])
    if f_master and f_master.exists():
        t_master = _pdbzn_read_table(f_master)
        for r in t_master.get("rows", []):
            rep = _pdbzn_rep(r)
            if rep:
                master_map[rep] = r
    kept = []
    removed = []
    for c in cluster_rows:
        rep = str(c.get("Representative", "") or "").strip().upper()
        if not rep:
            continue
        reasons = []
        mr = master_map.get(rep, {})
        l, o, m = _pdbzn_parse_len_oligomer_monomer((mr or {}).get("Length;Oligomer;Monomer"))
        per_monomer = None
        if m is not None and m > 0:
            per_monomer = float(m)
        elif (l is not None) and (o is not None) and o > 0:
            per_monomer = float(l) / float(o)
        elif l is not None:
            per_monomer = float(l)
        if (per_monomer is not None) and (per_monomer > max_residue_per_oligomer):
            reasons.append(f"残基数/几聚体={round(per_monomer, 2)} > {max_residue_per_oligomer}")
        tri_his = _pdbzn_to_num((mr or {}).get("TriHisCountMax"))
        if tri_his is not None and float(tri_his) < 3.0:
            reasons.append(f"TriHisCountMax={tri_his} < 3")
        structure_reasons, structure_stats = _pdbzn_step3_structure_reasons(rep, ligand_distance, metal_distance)
        reasons.extend(structure_reasons)
        max_his_coord = structure_stats.get("max_his_coord")
        if max_his_coord is not None and int(max_his_coord) < 3:
            if not any("组氨酸不足3个" in x for x in reasons):
                reasons.append(f"ZN周围配位组氨酸不足3个（检测到{int(max_his_coord)}）")
        row_out = {
            "Cluster_ID": c.get("Cluster_ID"),
            "Size": c.get("Size"),
            "Representative": rep,
            "Members": c.get("Members"),
            "Mean_Similarity_Score": c.get("Mean_Similarity_Score"),
            "Filter_Reasons": " | ".join(reasons),
        }
        if reasons:
            removed.append(row_out)
        else:
            kept.append(row_out)
    cols = ["Cluster_ID", "Size", "Representative", "Members", "Mean_Similarity_Score", "Filter_Reasons"]
    step_kept = _pdbzn_step_payload("step3_second_filter_kept", "Step 3: 第二次筛选（保留）", cols, kept, len(cluster_rows), row_limit)
    step_removed = _pdbzn_step_payload("step3_second_filter_removed", "Step 3: 第二次筛选（剔除）", cols, removed, len(cluster_rows), row_limit)
    return {
        "ok": True,
        "db_path": str(PDBZN_DB_PATH),
        "filters_used": {
            "row_limit": row_limit,
            "step1": step1_used,
            "step2": step2_used,
            "step3": {
                "ligand_distance": ligand_distance,
                "metal_distance": metal_distance,
                "max_residue_per_oligomer": max_residue_per_oligomer,
            },
        },
        "steps": [step_kept, step_removed],
        "summary": {
            "step1_count": int((step2_result.get("summary") or {}).get("step1_count", 0) or 0),
            "cluster_count_before_step3": len(cluster_rows),
            "removed_count": len(removed),
            "final_count": len(kept),
            "final_representatives": [str(x.get("Representative", "")) for x in kept if str(x.get("Representative", ""))][:200],
        },
    }


def _pdbzn_step4_validate_workflow(filters):
    f = filters or {}
    defaults = _pdbzn_workflow_defaults()
    try:
        row_limit = int(f.get("row_limit", defaults["row_limit"]) or defaults["row_limit"])
    except Exception:
        row_limit = defaults["row_limit"]
    if row_limit < 1:
        row_limit = defaults["row_limit"]
    step1_used = {**defaults["step1"], **(f.get("step1", {}) or {})}
    step2_used = {**defaults["step2"], **(f.get("step2", {}) or {})}
    step3_used = {**defaults.get("step3", {}), **(f.get("step3", {}) or {})}
    step4_used = {**defaults.get("step4", {}), **(f.get("step4", {}) or {})}
    require_his3 = bool(step4_used.get("require_his3", True))
    step3_result = _pdbzn_step3_filter_workflow(
        {
            "row_limit": max(row_limit, 1000000),
            "step1": step1_used,
            "step2": step2_used,
            "step3": step3_used,
        }
    )
    if not step3_result.get("ok"):
        return step3_result
    kept_rows = list((((step3_result.get("steps") or [{}])[0].get("rows")) or []))
    if not kept_rows:
        cols = [
            "Cluster_ID",
            "Size",
            "Representative",
            "Members",
            "TMAlign_Threshold",
            "TriHis_Recheck_Pass",
            "Symmetry_Spatial_Filter_Pass",
            "TriHis_Count",
            "TriHis_Residues",
            "ZN_Site",
            "His_ZN_Bond_Lengths_A",
            "His_ZN_Bond_Angles_Deg",
            "Residues_Within_5A_Count",
            "Residues_Within_5A",
            "Validation_Reasons",
        ]
        empty_step = _pdbzn_step_payload("step4_trihis_validate", "Step 4: 综合验证（仅3HIS几何复核）", cols, [], 0, row_limit)
        return {
            "ok": True,
            "db_path": str(PDBZN_DB_PATH),
            "filters_used": {
                "row_limit": row_limit,
                "step1": step1_used,
                "step2": step2_used,
                "step3": step3_used,
                "step4": {"require_his3": require_his3},
            },
            "steps": [empty_step],
            "summary": {
                "step1_count": int((step3_result.get("summary") or {}).get("step1_count", 0) or 0),
                "cluster_count_before_step3": int((step3_result.get("summary") or {}).get("cluster_count_before_step3", 0) or 0),
                "step3_kept_count": 0,
                "step4_cluster_count": 0,
                "step4_valid_count": 0,
                "step4_removed_by_geometry_count": 0,
                "final_count": 0,
                "final_representatives": [],
            },
        }
    out_rows = []
    for c in kept_rows:
        rep = str(c.get("Representative", "") or "").strip().upper()
        geo = _pdbzn_step4_geometry_report(rep, 5.0)
        tri_his_count = int(geo.get("tri_his_count", 0) or 0)
        pass_his = bool(geo.get("ok"))
        if (not require_his3) and tri_his_count >= 1 and bool(geo.get("symmetry_spatial_filter_pass")):
            pass_his = True
        reasons = []
        for rr in (geo.get("reasons") or []):
            reasons.append(str(rr))
        if not pass_his:
            reasons.append("3HIS几何复核失败")
        out_rows.append(
            {
                "Cluster_ID": c.get("Cluster_ID"),
                "Size": c.get("Size"),
                "Representative": rep,
                "Members": c.get("Members"),
                "TMAlign_Threshold": "",
                "TriHis_Recheck_Pass": pass_his,
                "Symmetry_Spatial_Filter_Pass": bool(geo.get("symmetry_spatial_filter_pass")),
                "TriHis_Count": tri_his_count,
                "TriHis_Residues": geo.get("tri_his_residues", ""),
                "ZN_Site": geo.get("zn_site", ""),
                "His_ZN_Bond_Lengths_A": geo.get("his_zn_bond_lengths", ""),
                "His_ZN_Bond_Angles_Deg": geo.get("his_zn_bond_angles", ""),
                "Residues_Within_5A_Count": geo.get("residues_within_5a_count", 0),
                "Residues_Within_5A": geo.get("residues_within_5a", ""),
                "Validation_Reasons": " | ".join(reasons),
            }
        )
    valid_rows = [r for r in out_rows if bool(r.get("TriHis_Recheck_Pass"))]
    cols = [
        "Cluster_ID",
        "Size",
        "Representative",
        "Members",
        "TMAlign_Threshold",
        "TriHis_Recheck_Pass",
        "Symmetry_Spatial_Filter_Pass",
        "TriHis_Count",
        "TriHis_Residues",
        "ZN_Site",
        "His_ZN_Bond_Lengths_A",
        "His_ZN_Bond_Angles_Deg",
        "Residues_Within_5A_Count",
        "Residues_Within_5A",
        "Validation_Reasons",
    ]
    step4 = _pdbzn_step_payload("step4_trihis_validate", "Step 4: 综合验证（仅3HIS几何复核）", cols, valid_rows, len(out_rows), row_limit)
    return {
        "ok": True,
        "db_path": str(PDBZN_DB_PATH),
        "filters_used": {
            "row_limit": row_limit,
            "step1": step1_used,
            "step2": step2_used,
            "step3": step3_used,
            "step4": {"require_his3": require_his3},
        },
        "steps": [step4],
        "summary": {
            "step1_count": int((step3_result.get("summary") or {}).get("step1_count", 0) or 0),
            "cluster_count_before_step3": int((step3_result.get("summary") or {}).get("cluster_count_before_step3", 0) or 0),
            "step3_kept_count": len(kept_rows),
            "step4_cluster_count": len(out_rows),
            "step4_valid_count": len(valid_rows),
            "step4_removed_by_geometry_count": max(0, len(out_rows) - len(valid_rows)),
            "final_count": len(valid_rows),
            "final_representatives": [str(x.get("Representative", "")) for x in valid_rows if str(x.get("Representative", ""))][:200],
        },
    }


def _pdbzn_write_table(path: Path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(list(columns or []))
        for row in (rows or []):
            w.writerow([str((row or {}).get(c, "") or "") for c in columns])


def _pdbzn_write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _pdbzn_read_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _pdbzn_table_rep_or_id(row):
    rep = _pdbzn_rep(row)
    if rep:
        return rep
    rid = str((row or {}).get("_id", "") or "").strip().upper()
    if rid:
        return rid
    receptor = str((row or {}).get("Receptor_PDB", "") or "").strip()
    if receptor:
        return _pdbzn_guess_pdb_id(Path(receptor).name)
    return ""


def _pdbzn_json_table_rows():
    for path in [DATA_DIR / "table1_master_full.json", DATA_DIR / "master_full.json"]:
        hit = _safe_existing_file(path)
        if hit is None:
            continue
        obj = _pdbzn_read_json(hit)
        if not isinstance(obj, dict):
            continue
        rows = list(obj.get("rows") or [])
        header = list(obj.get("header") or [])
        if rows or header:
            return header, rows
    return [], []


def _pdbzn_base_table_rows():
    header, rows = _pdbzn_json_table_rows()
    if rows or header:
        return header, rows
    master_path = _pdbzn_existing_file(["zn_his_master_table7.csv", "zn_his_master_table.csv"])
    if not master_path:
        return [], []
    data = _pdbzn_read_table(master_path)
    header = list(data.get("columns") or [])
    rows = list(data.get("rows") or [])
    return header, rows


def _pdbzn_registry_sort_key(entry):
    table_id = str((entry or {}).get("id", "") or "")
    if table_id == "table1":
        return (0, 0, table_id)
    updated = int((entry or {}).get("updated_at", 0) or 0)
    return (1, -updated, table_id)


def _pdbzn_load_table_registry():
    tables = []
    hit = _safe_existing_file(TABLE_REGISTRY_PATH)
    if hit is not None:
        obj = _pdbzn_read_json(hit)
        if isinstance(obj, dict) and isinstance(obj.get("tables"), list):
            tables = [t for t in obj.get("tables") if isinstance(t, dict)]
    if not any(str((t or {}).get("id", "") or "") == "table1" for t in tables):
        header, rows = _pdbzn_json_table_rows()
        row_count = len(rows)
        tables.append(
            {
                "id": "table1",
                "label": f"Table 1 ({row_count} proteins)" if row_count else "Table 1",
                "full": "table1_master_full.json" if (DATA_DIR / "table1_master_full.json").exists() else "master_full.json",
                "data": "table1_data.json" if (DATA_DIR / "table1_data.json").exists() else "data.json",
                "csv": "table1_master_table.csv" if (DATA_DIR / "table1_master_table.csv").exists() else "master_table.csv",
                "row_count": row_count,
                "description": "主表快照（Table 1）",
            }
        )
    tables.sort(key=_pdbzn_registry_sort_key)
    return tables


def _pdbzn_save_table_registry(tables):
    clean = [t for t in (tables or []) if isinstance(t, dict) and str((t or {}).get("id", "") or "").strip()]
    clean.sort(key=_pdbzn_registry_sort_key)
    _pdbzn_write_json(TABLE_REGISTRY_PATH, {"tables": clean})


def _pdbzn_upsert_registry_entry(entry):
    tables = _pdbzn_load_table_registry()
    table_id = str((entry or {}).get("id", "") or "").strip()
    tables = [t for t in tables if str((t or {}).get("id", "") or "").strip() != table_id]
    tables.append(entry)
    _pdbzn_save_table_registry(tables)


def _pdbzn_label_slug(label):
    raw = safe_name(str(label or "").strip().lower()).strip("._-")
    return raw or "workflow_table"


def _pdbzn_unique_table_identity(label):
    base_label = str(label or "").strip() or time.strftime("Step5 %Y-%m-%d %H:%M:%S", time.localtime())
    base_slug = _pdbzn_label_slug(base_label)
    tables = _pdbzn_load_table_registry()
    used_ids = {str((t or {}).get("id", "") or "").strip().lower() for t in tables}
    used_labels = {str((t or {}).get("label", "") or "").strip() for t in tables}
    if base_slug not in used_ids and base_label not in used_labels:
        return base_slug, base_label
    idx = 2
    while True:
        cand_id = f"{base_slug}_{idx}"
        cand_label = f"{base_label} ({idx})"
        if cand_id not in used_ids and cand_label not in used_labels:
            return cand_id, cand_label
        idx += 1


def _pdbzn_receptor_rel(rep, receptor_path):
    rid = str(rep or "").strip().upper()
    receptor_file = _pdbzn_resolve_receptor_pdb_path(rid, receptor_path)
    local = _safe_existing_file(DATA_DIR / "structures" / f"{rid}.pdb")
    if local is None:
        local = _safe_existing_file(DATA_DIR / "structures" / f"{rid}.cif")
    if local is not None:
        return f"../backend/data/structures/{local.name}"
    if receptor_file is not None:
        try:
            rel = receptor_file.resolve().relative_to(DATA_DIR.resolve())
            return "../backend/data/" + str(rel).replace("\\", "/")
        except Exception:
            return ""
    return ""


def _pdbzn_dataset_item_from_row(row):
    rep = _pdbzn_table_rep_or_id(row)
    receptor_path = str((row or {}).get("Receptor_PDB", "") or "").strip()
    length, oligomer, monomer = _pdbzn_parse_len_oligomer_monomer((row or {}).get("Length;Oligomer;Monomer"))
    return {
        "id": rep,
        "name": str((row or {}).get("Protein_Name", "") or ""),
        "cluster": str((row or {}).get("Cluster_ID", "") or (row or {}).get("Step5_Cluster_ID", "") or ""),
        "receptor_pdb": receptor_path,
        "receptor_rel": _pdbzn_receptor_rel(rep, receptor_path),
        "best_sdf": str((row or {}).get("Best_SDF", "") or ""),
        "species": str((row or {}).get("Species", "") or ""),
        "monomer_seq": str((row or {}).get("MonomerSeq", "") or ""),
        "residue_length": "" if length is None else int(length),
        "oligomer": "" if oligomer is None else int(oligomer),
        "monomer_length": "" if monomer is None else int(monomer),
        "zn_coord_residue_count": str((row or {}).get("Zn_CoordResidueCount", "") or (row or {}).get("Residues_Within_5A_Count", "") or ""),
        "zn_coord_his_count": str((row or {}).get("Zn_CoordHisCount", "") or (row or {}).get("TriHisCountMax", "") or (row or {}).get("Step5_TriHis_Count", "") or ""),
        "zn_coord_non_his_count": str((row or {}).get("Zn_CoordNonHisCount", "") or ""),
    }


def _pdbzn_write_named_table(label, header, rows, stage, description):
    with _pdbzn_table_registry_lock:
        safe_id, unique_label = _pdbzn_unique_table_identity(label)
        clean_rows = []
        ids = []
        for src in rows or []:
            row = {c: (src or {}).get(c, "") for c in header}
            clean_rows.append(row)
            rep = _pdbzn_table_rep_or_id(row)
            if rep and rep not in ids:
                ids.append(rep)
        full_name = f"{safe_id}_master_full.json"
        data_name = f"{safe_id}_data.json"
        csv_name = f"{safe_id}_master_table.csv"
        ids_name = f"{safe_id}_ids.json"
        data_items = [_pdbzn_dataset_item_from_row(r) for r in clean_rows if _pdbzn_table_rep_or_id(r)]
        _pdbzn_write_json(DATA_DIR / full_name, {"rows": clean_rows, "header": list(header or [])})
        _pdbzn_write_json(DATA_DIR / data_name, {"items": data_items})
        _pdbzn_write_json(DATA_DIR / ids_name, {"table_id": safe_id, "label": unique_label, "ids": ids})
        _pdbzn_write_table(DATA_DIR / csv_name, ["_id"] + list(header or []), [{"_id": _pdbzn_table_rep_or_id(r), **r} for r in clean_rows])
        entry = {
            "id": safe_id,
            "label": unique_label,
            "full": full_name,
            "data": data_name,
            "csv": csv_name,
            "ids": ids_name,
            "row_count": len(clean_rows),
            "stage": stage,
            "description": description,
            "updated_at": now_ts(),
            "updated_at_text": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        }
        _pdbzn_upsert_registry_entry(entry)
        return entry


def _pdbzn_table_by_rep(path: Path):
    out = {}
    if not path or not path.exists():
        return out
    data = _pdbzn_read_table(path)
    for r in data.get("rows", []):
        rep = _pdbzn_rep(r)
        if rep:
            out[rep] = r
    return out


def _pdbzn_guess_diffdock_for_rep(rep):
    pid = str(rep or "").strip().upper()
    if not pid:
        return {}
    bases = [
        DATA_DIR / "diffdock" / pid,
        PDBZN_BASE_DIR / "diffdock_outputs" / pid / pid,
        PDBZN_BASE_DIR / "diffdock" / pid,
    ]
    base = None
    for cand in bases:
        hit = _safe_existing_dir(cand)
        if hit is not None:
            base = hit
            break
    if base is None:
        return {}
    best_file = None
    best_conf = None
    best_rank = None
    for sdf in sorted(base.glob("rank*_confidence-*.sdf")):
        m_rank = re.search(r"rank(\d+)", sdf.name, flags=re.IGNORECASE)
        m_conf = re.search(r"confidence-([-+]?[0-9]*\.?[0-9]+)", sdf.name, flags=re.IGNORECASE)
        rank_num = int(m_rank.group(1)) if m_rank else None
        conf_num = _pdbzn_to_num(m_conf.group(1)) if m_conf else None
        if best_file is None:
            best_file = sdf
            best_conf = conf_num
            best_rank = rank_num
            continue
        if conf_num is not None and (best_conf is None or float(conf_num) > float(best_conf)):
            best_file = sdf
            best_conf = conf_num
            best_rank = rank_num
            continue
        if conf_num is None and best_conf is None and rank_num is not None and (best_rank is None or rank_num < best_rank):
            best_file = sdf
            best_rank = rank_num
    if best_file is None:
        return {}
    rec = _pdbzn_find_tmalign_pdb_path(pid)
    out = {"Best_SDF": str(best_file), "Best_Confidence": best_conf, "Best_Rank": best_rank}
    if rec is not None:
        out["Receptor_PDB"] = str(rec)
    return out


def _pdbzn_resolve_best_sdf_path(rep, best_sdf_path):
    direct = _safe_existing_file(best_sdf_path)
    if direct is not None:
        return direct
    pid = str(rep or "").strip().upper()
    if not pid:
        return None
    name = Path(str(best_sdf_path or "").strip()).name
    if not name:
        return None
    candidates = [
        DATA_DIR / "diffdock" / pid / name,
        PDBZN_BASE_DIR / "diffdock_outputs" / pid / pid / name,
        PDBZN_BASE_DIR / "diffdock" / pid / name,
    ]
    for cand in candidates:
        hit = _safe_existing_file(cand)
        if hit is not None:
            return hit
    return None


def _pdbzn_parse_sdf_atoms_text(text):
    lines = str(text or "").splitlines()
    if len(lines) < 4:
        return []
    try:
        n_atoms = int(str(lines[3] or "")[:3].strip())
    except Exception:
        return []
    if n_atoms <= 0:
        return []
    atoms = []
    for i in range(4, min(4 + n_atoms, len(lines))):
        line = str(lines[i] or "")
        try:
            x = float(line[0:10].strip())
            y = float(line[10:20].strip())
            z = float(line[20:30].strip())
        except Exception:
            continue
        elem = str(line[31:34] or "").strip().upper()
        atoms.append({"x": x, "y": y, "z": z, "elem": elem})
    return atoms


def _pdbzn_best_his_zn_site(atoms):
    if not atoms:
        return {"zn": None, "zn_atoms": [], "tri_his_zn_atoms": []}
    zn_atoms = [
        a for a in atoms
        if str(a.get("comp") or "").upper() == "ZN" or str(a.get("element") or "").upper() == "ZN"
    ]
    if not zn_atoms:
        return {"zn": None, "zn_atoms": [], "tri_his_zn_atoms": []}
    his_atoms = [
        a for a in atoms
        if str(a.get("comp") or "").upper() in {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}
        and str(a.get("atom") or "").upper() in {"ND1", "NE2"}
    ]
    best_any = None
    best_3 = None
    tri_his_zn_atoms = []
    for z in zn_atoms:
        nearest_by_res = {}
        for ha in his_atoms:
            d = _pdbzn_dist(z, ha)
            if d > 3.0:
                continue
            rk = _pdbzn_residue_key(ha)
            old = nearest_by_res.get(rk)
            if old is None or d < old["distance"]:
                nearest_by_res[rk] = {"atom": ha, "distance": float(d)}
        picked = sorted(nearest_by_res.items(), key=lambda kv: kv[1]["distance"])
        tri = picked[:3]
        candidate = {
            "zn": z,
            "his_count": len(picked),
            "score_sum3": sum(x[1]["distance"] for x in tri) if len(tri) >= 3 else 9999.0,
        }
        if best_any is None or candidate["his_count"] > best_any["his_count"] or (
            candidate["his_count"] == best_any["his_count"] and candidate["score_sum3"] < best_any["score_sum3"]
        ):
            best_any = candidate
        if candidate["his_count"] >= 3:
            tri_his_zn_atoms.append(z)
            if best_3 is None or candidate["his_count"] > best_3["his_count"] or (
                candidate["his_count"] == best_3["his_count"] and candidate["score_sum3"] < best_3["score_sum3"]
            ):
                best_3 = candidate
    best = best_3 or best_any
    return {"zn": (best or {}).get("zn"), "zn_atoms": zn_atoms, "tri_his_zn_atoms": tri_his_zn_atoms}


def _pdbzn_compute_diffdock_zn_metrics(rep, best_sdf_path, receptor_file):
    sdf_hit = _pdbzn_resolve_best_sdf_path(rep, best_sdf_path)
    receptor_hit = _safe_existing_file(receptor_file)
    if sdf_hit is None or receptor_hit is None:
        return {}
    try:
        lig_atoms = [a for a in _pdbzn_parse_sdf_atoms_text(sdf_hit.read_text(encoding="utf-8", errors="ignore")) if str(a.get("elem") or "").upper() != "H"]
    except Exception:
        return {}
    if not lig_atoms:
        return {}
    try:
        receptor_atoms = _pdbzn_structure_atom_rows(receptor_hit)
    except Exception:
        return {}
    best_site = _pdbzn_best_his_zn_site(receptor_atoms)
    target_zn = best_site.get("zn")
    zn_atoms = list(best_site.get("tri_his_zn_atoms") or [])
    if not zn_atoms:
        zn_atoms = [target_zn] if target_zn is not None else list(best_site.get("zn_atoms") or [])
    if target_zn is None or not zn_atoms:
        return {}
    nearest = None
    for a in lig_atoms:
        d = _pdbzn_dist(a, target_zn)
        if nearest is None or d < nearest:
            nearest = float(d)
    all_zn = []
    for z in zn_atoms:
        zmin = None
        for a in lig_atoms:
            d = _pdbzn_dist(a, z)
            if zmin is None or d < zmin:
                zmin = float(d)
        if zmin is not None:
            all_zn.append(zmin)
    all_zn.sort()
    out = {}
    if nearest is not None:
        out["NearestDistanceTo3HisZn"] = round(nearest, 12)
    if all_zn:
        out["AllZnDistancesSorted"] = ";".join(f"{x:.3f}" for x in all_zn)
    return out


def _pdbzn_fpocket_num(params, key_contains):
    if not isinstance(params, dict):
        return None
    for k, v in params.items():
        key = str(k or "").strip().lower()
        if key_contains in key:
            n = _pdbzn_to_num(v)
            if n is not None:
                return float(n)
    return None


def _pdbzn_fpocket_param_num(params, keywords):
    if not isinstance(params, dict):
        return None
    wanted = [str(x or "").strip().lower() for x in (keywords or []) if str(x or "").strip()]
    if not wanted:
        return None
    for k, v in params.items():
        key = str(k or "").strip().lower()
        if all(token in key for token in wanted):
            n = _pdbzn_to_num(v)
            if n is not None:
                return float(n)
    return None


def _pdbzn_fpocket_best_from_rows(rows):
    best = None
    for row in rows or []:
        params = row.get("params", {}) if isinstance(row, dict) else {}
        drug = _pdbzn_fpocket_param_num(params, ["drug", "score"])
        pocket_score = _pdbzn_fpocket_param_num(params, ["pocket", "score"])
        vol = _pdbzn_fpocket_param_num(params, ["volume"])
        alpha = _pdbzn_fpocket_param_num(params, ["number", "alpha", "sphere"])
        item = {
            "pocket_id": row.get("pocket_id"),
            "druggability": drug,
            "score": pocket_score,
            "volume": vol,
            "alpha_spheres": alpha,
        }
        if best is None:
            best = item
            continue
        best_drug = best.get("druggability")
        best_score = best.get("score")
        best_volume = best.get("volume")
        if drug is not None and (best_drug is None or drug > best_drug):
            best = item
            continue
        if drug is None and best_drug is None and pocket_score is not None and (best_score is None or pocket_score > best_score):
            best = item
            continue
        if drug is None and best_drug is None and pocket_score is None and best_score is None and vol is not None and (best_volume is None or vol > best_volume):
            best = item
    return best


def _pdbzn_num_or(value, default):
    n = _pdbzn_to_num(value)
    return default if n is None else n


def _pdbzn_load_pocket_rows_by_rep():
    global _pdbzn_pocket_table_cache
    with _pdbzn_pocket_table_cache_lock:
        if _pdbzn_pocket_table_cache is not None:
            return _pdbzn_pocket_table_cache
        path = _safe_existing_file(DATA_DIR / "pocket.csv")
        grouped = {}
        if path is not None:
            data = _pdbzn_read_table(path)
            for row in data.get("rows", []):
                rep = str((row or {}).get("Representative", "") or "").strip().upper()
                if rep:
                    grouped.setdefault(rep, []).append(row)
        _pdbzn_pocket_table_cache = grouped
        return grouped


def _pdbzn_best_pocket_catalog_row(rep):
    grouped = _pdbzn_load_pocket_rows_by_rep()
    rows = list(grouped.get(str(rep or "").strip().upper(), []) or [])
    if not rows:
        return None
    matching = [r for r in rows if _pdbzn_is_truthy((r or {}).get("zn_match", ""))]
    if matching:
        best = sorted(
            matching,
            key=lambda r: (
                float(_pdbzn_num_or((r or {}).get("min_dist_to_any_zn", ""), 10**9)),
                -float(_pdbzn_num_or((r or {}).get("druggability_score", ""), -10**9)),
                -float(_pdbzn_num_or((r or {}).get("score", ""), -10**9)),
                -float(_pdbzn_num_or((r or {}).get("volume", ""), -10**9)),
                int(_pdbzn_num_or((r or {}).get("pocket_id", ""), 10**9)),
            ),
        )[0]
        rule = "zn_priority"
    else:
        best = sorted(
            rows,
            key=lambda r: (
                -float(_pdbzn_num_or((r or {}).get("volume", ""), -10**9)),
                -float(_pdbzn_num_or((r or {}).get("druggability_score", ""), -10**9)),
                -float(_pdbzn_num_or((r or {}).get("score", ""), -10**9)),
                int(_pdbzn_num_or((r or {}).get("pocket_id", ""), 10**9)),
            ),
        )[0]
        rule = "largest_volume"
    return {
        "Representative": str(rep or "").strip().upper(),
        "Protein_Name": best.get("Protein_Name", ""),
        "Protein_Category": best.get("Protein_Category", ""),
        "Species": best.get("Species", ""),
        "Receptor_PDB": best.get("Receptor_PDB", ""),
        "BestPocket_ID": best.get("pocket_id", ""),
        "BestPocket_Score": best.get("score", ""),
        "BestPocket_Druggability": best.get("druggability_score", ""),
        "BestPocket_Volume": best.get("volume", ""),
        "BestPocket_TotalSASA": best.get("total_sasa", ""),
        "BestPocket_PolarSASA": best.get("polar_sasa", ""),
        "BestPocket_ApolarSASA": best.get("apolar_sasa", ""),
        "BestPocket_AlphaSpheres": best.get("alpha_spheres", ""),
        "BestPocket_HisCount": best.get("pocket_his_count", ""),
        "BestPocket_ZnCount": best.get("pocket_zn_count", ""),
        "BestPocket_MinDistToZn": best.get("min_dist_to_any_zn", ""),
        "BestPocket_ZnMatch": best.get("zn_match", ""),
        "BestPocket_SelectRule": rule,
    }


def _pdbzn_geometry_summary_from_step4(row):
    if not isinstance(row, dict):
        return ""
    lengths = []
    for part in str(row.get("His_ZN_Bond_Lengths_A", "") or "").split(";"):
        seg = str(part or "").strip()
        if not seg:
            continue
        tail = seg.rsplit(":", 1)[-1]
        val = _pdbzn_to_num(tail)
        if val is not None:
            lengths.append(float(val))
    angles = []
    for part in str(row.get("His_ZN_Bond_Angles_Deg", "") or "").split(";"):
        seg = str(part or "").strip()
        if not seg:
            continue
        tail = seg.rsplit(":", 1)[-1]
        val = _pdbzn_to_num(tail)
        if val is not None:
            angles.append(float(val))
    if not lengths and not angles:
        return ""
    pieces = []
    if lengths:
        pieces.append("D=[" + ",".join(f"{x:.2f}" for x in lengths) + "]")
    if angles:
        pieces.append("A=[" + ",".join(f"{x:.1f}" for x in angles) + "]")
    return ",".join(pieces)


def _pdbzn_run_fpocket_quick(pdb_path: Path, timeout_sec: float):
    fpocket_bin = discover_fpocket_bin()
    if not fpocket_bin:
        return {"status": "fpocket_not_found"}
    pdb_hit = _safe_existing_file(pdb_path)
    if pdb_hit is None:
        return {"status": "pdb_not_found"}
    pdb_path = pdb_hit
    work_root = FPOCKET_WORK_DIR / f"step5_{uuid.uuid4().hex[:10]}"
    work_root.mkdir(parents=True, exist_ok=True)
    try:
        local_pdb = work_root / pdb_path.name
        shutil.copy2(str(pdb_path), str(local_pdb))
        proc = subprocess.run(
            [fpocket_bin, "-f", str(local_pdb)],
            cwd=str(work_root),
            capture_output=True,
            text=True,
            timeout=max(10.0, float(timeout_sec or 120.0)),
            env=os.environ.copy(),
        )
        if proc.returncode != 0:
            return {"status": "fpocket_failed"}
        out_dir = work_root / f"{local_pdb.stem}_out"
        if not out_dir.exists() or not out_dir.is_dir():
            return {"status": "fpocket_no_output"}
        info_files = sorted(out_dir.glob("*_info.txt"))
        if not info_files:
            return {"status": "fpocket_no_info"}
        info = parse_fpocket_info_file(info_files[0])
        rows = list(info.get("rows", []))
        if not rows:
            return {"status": "fpocket_empty"}
        best = _pdbzn_fpocket_best_from_rows(rows)
        if not best:
            return {"status": "fpocket_empty"}
        return {
            "status": "ok",
            "best_pocket_id": best.get("pocket_id"),
            "best_score": best.get("score"),
            "best_druggability": best.get("druggability"),
            "best_volume": best.get("volume"),
            "best_alpha_spheres": best.get("alpha_spheres"),
        }
    except subprocess.TimeoutExpired:
        return {"status": "fpocket_timeout"}
    except Exception:
        return {"status": "fpocket_error"}
    finally:
        try:
            shutil.rmtree(work_root, ignore_errors=True)
        except Exception:
            pass


def _pdbzn_resolve_receptor_pdb_path(rep, receptor_path):
    direct = _safe_existing_file(receptor_path)
    if direct is not None:
        return direct
    return _pdbzn_find_tmalign_pdb_path(rep)


def _pdbzn_step5_final_score(row):
    conf = _pdbzn_to_num((row or {}).get("Best_Confidence"))
    near = _pdbzn_to_num((row or {}).get("NearestDistanceTo3HisZn"))
    his = _pdbzn_to_num((row or {}).get("TriHisCountMax"))
    fp = _pdbzn_to_num((row or {}).get("BestPocket_Druggability"))
    if fp is None:
        fp = _pdbzn_to_num((row or {}).get("BestPocket_Score"))
    if fp is None:
        fp = _pdbzn_to_num((row or {}).get("Step5_FPocket_BestScore"))
    score = 0.0
    if conf is not None:
        score += float(conf)
    if near is not None:
        score += max(0.0, 10.0 - float(near))
    if his is not None:
        score += min(3.0, max(0.0, float(his)))
    if fp is not None:
        score += float(fp)
    return round(score, 4)


def _pdbzn_step5_finalize_workflow(filters):
    f = filters or {}
    defaults = _pdbzn_workflow_defaults()
    try:
        row_limit = int(f.get("row_limit", defaults["row_limit"]) or defaults["row_limit"])
    except Exception:
        row_limit = defaults["row_limit"]
    if row_limit < 1:
        row_limit = defaults["row_limit"]
    step1_used = {**defaults["step1"], **(f.get("step1", {}) or {})}
    step2_used = {**defaults["step2"], **(f.get("step2", {}) or {})}
    step3_used = {**defaults["step3"], **(f.get("step3", {}) or {})}
    step4_used = {**defaults["step4"], **(f.get("step4", {}) or {})}
    step5_used = {**defaults.get("step5", {}), **(f.get("step5", {}) or {})}
    run_fpocket = bool(step5_used.get("run_fpocket", True))
    fpocket_max_runs = _pdbzn_to_num(step5_used.get("fpocket_max_runs", 20))
    if fpocket_max_runs is None:
        fpocket_max_runs = 20
    fpocket_max_runs = int(max(0, fpocket_max_runs))
    fpocket_timeout_sec = _pdbzn_to_num(step5_used.get("fpocket_timeout_sec", 120))
    if fpocket_timeout_sec is None:
        fpocket_timeout_sec = 120.0
    fpocket_timeout_sec = float(max(10.0, fpocket_timeout_sec))
    write_to_master = bool(step5_used.get("write_to_master", True))
    table_label = str(step5_used.get("table_label", "") or "").strip() or time.strftime("Step5 %Y-%m-%d %H:%M:%S", time.localtime())
    output_name = str(step5_used.get("output_file") or "zn_his_master_table_step5.csv").strip() or "zn_his_master_table_step5.csv"
    output_path = PDBZN_EXPORT_DIR / output_name
    step4_result = _pdbzn_step4_validate_workflow(
        {
            "row_limit": max(row_limit, 1000000),
            "step1": step1_used,
            "step2": step2_used,
            "step3": step3_used,
            "step4": step4_used,
        }
    )
    if not step4_result.get("ok"):
        return step4_result
    step4_rows = list((((step4_result.get("steps") or [{}])[0].get("rows")) or []))
    master_path = _pdbzn_existing_file(["zn_his_master_table7.csv", "zn_his_master_table.csv"])
    base_header, base_master_rows = _pdbzn_base_table_rows()
    master_cols = [c for c in (base_header or []) if c and c != "_id" and not str(c).startswith("Step5_")]
    if not master_cols and base_master_rows:
        master_cols = [c for c in base_master_rows[0].keys() if c and c != "_id" and not str(c).startswith("Step5_")]
    required_cols = [
        "Representative",
        "Similarity_Score",
        "Best_Confidence",
        "NearestDistanceTo3HisZn",
        "ZN_Depth",
        "ZN_SASA",
        "Length;Oligomer;Monomer",
        "Protein_Name",
        "Protein_Category",
        "Cluster_ID",
        "Details",
        "ZN_Depth_Rounded",
        "Neighbor_Count",
        "Neighbor_Profile_Similarity",
        "Neighbor_List",
        "Best_Rank",
        "AllZnDistancesSorted",
        "TriHisSatisfied",
        "TriHisCountMax",
        "Receptor_PDB",
        "Best_SDF",
        "Species",
        "MonomerSeq",
        "BestPocket_ID",
        "BestPocket_Score",
        "BestPocket_Druggability",
        "BestPocket_Volume",
        "BestPocket_TotalSASA",
        "BestPocket_PolarSASA",
        "BestPocket_ApolarSASA",
        "BestPocket_AlphaSpheres",
        "BestPocket_HisCount",
        "BestPocket_ZnCount",
        "BestPocket_MinDistToZn",
        "BestPocket_ZnMatch",
        "BestPocket_SelectRule",
        "Zn_CoordResidueCount",
        "Zn_CoordHisCount",
        "Zn_CoordNonHisCount",
        "Zn_CoordResidues",
        "Zn_CoordSite",
        "Zn_CoordNote",
    ]
    for c in required_cols:
        if c not in master_cols:
            master_cols.append(c)
    base_rows = []
    for r in base_master_rows:
        source_tag = str((r or {}).get("Step5_Source", "") or "").strip()
        if source_tag in {"step4_finalize", "step4_validate", "step5_finalize"}:
            continue
        base_rows.append({c: (r or {}).get(c, "") for c in master_cols})
    master_by_rep = {}
    for r in base_master_rows:
        rep = _pdbzn_table_rep_or_id(r)
        if rep:
            master_by_rep[rep] = r
    diffdock_path = _pdbzn_existing_file(["diffdock_screen5.csv", "diffdock_screen.csv"])
    diffdock_by_rep = _pdbzn_table_by_rep(diffdock_path) if diffdock_path else {}
    fpocket_run_count = 0
    fpocket_status_count = {}
    appended_rows = []
    step_result_rows = []
    for s4 in step4_rows:
        rep = str((s4 or {}).get("Representative", "") or "").strip().upper()
        if not rep:
            continue
        row = {c: "" for c in master_cols}
        src = master_by_rep.get(rep, {})
        for c in master_cols:
            v = src.get(c, "")
            if v is not None and str(v) != "":
                row[c] = v
        row["Representative"] = rep
        pocket_catalog = _pdbzn_best_pocket_catalog_row(rep) or {}
        for k, v in pocket_catalog.items():
            if k in row and str(row.get(k, "")).strip() == "" and str(v or "").strip() != "":
                row[k] = v
        if str(s4.get("Cluster_ID", "")).strip() != "":
            row["Cluster_ID"] = s4.get("Cluster_ID", "")
        row["TriHisSatisfied"] = "true" if _pdbzn_is_truthy(s4.get("TriHis_Recheck_Pass")) else "false"
        s4_his = _pdbzn_to_num(s4.get("TriHis_Count"))
        if s4_his is not None:
            row["TriHisCountMax"] = int(s4_his)
        else:
            tri_his = _pdbzn_to_num(row.get("TriHisCountMax"))
            if tri_his is None:
                reasons, stats = _pdbzn_step3_structure_reasons(rep, step3_used.get("ligand_distance", 5.0), step3_used.get("metal_distance", 8.0))
                max_his = stats.get("max_his_coord")
                if max_his is not None:
                    row["TriHisCountMax"] = int(max_his)
        details_summary = _pdbzn_geometry_summary_from_step4(s4)
        if details_summary and str(row.get("Details", "")).strip() in {"", "STEP5_FROM_STEP4"}:
            row["Details"] = details_summary
        drow = diffdock_by_rep.get(rep, {})
        if not drow:
            drow = _pdbzn_guess_diffdock_for_rep(rep)
        for k in ["Receptor_PDB", "Best_SDF", "Best_Rank", "Best_Confidence", "NearestDistanceTo3HisZn", "AllZnDistancesSorted", "TriHisSatisfied", "TriHisCountMax"]:
            if str(row.get(k, "")).strip() == "" and str(drow.get(k, "")).strip() != "":
                row[k] = drow.get(k)
        receptor_path = str(row.get("Receptor_PDB", "") or "").strip()
        receptor_file = _pdbzn_resolve_receptor_pdb_path(rep, receptor_path)
        if receptor_file is not None:
            row["Receptor_PDB"] = str(receptor_file)
        diffdock_geom = _pdbzn_compute_diffdock_zn_metrics(rep, row.get("Best_SDF", ""), receptor_file)
        for k, v in diffdock_geom.items():
            if v is not None and str(v).strip() != "":
                row[k] = v
        coord = _pdbzn_coordination_metrics(rep, receptor_file)
        if coord:
            for k, v in coord.items():
                row[k] = v
        fp_status = "not_run"
        pocket_minimal_missing = any(
            str(row.get(k, "")).strip() == ""
            for k in ["BestPocket_ID", "BestPocket_Score", "BestPocket_Druggability", "BestPocket_Volume"]
        )
        if run_fpocket and receptor_file and fpocket_run_count < fpocket_max_runs and pocket_minimal_missing:
            fp_res = _pdbzn_run_fpocket_quick(receptor_file, fpocket_timeout_sec)
            fp_status = str(fp_res.get("status", "fpocket_error"))
            if fp_status == "ok":
                best_pocket_id = fp_res.get("best_pocket_id", "")
                best_score = fp_res.get("best_score")
                best_druggability = fp_res.get("best_druggability")
                best_volume = fp_res.get("best_volume")
                best_alpha_spheres = fp_res.get("best_alpha_spheres")
                if str(row.get("BestPocket_ID", "")).strip() == "" and str(best_pocket_id or "").strip() != "":
                    row["BestPocket_ID"] = best_pocket_id
                if str(row.get("BestPocket_Score", "")).strip() == "" and best_score is not None:
                    row["BestPocket_Score"] = best_score
                if str(row.get("BestPocket_Druggability", "")).strip() == "" and best_druggability is not None:
                    row["BestPocket_Druggability"] = best_druggability
                if str(row.get("BestPocket_Volume", "")).strip() == "" and best_volume is not None:
                    row["BestPocket_Volume"] = best_volume
                if str(row.get("BestPocket_AlphaSpheres", "")).strip() == "" and best_alpha_spheres is not None:
                    row["BestPocket_AlphaSpheres"] = best_alpha_spheres
            fpocket_run_count += 1
        elif run_fpocket and not receptor_file:
            fp_status = "pdb_unavailable"
        fpocket_status_count[fp_status] = int(fpocket_status_count.get(fp_status, 0) or 0) + 1
        final_score = _pdbzn_step5_final_score(row)
        appended_rows.append(row)
        step_result_rows.append(
            {
                "Representative": rep,
                "Cluster_ID": s4.get("Cluster_ID", ""),
                "Cluster_Size": s4.get("Size", ""),
                "TriHis_Recheck_Pass": s4.get("TriHis_Recheck_Pass", ""),
                "Symmetry_Spatial_Filter_Pass": s4.get("Symmetry_Spatial_Filter_Pass", ""),
                "TriHis_Count": row.get("TriHisCountMax", ""),
                "Zn_CoordSite": row.get("Zn_CoordSite", "") or s4.get("ZN_Site", ""),
                "Zn_CoordResidues": row.get("Zn_CoordResidues", ""),
                "Best_Confidence": row.get("Best_Confidence", ""),
                "NearestDistanceTo3HisZn": row.get("NearestDistanceTo3HisZn", ""),
                "BestPocket_ID": row.get("BestPocket_ID", ""),
                "BestPocket_Druggability": row.get("BestPocket_Druggability", ""),
                "BestPocket_Score": row.get("BestPocket_Score", ""),
                "BestPocket_Volume": row.get("BestPocket_Volume", ""),
                "FPocket_Status": fp_status,
                "FinalScore": final_score,
                "Receptor_PDB": row.get("Receptor_PDB", ""),
                "Best_SDF": row.get("Best_SDF", ""),
            }
        )
    generated_table = _pdbzn_write_named_table(
        table_label,
        master_cols,
        appended_rows,
        "step5_finalize",
        "PDB 工作流生成结果（主表标准列）",
    )
    _pdbzn_write_table(output_path, master_cols, appended_rows)
    written_master = False
    if write_to_master and master_path:
        out_rows = list(base_rows)
        if appended_rows:
            out_rows.append({c: "" for c in master_cols})
            out_rows.extend(appended_rows)
        _pdbzn_write_table(master_path, master_cols, out_rows)
        written_master = True
    step_cols = [
        "Representative",
        "Cluster_ID",
        "Cluster_Size",
        "TriHis_Recheck_Pass",
        "Symmetry_Spatial_Filter_Pass",
        "TriHis_Count",
        "Zn_CoordSite",
        "Zn_CoordResidues",
        "Best_Confidence",
        "NearestDistanceTo3HisZn",
        "BestPocket_ID",
        "BestPocket_Druggability",
        "BestPocket_Score",
        "BestPocket_Volume",
        "FPocket_Status",
        "FinalScore",
        "Receptor_PDB",
        "Best_SDF",
    ]
    step5 = _pdbzn_step_payload("step5_master_finalize", "Step 5: 生成主表列并补全 DiffDock / fpocket / Zn 配位信息", step_cols, step_result_rows, len(step4_rows), row_limit)
    return {
        "ok": True,
        "db_path": str(PDBZN_DB_PATH),
        "filters_used": {
            "row_limit": row_limit,
            "step1": step1_used,
            "step2": step2_used,
            "step3": step3_used,
            "step4": step4_used,
            "step5": {
                "run_fpocket": run_fpocket,
                "fpocket_max_runs": fpocket_max_runs,
                "fpocket_timeout_sec": fpocket_timeout_sec,
                "write_to_master": write_to_master,
                "table_label": table_label,
                "output_file": output_name,
            },
        },
        "steps": [step5],
        "summary": {
            "step4_valid_count": len(step4_rows),
            "appended_count": len(appended_rows),
            "step5_trihis_pass_count": len([x for x in step4_rows if _pdbzn_is_truthy((x or {}).get("TriHis_Recheck_Pass", ""))]),
            "step5_symmetry_pass_count": len([x for x in step4_rows if _pdbzn_is_truthy((x or {}).get("Symmetry_Spatial_Filter_Pass", ""))]),
            "fpocket_runs": fpocket_run_count,
            "fpocket_status_count": fpocket_status_count,
            "final_count": len(appended_rows),
            "final_representatives": [str(x.get("Representative", "")) for x in appended_rows if str(x.get("Representative", ""))][:200],
            "master_file": str(master_path) if master_path else "",
            "output_file": str(output_path),
            "written_master": written_master,
            "table_label": generated_table.get("label", table_label),
            "generated_table": generated_table,
            "table_registry_file": str(TABLE_REGISTRY_PATH),
        },
    }


def resolve_under(base: Path, rel_path: str):
    cleaned = (rel_path or "").replace("\\", "/").lstrip("/")
    target = (base / cleaned).resolve()
    base_real = base.resolve()
    if not str(target).startswith(str(base_real)):
        return None
    return target


def recover_fpocket_job(job_id: str):
    safe_job_id = safe_name(job_id)
    if safe_job_id != job_id:
        return None
    job_dir = FPOCKET_WORK_DIR / job_id
    if (not job_dir.exists()) or (not job_dir.is_dir()):
        return None
    in_dir = job_dir / "input"
    pdb_name = ""
    if in_dir.exists():
        pdb_files = sorted([p for p in in_dir.iterdir() if p.is_file() and p.suffix.lower() in {".pdb", ".ent"}], key=lambda p: p.name)
        if pdb_files:
            pdb_name = pdb_files[0].name
    outputs = list_fpocket_outputs(job_id)
    log_path = job_dir / "run.log"
    fpocket_bin = ""
    if log_path.exists():
        try:
            first_line = next((ln for ln in log_path.read_text(encoding="utf-8", errors="ignore").splitlines() if ln.startswith("CMD: ")), "")
            if first_line:
                parts = shlex.split(first_line[5:].strip())
                if parts:
                    fpocket_bin = parts[0]
        except Exception:
            fpocket_bin = ""
    created_at = int(job_dir.stat().st_mtime)
    status = "done" if outputs else ("failed" if log_path.exists() else "queued")
    error = "" if outputs else ("fpocket输出不存在，可能执行失败或后端重启后状态丢失" if log_path.exists() else "")
    exit_code = 0 if outputs else None
    job = {
        "id": job_id,
        "status": status,
        "created_at": created_at,
        "started_at": created_at,
        "ended_at": created_at,
        "pdb_name": pdb_name,
        "outputs": outputs,
        "error": error,
        "pid": None,
        "fpocket_bin": fpocket_bin,
        "phase": "done_recovered" if outputs else "recovered",
        "phase_at": now_ts(),
        "exit_code": exit_code,
    }
    with fpocket_jobs_lock:
        fpocket_jobs[job_id] = job
    return job


def run_job(job_id: str):
    with jobs_lock:
        job = jobs.get(job_id)
        if not job:
            return
        job["status"] = "running"
        job["started_at"] = now_ts()
        job["phase"] = "starting"
        job["phase_at"] = now_ts()
    job_dir = WORK_DIR / job_id
    in_dir = job_dir / "input"
    out_dir = job_dir / "output"
    log_path = job_dir / "run.log"
    with jobs_lock:
        job["phase"] = "prepare_dirs"
        job["phase_at"] = now_ts()
    out_dir.mkdir(parents=True, exist_ok=True)
    receptor = in_dir / job["receptor_name"]
    ligand = in_dir / job["ligand_name"]
    inference_steps = str(job.get("inference_steps", 20))
    samples_per_complex = str(job.get("samples_per_complex", 10))
    with jobs_lock:
        job["phase"] = "open_log"
        job["phase_at"] = now_ts()
    with open(log_path, "a", encoding="utf-8") as logf:
        if not INFER_PY:
            msg = "未找到DiffDock inference.py，请设置环境变量 DIFFDOCK_INFER_PY"
            logf.write(msg + "\n")
            with jobs_lock:
                job["status"] = "failed"
                job["error"] = msg
                job["ended_at"] = now_ts()
                job["phase"] = "failed_no_infer"
                job["phase_at"] = now_ts()
            return
        cmd = [
            PYTHON_BIN,
            INFER_PY,
            "--protein_path",
            str(receptor),
            "--ligand",
            str(ligand),
            "--out_dir",
            str(out_dir),
            "--inference_steps",
            inference_steps,
            "--samples_per_complex",
            samples_per_complex,
        ]
        logf.write("CMD: " + " ".join(shlex.quote(x) for x in cmd) + "\n")
        logf.flush()
        try:
            with jobs_lock:
                job["phase"] = "validate_env"
                job["phase_at"] = now_ts()
            proc = subprocess.Popen(
                cmd,
                cwd=str(Path(INFER_PY).parent),
                stdout=logf,
                stderr=subprocess.STDOUT,
                env=os.environ.copy(),
            )
            with jobs_lock:
                job["pid"] = proc.pid
                job["phase"] = "spawn"
                job["phase_at"] = now_ts()
                job["phase"] = "waiting"
                job["phase_at"] = now_ts()
            code = proc.wait()
            with jobs_lock:
                job["phase"] = "collect_results"
                job["phase_at"] = now_ts()
            poses = list_pose_files(out_dir)
            with jobs_lock:
                job["exit_code"] = code
                job["poses"] = poses
                job["status"] = "done" if code == 0 and poses else "failed"
                if code == 0 and not poses:
                    job["error"] = "DiffDock完成但未产出SDF结果"
                    job["phase"] = "failed_no_pose"
                if code != 0:
                    job["error"] = f"DiffDock执行失败，退出码 {code}"
                    job["phase"] = "failed"
                if code == 0 and poses:
                    job["phase"] = "done"
                job["ended_at"] = now_ts()
                job["phase_at"] = now_ts()
        except Exception as e:
            traceback.print_exc(file=logf)
            with jobs_lock:
                job["status"] = "failed"
                job["error"] = str(e)
                job["ended_at"] = now_ts()
                job["phase"] = "failed_exception"
                job["phase_at"] = now_ts()


def run_fpocket_job(job_id: str):
    try:
        with fpocket_jobs_lock:
            job = fpocket_jobs.get(job_id)
            if not job:
                return
            job["status"] = "running"
            job["started_at"] = now_ts()
            job["phase"] = "starting"
            job["phase_at"] = now_ts()
        job_dir = FPOCKET_WORK_DIR / job_id
        with fpocket_jobs_lock:
            job["phase"] = "prepare_dirs"
            job["phase_at"] = now_ts()
        in_dir = job_dir / "input"
        out_dir = job_dir / "output"
        log_path = job_dir / "run.log"
        out_dir.mkdir(parents=True, exist_ok=True)
        pdb_path = in_dir / job["pdb_name"]
        with fpocket_jobs_lock:
            job["phase"] = "discover_bin"
            job["phase_at"] = now_ts()
        fpocket_bin = discover_fpocket_bin()
        with fpocket_jobs_lock:
            job["phase"] = "open_log"
            job["phase_at"] = now_ts()
        with open(log_path, "a", encoding="utf-8") as logf:
            if not fpocket_bin:
                msg = "未找到fpocket，请先安装或设置环境变量 FPOCKET_BIN"
                logf.write(msg + "\n")
                with fpocket_jobs_lock:
                    job["status"] = "failed"
                    job["error"] = msg
                    job["ended_at"] = now_ts()
                    job["phase"] = "failed_no_bin"
                    job["phase_at"] = now_ts()
                return
            cmd = [fpocket_bin, "-f", str(pdb_path)]
            logf.write("CMD: " + " ".join(shlex.quote(x) for x in cmd) + "\n")
            logf.flush()
            try:
                with fpocket_jobs_lock:
                    job["phase"] = "spawn"
                    job["phase_at"] = now_ts()
                proc = subprocess.Popen(
                    cmd,
                    cwd=str(out_dir),
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    env=os.environ.copy(),
                )
                with fpocket_jobs_lock:
                    job["pid"] = proc.pid
                    job["fpocket_bin"] = fpocket_bin
                    job["phase"] = "waiting"
                    job["phase_at"] = now_ts()
                code = proc.wait()
                outputs = list_fpocket_outputs(job_id)
                with fpocket_jobs_lock:
                    job["exit_code"] = code
                    job["status"] = "done" if code == 0 else "failed"
                    job["outputs"] = outputs
                    if code != 0:
                        job["error"] = f"fpocket执行失败，退出码 {code}"
                    job["ended_at"] = now_ts()
                    job["phase"] = "done" if code == 0 else "failed_exit"
                    job["phase_at"] = now_ts()
            except Exception as e:
                traceback.print_exc(file=logf)
                with fpocket_jobs_lock:
                    job["status"] = "failed"
                    job["error"] = str(e)
                    job["ended_at"] = now_ts()
                    job["phase"] = "failed_exception"
                    job["phase_at"] = now_ts()
    except Exception as e:
        with fpocket_jobs_lock:
            job = fpocket_jobs.get(job_id)
            if job:
                job["status"] = "failed"
                job["error"] = f"fpocket任务启动失败: {e}"
                job["ended_at"] = now_ts()
                job["phase"] = "failed_bootstrap"
                job["phase_at"] = now_ts()


def run_tmalign_job(job_id: str):
    with tmalign_jobs_lock:
        job = tmalign_jobs.get(job_id)
        if not job:
            return
        job["status"] = "running"
        job["started_at"] = now_ts()
    job_dir = TMALIGN_WORK_DIR / job_id
    in_dir = job_dir / "input"
    log_path = job_dir / "run.log"
    pdb1 = in_dir / job["pdb1_name"]
    pdb2 = in_dir / job["pdb2_name"]
    tmalign_bin = discover_tmalign_bin()
    with open(log_path, "a", encoding="utf-8") as logf:
        if not tmalign_bin:
            msg = "未找到TMalign，请先安装或设置环境变量 TMALIGN_BIN"
            logf.write(msg + "\n")
            with tmalign_jobs_lock:
                job["status"] = "failed"
                job["error"] = msg
                job["ended_at"] = now_ts()
            return
        cmd = [tmalign_bin, str(pdb1), str(pdb2)]
        logf.write("CMD: " + " ".join(shlex.quote(x) for x in cmd) + "\n")
        logf.flush()
        try:
            proc = subprocess.Popen(
                cmd,
                cwd=str(job_dir),
                stdout=logf,
                stderr=subprocess.STDOUT,
                env=os.environ.copy(),
            )
            with tmalign_jobs_lock:
                job["pid"] = proc.pid
                job["tmalign_bin"] = tmalign_bin
            code = proc.wait()
            log_text = log_path.read_text(encoding="utf-8", errors="ignore") if log_path.exists() else ""
            metrics = parse_tmalign_log(log_text)
            with tmalign_jobs_lock:
                job["exit_code"] = code
                job["status"] = "done" if code == 0 else "failed"
                job["metrics"] = metrics
                if code != 0:
                    job["error"] = f"TMalign执行失败，退出码 {code}"
                job["ended_at"] = now_ts()
        except Exception as e:
            traceback.print_exc(file=logf)
            with tmalign_jobs_lock:
                job["status"] = "failed"
                job["error"] = str(e)
                job["ended_at"] = now_ts()


class Handler(BaseHTTPRequestHandler):
    def _write_json(self, code, payload):
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET,POST,OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()
        self.wfile.write(body)

    def _write_file(self, path: Path):
        if not path.exists() or not path.is_file():
            self._write_json(404, {"ok": False, "error": "file not found"})
            return
        data = path.read_bytes()
        ctype = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(data)

    def do_OPTIONS(self):
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET,POST,OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def do_GET(self):
        parsed = urlparse(self.path)
        path = parsed.path
        if path == "/api/diffdock/ping":
            self._write_json(200, {"ok": True, "service": "diffdock-api", "infer_py": INFER_PY})
            return
        if path == "/api/fpocket/ping":
            self._write_json(200, {
                "ok": True,
                "service": "fpocket-api",
                "fpocket_bin": discover_fpocket_bin(),
                "api_version": "2",
                "features": {
                    "fpocket_phase": True,
                    "fpocket_watchdog": True,
                },
            })
            return
        if path == "/api/tmalign/ping":
            self._write_json(200, {"ok": True, "service": "tmalign-api", "tmalign_bin": discover_tmalign_bin()})
            return
        if path == "/api/pdbzn/workflow/config":
            selection = _pdbzn_active_import_selection()
            self._write_json(200, {
                "ok": True,
                "defaults": _pdbzn_workflow_defaults(),
                "db_path": str(PDBZN_DB_PATH),
                "db_stats": _pdbzn_db_stats(),
                "source_dir": str(PDBZN_BASE_DIR),
                "source_total": _pdbzn_source_total_count(),
                "import_selection": {
                    "active": bool(selection and selection.get("ids")),
                    "source_name": str((selection or {}).get("source_name", "") or ""),
                    "count": len((selection or {}).get("ids") or []),
                    "updated_at_text": str((selection or {}).get("updated_at_text", "") or ""),
                    "ids_preview": list((selection or {}).get("ids") or [])[:20],
                },
                "files": {k: (str(v) if v else "") for k, v in _pdbzn_file_map().items()},
            })
            return
        if path == "/api/data/file":
            rel_path = str((parse_qs(parsed.query or "").get("path") or [""])[0] or "")
            file_path = resolve_under(DATA_DIR, rel_path)
            if (not file_path) or (not file_path.exists()) or (not file_path.is_file()):
                self._write_json(404, {"ok": False, "error": "data file not found"})
                return
            self._write_file(file_path)
            return
        if path.startswith("/api/pdbzn/structure/"):
            pdb_id = _pdbzn_guess_pdb_id(path.split("/")[-1])
            if not pdb_id:
                self._write_json(400, {"ok": False, "error": "invalid pdb id"})
                return
            file_path = _pdbzn_find_structure_path(pdb_id) or _pdbzn_find_tmalign_pdb_path(pdb_id)
            if file_path is None:
                self._write_json(404, {"ok": False, "error": "structure not found"})
                return
            self._write_file(file_path)
            return
        if path == "/api/diffdock/jobs":
            with jobs_lock:
                arr = list(jobs.values())
            arr.sort(key=lambda x: x.get("created_at", 0), reverse=True)
            self._write_json(200, {"ok": True, "jobs": arr[:50]})
            return
        if path == "/api/fpocket/jobs":
            with fpocket_jobs_lock:
                arr = list(fpocket_jobs.values())
            arr.sort(key=lambda x: x.get("created_at", 0), reverse=True)
            self._write_json(200, {"ok": True, "jobs": arr[:50]})
            return
        if path == "/api/tmalign/jobs":
            with tmalign_jobs_lock:
                arr = list(tmalign_jobs.values())
            arr.sort(key=lambda x: x.get("created_at", 0), reverse=True)
            self._write_json(200, {"ok": True, "jobs": arr[:50]})
            return
        if path.startswith("/api/diffdock/status/"):
            job_id = path.split("/")[-1]
            with jobs_lock:
                job = jobs.get(job_id)
            if not job:
                self._write_json(404, {"ok": False, "error": "job not found"})
                return
            self._write_json(200, {"ok": True, "job": job})
            return
        if path.startswith("/api/fpocket/status/"):
            job_id = path.split("/")[-1]
            with fpocket_jobs_lock:
                job = fpocket_jobs.get(job_id)
            if not job:
                job = recover_fpocket_job(job_id)
                if not job:
                    self._write_json(404, {"ok": False, "error": "job not found"})
                    return
            if job.get("status") == "running":
                started = int(job.get("started_at") or 0)
                phase = str(job.get("phase") or "")
                now = now_ts()
                if (job.get("pid") is None) and started and (now - started > 30) and phase not in {"waiting", "done"}:
                    with fpocket_jobs_lock:
                        j2 = fpocket_jobs.get(job_id) or job
                        if j2.get("status") == "running" and j2.get("pid") is None:
                            j2["status"] = "failed"
                            j2["ended_at"] = now
                            j2["phase"] = "failed_watchdog"
                            j2["phase_at"] = now
                            j2["error"] = f"fpocket任务卡住（phase={phase or '-'}，30s无pid），请重启后端后重试"
                            job = j2
            self._write_json(200, {"ok": True, "job": job})
            return
        if path.startswith("/api/fpocket/pockets/"):
            job_id = path.split("/")[-1]
            with fpocket_jobs_lock:
                job = fpocket_jobs.get(job_id)
            if not job:
                job = recover_fpocket_job(job_id)
                if not job:
                    self._write_json(404, {"ok": False, "error": "job not found"})
                    return
            payload = build_fpocket_pocket_payload(job_id, job)
            if not payload.get("ok"):
                self._write_json(400, payload)
                return
            self._write_json(200, payload)
            return
        if path.startswith("/api/tmalign/status/"):
            job_id = path.split("/")[-1]
            with tmalign_jobs_lock:
                job = tmalign_jobs.get(job_id)
            if not job:
                self._write_json(404, {"ok": False, "error": "job not found"})
                return
            self._write_json(200, {"ok": True, "job": job})
            return
        if path.startswith("/api/diffdock/file/"):
            parts = path.split("/")
            if len(parts) < 6:
                self._write_json(400, {"ok": False, "error": "invalid file path"})
                return
            job_id = parts[4]
            fname = safe_name(parts[5])
            file_path = WORK_DIR / job_id / "output" / fname
            self._write_file(file_path)
            return
        if path.startswith("/api/fpocket/file/"):
            prefix = "/api/fpocket/file/"
            tail = path[len(prefix):]
            seg = tail.split("/", 1)
            if len(seg) != 2:
                self._write_json(400, {"ok": False, "error": "invalid file path"})
                return
            job_id, rel_path = seg[0], seg[1]
            job_dir = FPOCKET_WORK_DIR / job_id
            file_path = resolve_under(job_dir, rel_path)
            if (not file_path) or (not file_path.exists()) or (not file_path.is_file()):
                self._write_json(404, {"ok": False, "error": "file not found"})
                return
            input_dir = (job_dir / "input").resolve()
            output_dir = (job_dir / "output").resolve()
            in_allowed = str(file_path).startswith(str(input_dir))
            out_allowed = str(file_path).startswith(str(output_dir))
            if not (in_allowed or out_allowed):
                self._write_json(400, {"ok": False, "error": "invalid file path"})
                return
            self._write_file(file_path)
            return
        if path.startswith("/api/diffdock/log/"):
            job_id = path.split("/")[-1]
            log_path = WORK_DIR / job_id / "run.log"
            if not log_path.exists():
                self._write_json(404, {"ok": False, "error": "log not found"})
                return
            self._write_json(200, {"ok": True, "log": log_path.read_text(encoding="utf-8", errors="ignore")})
            return
        if path.startswith("/api/fpocket/log/"):
            job_id = path.split("/")[-1]
            log_path = FPOCKET_WORK_DIR / job_id / "run.log"
            if not log_path.exists():
                self._write_json(404, {"ok": False, "error": "log not found"})
                return
            self._write_json(200, {"ok": True, "log": log_path.read_text(encoding="utf-8", errors="ignore")})
            return
        if path.startswith("/api/tmalign/log/"):
            job_id = path.split("/")[-1]
            log_path = TMALIGN_WORK_DIR / job_id / "run.log"
            if not log_path.exists():
                self._write_json(404, {"ok": False, "error": "log not found"})
                return
            self._write_json(200, {"ok": True, "log": log_path.read_text(encoding="utf-8", errors="ignore")})
            return
        self._write_json(404, {"ok": False, "error": "not found"})

    def do_POST(self):
        parsed = urlparse(self.path)
        ctype = (self.headers.get("content-type", "") or "").split(";")[0].strip().lower()
        if ctype != "application/json":
            self._write_json(400, {"ok": False, "error": "content-type must be application/json"})
            return
        try:
            n = int(self.headers.get("content-length", "0"))
            raw = self.rfile.read(n)
            payload = json.loads(raw.decode("utf-8"))
        except Exception:
            self._write_json(400, {"ok": False, "error": "invalid json body"})
            return
        if parsed.path == "/api/diffdock/submit":
            receptor_name = safe_name(payload.get("receptor_name", "receptor.pdb"))
            ligand_name = safe_name(payload.get("ligand_name", "ligand.sdf"))
            receptor_text = payload.get("receptor_text", "")
            ligand_text = payload.get("ligand_text", "")
            if not receptor_text or not ligand_text:
                self._write_json(400, {"ok": False, "error": "missing receptor_text or ligand_text"})
                return
            job_id = uuid.uuid4().hex[:12]
            job_dir = WORK_DIR / job_id
            in_dir = job_dir / "input"
            in_dir.mkdir(parents=True, exist_ok=True)
            receptor_path = in_dir / receptor_name
            ligand_path = in_dir / ligand_name
            receptor_path.write_text(receptor_text, encoding="utf-8")
            ligand_path.write_text(ligand_text, encoding="utf-8")
            inf_steps = 20
            sample_n = 10
            try:
                inf_steps = int(payload.get("inference_steps", 20))
            except Exception:
                pass
            try:
                sample_n = int(payload.get("samples_per_complex", 10))
            except Exception:
                pass
            job = {
                "id": job_id,
                "status": "queued",
                "created_at": now_ts(),
                "receptor_name": receptor_name,
                "ligand_name": ligand_name,
                "inference_steps": inf_steps,
                "samples_per_complex": sample_n,
                "poses": [],
                "error": "",
                "pid": None,
                "phase": "queued",
                "phase_at": now_ts(),
            }
            with jobs_lock:
                jobs[job_id] = job
            t = threading.Thread(target=run_job, args=(job_id,), daemon=True)
            t.start()
            self._write_json(200, {"ok": True, "job_id": job_id, "status_url": f"/api/diffdock/status/{job_id}"})
            return
        if parsed.path == "/api/fpocket/submit":
            pdb_name = safe_name(payload.get("pdb_name", "receptor.pdb"))
            pdb_text = payload.get("pdb_text", "")
            if not pdb_text:
                self._write_json(400, {"ok": False, "error": "missing pdb_text"})
                return
            if not (pdb_name.lower().endswith(".pdb") or pdb_name.lower().endswith(".ent")):
                self._write_json(400, {"ok": False, "error": "pdb_name必须是.pdb或.ent"})
                return
            job_id = uuid.uuid4().hex[:12]
            job_dir = FPOCKET_WORK_DIR / job_id
            in_dir = job_dir / "input"
            in_dir.mkdir(parents=True, exist_ok=True)
            pdb_path = in_dir / pdb_name
            pdb_path.write_text(pdb_text, encoding="utf-8")
            job = {
                "id": job_id,
                "status": "queued",
                "created_at": now_ts(),
                "pdb_name": pdb_name,
                "outputs": [],
                "error": "",
                "pid": None,
                "fpocket_bin": "",
                "phase": "queued",
                "phase_at": now_ts(),
            }
            with fpocket_jobs_lock:
                fpocket_jobs[job_id] = job
            t = threading.Thread(target=run_fpocket_job, args=(job_id,), daemon=True)
            t.start()
            self._write_json(200, {"ok": True, "job_id": job_id, "status_url": f"/api/fpocket/status/{job_id}"})
            return
        if parsed.path == "/api/tmalign/submit":
            pdb1_name = safe_name(payload.get("pdb1_name", "struct1.pdb"))
            pdb2_name = safe_name(payload.get("pdb2_name", "struct2.pdb"))
            pdb1_text = payload.get("pdb1_text", "")
            pdb2_text = payload.get("pdb2_text", "")
            if not pdb1_text or not pdb2_text:
                self._write_json(400, {"ok": False, "error": "missing pdb1_text or pdb2_text"})
                return
            if not (pdb1_name.lower().endswith(".pdb") or pdb1_name.lower().endswith(".ent")):
                self._write_json(400, {"ok": False, "error": "pdb1_name必须是.pdb或.ent"})
                return
            if not (pdb2_name.lower().endswith(".pdb") or pdb2_name.lower().endswith(".ent")):
                self._write_json(400, {"ok": False, "error": "pdb2_name必须是.pdb或.ent"})
                return
            job_id = uuid.uuid4().hex[:12]
            job_dir = TMALIGN_WORK_DIR / job_id
            in_dir = job_dir / "input"
            in_dir.mkdir(parents=True, exist_ok=True)
            pdb1_path = in_dir / pdb1_name
            pdb2_path = in_dir / pdb2_name
            pdb1_path.write_text(pdb1_text, encoding="utf-8")
            pdb2_path.write_text(pdb2_text, encoding="utf-8")
            job = {
                "id": job_id,
                "status": "queued",
                "created_at": now_ts(),
                "pdb1_name": pdb1_name,
                "pdb2_name": pdb2_name,
                "metrics": {},
                "error": "",
                "pid": None,
                "tmalign_bin": "",
            }
            with tmalign_jobs_lock:
                tmalign_jobs[job_id] = job
            t = threading.Thread(target=run_tmalign_job, args=(job_id,), daemon=True)
            t.start()
            self._write_json(200, {"ok": True, "job_id": job_id, "status_url": f"/api/tmalign/status/{job_id}"})
            return
        if parsed.path == "/api/pdbzn/workflow/run":
            result = _pdbzn_run_workflow(payload)
            if result.get("ok"):
                self._write_json(200, result)
                return
            self._write_json(400, result)
            return
        if parsed.path == "/api/pdbzn/workflow/cluster":
            result = _pdbzn_cluster_workflow(payload)
            if result.get("ok"):
                self._write_json(200, result)
                return
            self._write_json(400, result)
            return
        if parsed.path == "/api/pdbzn/workflow/filter":
            if isinstance(payload, dict) and ("step4" in payload):
                result = _pdbzn_step4_validate_workflow(payload)
            else:
                result = _pdbzn_step3_filter_workflow(payload)
            if result.get("ok"):
                self._write_json(200, result)
                return
            self._write_json(400, result)
            return
        if parsed.path == "/api/pdbzn/workflow/validate":
            result = _pdbzn_step4_validate_workflow(payload)
            if result.get("ok"):
                self._write_json(200, result)
                return
            self._write_json(400, result)
            return
        if parsed.path == "/api/pdbzn/workflow/finalize":
            try:
                result = _pdbzn_step5_finalize_workflow(payload)
            except Exception as e:
                traceback.print_exc()
                self._write_json(500, {"ok": False, "error": f"step5_internal_error: {e}"})
                return
            if result.get("ok"):
                self._write_json(200, result)
                return
            self._write_json(400, result)
            return
        if parsed.path == "/api/pdbzn/workflow/import":
            result = _pdbzn_import_database(payload)
            if result.get("ok"):
                self._write_json(200, result)
                return
            self._write_json(400, result)
            return
        self._write_json(404, {"ok": False, "error": "not found"})


def run():
    server = ThreadingHTTPServer((HOST, PORT), Handler)
    print(f"DiffDock API listening on http://{HOST}:{PORT}")
    server.serve_forever()


if __name__ == "__main__":
    run()
