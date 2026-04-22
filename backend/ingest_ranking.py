import csv
import json
import os
import shutil
import urllib.request
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = BASE_DIR.parent
DATA_DIR = BASE_DIR / "data"
RANKING_CSV = Path(os.environ.get("PDBZN_RANKING_CSV", PROJECT_DIR / "PDB_ZN" / "zn_his_similarity_ranking2.csv"))
OUT_STRUCT_DIR = DATA_DIR / "structures"
OUT_JSON = DATA_DIR / "ranking.json"

def read_ids():
    ids = []
    if not RANKING_CSV.exists():
        raise SystemExit("ranking csv not found")
    with RANKING_CSV.open(newline="", encoding="utf-8", errors="ignore") as f:
        rdr = csv.reader(f)
        for row in rdr:
            if not row: continue
            # heuristics: second column is id; fallback search token len==4 alnum
            cand = None
            if len(row) >= 2 and row[1] and len(row[1]) == 4:
                cand = row[1].strip()
            else:
                for t in row:
                    t = (t or "").strip()
                    if len(t) == 4 and t.isalnum():
                        cand = t
                        break
            if not cand: continue
            ids.append((cand, row[2] if len(row) > 2 else ""))
    return ids

def ensure_structure(pdbid):
    OUT_STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    root = PROJECT_DIR / "PDB_ZN" / "diffdock_inputs" / "receptors"
    local_pdb = root / f"{pdbid}.pdb"
    local_cif = root / f"{pdbid}.cif"
    if local_pdb.exists():
        dst = OUT_STRUCT_DIR / f"{pdbid}.pdb"
        shutil.copy(str(local_pdb), str(dst))
        return f"../backend/data/structures/{pdbid}.pdb"
    if local_cif.exists():
        dst = OUT_STRUCT_DIR / f"{pdbid}.cif"
        shutil.copy(str(local_cif), str(dst))
        return f"../backend/data/structures/{pdbid}.cif"
    # try download cif
    url = f"https://files.rcsb.org/download/{pdbid}.cif"
    dst = OUT_STRUCT_DIR / f"{pdbid}.cif"
    try:
        urllib.request.urlretrieve(url, str(dst))
        return f"../backend/data/structures/{pdbid}.cif"
    except Exception:
        return ""

def main():
    items = []
    for pdbid, score in read_ids():
        rel = ensure_structure(pdbid)
        items.append({"id": pdbid, "score": score, "receptor_rel": rel})
    with OUT_JSON.open("w", encoding="utf-8") as f:
        json.dump({"items": items}, f, ensure_ascii=False, indent=2)
    print("OK", OUT_JSON, len(items))

if __name__ == "__main__":
    main()
