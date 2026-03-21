import os, csv, json, shutil, urllib.request

RANKING_CSV = "/root/PDB_ZN/zn_his_similarity_ranking2.csv"
OUT_STRUCT_DIR = "/root/data-platform/backend/data/structures"
OUT_JSON = "/root/data-platform/backend/data/ranking.json"

def read_ids():
    ids = []
    if not os.path.exists(RANKING_CSV):
        raise SystemExit("ranking csv not found")
    with open(RANKING_CSV, newline="") as f:
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
    os.makedirs(OUT_STRUCT_DIR, exist_ok=True)
    local_pdb = f"/root/PDB_ZN/diffdock_inputs/receptors/{pdbid}.pdb"
    local_cif = f"/root/PDB_ZN/diffdock_inputs/receptors/{pdbid}.cif"
    if os.path.exists(local_pdb):
        dst = os.path.join(OUT_STRUCT_DIR, f"{pdbid}.pdb")
        shutil.copy(local_pdb, dst)
        return f"../backend/data/structures/{pdbid}.pdb"
    if os.path.exists(local_cif):
        dst = os.path.join(OUT_STRUCT_DIR, f"{pdbid}.cif")
        shutil.copy(local_cif, dst)
        return f"../backend/data/structures/{pdbid}.cif"
    # try download cif
    url = f"https://files.rcsb.org/download/{pdbid}.cif"
    dst = os.path.join(OUT_STRUCT_DIR, f"{pdbid}.cif")
    try:
        urllib.request.urlretrieve(url, dst)
        return f"../backend/data/structures/{pdbid}.cif"
    except Exception:
        return ""

def main():
    items = []
    for pdbid, score in read_ids():
        rel = ensure_structure(pdbid)
        items.append({"id": pdbid, "score": score, "receptor_rel": rel})
    with open(OUT_JSON, "w") as f:
        json.dump({"items": items}, f, ensure_ascii=False, indent=2)
    print("OK", OUT_JSON, len(items))

if __name__ == "__main__":
    main()
