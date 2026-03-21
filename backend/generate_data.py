import os, csv, json

def load_master():
    paths = [
        "/root/PDB_ZN/zn_his_master_table7.csv",
        "/root/zn_his_master_table.csv",
    ]
    for p in paths:
        if os.path.exists(p):
            with open(p, newline="") as f:
                rows = list(csv.reader(f))
            return p, rows
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

def main():
    master_path, rows = load_master()
    header = rows[0]
    idx = {name: i for i, name in enumerate(header)}
    out = []
    for r in rows[1:]:
        def get(name):
            i = idx.get(name)
            return r[i] if i is not None and i < len(r) else ""
        rep = get("Representative") or os.path.splitext(os.path.basename(get("Receptor_PDB")))[0]
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
        }
        out.append(item)
    data_path = "/root/data-platform/backend/data/data.json"
    os.makedirs(os.path.dirname(data_path), exist_ok=True)
    with open(data_path, "w") as f:
        json.dump({"items": out}, f, ensure_ascii=False, indent=2)
    print("OK", data_path, len(out))
    full = []
    for r in rows[1:]:
        obj = {}
        for i, k in enumerate(header):
            obj[k] = r[i] if i < len(r) else ""
        rep_idx = idx.get("Representative")
        if rep_idx is not None and rep_idx < len(r):
            obj["_id"] = r[rep_idx]
        else:
            rp_idx = idx.get("Receptor_PDB")
            if rp_idx is not None and rp_idx < len(r):
                base = os.path.basename(r[rp_idx]) if r[rp_idx] else ""
                obj["_id"] = os.path.splitext(base)[0]
        full.append(obj)
    full_path = "/root/data-platform/backend/data/master_full.json"
    with open(full_path, "w") as f:
        json.dump({"rows": full, "header": header}, f, ensure_ascii=False, indent=2)
    print("OK", full_path, len(full))

if __name__ == "__main__":
    main()
