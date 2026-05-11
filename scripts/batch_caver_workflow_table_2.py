#!/usr/bin/env python3
import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "backend" / "data"
RUNTIME_DIR = PROJECT_DIR / "backend" / "runtime"
JOB_ROOT = RUNTIME_DIR / "caver_batch_workflow_table_2"
REPORT_DIR = RUNTIME_DIR / "reports"

CAVER_HOME = PROJECT_DIR / ".tools" / "src" / "caver_3.0" / "caver_3.0" / "caver"
CAVER_JAR = CAVER_HOME / "caver.jar"
JAVA_BIN = shutil.which("java") or "java"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def parse_first_zn_xyz(pdb_text: str):
    for line in str(pdb_text or "").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        atom_name = str(line[12:16] or "").strip().upper()
        resn = str(line[17:20] or "").strip().upper()
        elem = str(line[76:78] or "").strip().upper() or atom_name
        if not (elem == "ZN" or resn == "ZN" or atom_name == "ZN"):
            continue
        x = float(str(line[30:38] or "").strip())
        y = float(str(line[38:46] or "").strip())
        z = float(str(line[46:54] or "").strip())
        return (x, y, z)
    raise ValueError("No ZN atom found")


def parse_best_pocket_center(pid: str):
    pocket_path = DATA_DIR / "pockets" / pid / "best_pocket_atm.pdb"
    xs = []
    ys = []
    zs = []
    for line in pocket_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            xs.append(float(str(line[30:38] or "").strip()))
            ys.append(float(str(line[38:46] or "").strip()))
            zs.append(float(str(line[46:54] or "").strip()))
        except Exception:
            continue
    if not xs:
        raise ValueError(f"No pocket atoms found for {pid}")
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def write_config(path: Path, x: float, y: float, z: float):
    path.write_text(
        "\n".join(
            [
                f"starting_point_coordinates {x:.3f} {y:.3f} {z:.3f}",
                "probe_radius 0.9",
                "shell_radius 3",
                "shell_depth 4",
                "frame_clustering_threshold 1",
                "clustering_threshold 3.5",
                "seed 1",
                "",
            ]
        ),
        encoding="utf-8",
    )


def run_caver_for_pid(pid: str, seed_mode: str, probe_radius: float, shell_radius: float, shell_depth: float, report_tag: str):
    pdb_path = DATA_DIR / "structures" / f"{pid}.pdb"
    if not pdb_path.exists():
        raise FileNotFoundError(str(pdb_path))
    pdb_text = pdb_path.read_text(encoding="utf-8", errors="ignore")
    if seed_mode == "pocket":
        x, y, z = parse_best_pocket_center(pid)
    else:
        x, y, z = parse_first_zn_xyz(pdb_text)

    job_dir = JOB_ROOT / f"{report_tag}_{pid}"
    input_dir = job_dir / "input"
    output_dir = job_dir / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    input_pdb = input_dir / f"{pid}.pdb"
    input_pdb.write_text(pdb_text, encoding="utf-8")
    config_path = input_dir / "config.txt"
    config_path.write_text(
        "\n".join(
            [
                f"starting_point_coordinates {x:.3f} {y:.3f} {z:.3f}",
                f"probe_radius {probe_radius}",
                f"shell_radius {shell_radius}",
                f"shell_depth {shell_depth}",
                "frame_clustering_threshold 1",
                "clustering_threshold 3.5",
                "seed 1",
                "",
            ]
        ),
        encoding="utf-8",
    )

    summary_path = output_dir / "summary_precise_numbers.csv"
    if summary_path.exists():
        return {
            "pid": pid,
            "seed_x": round(x, 3),
            "seed_y": round(y, 3),
            "seed_z": round(z, 3),
            "job_dir": str(job_dir),
        }

    cmd = [
        JAVA_BIN,
        "-Xmx2g",
        "-jar",
        str(CAVER_JAR.resolve()),
        "-home",
        str(CAVER_HOME.resolve()),
        "-pdb",
        str(input_dir.resolve()),
        "-conf",
        str(config_path.resolve()),
        "-out",
        str(output_dir.resolve()),
    ]
    log_path = job_dir / "run.log"
    with log_path.open("w", encoding="utf-8") as logf:
        proc = subprocess.run(cmd, cwd=str(job_dir), stdout=logf, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"CAVER failed for {pid}, see {log_path}")
    whole = job_dir / "whole_clustering.csv"
    if whole.exists():
        target = output_dir / "analysis" / "whole_clustering.csv"
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(whole), str(target))
    return {
        "pid": pid,
        "seed_x": round(x, 3),
        "seed_y": round(y, 3),
        "seed_z": round(z, 3),
        "seed_mode": seed_mode,
        "job_dir": str(job_dir),
    }


def read_csv_rows(path: Path):
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = [str(x or "").strip() for x in (reader.fieldnames or [])]
        rows = []
        for raw in reader:
            clean = {}
            for key, value in raw.items():
                clean[str(key or "").strip()] = value
            rows.append(clean)
        return fieldnames, rows


def classify_radius(radius, pass_radius: float, borderline_radius: float):
    if radius >= pass_radius:
        return "yes"
    if radius >= borderline_radius:
        return "borderline"
    return "no"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed-mode", choices=["zn", "pocket"], default="zn")
    parser.add_argument("--probe-radius", type=float, default=0.9)
    parser.add_argument("--shell-radius", type=float, default=3.0)
    parser.add_argument("--shell-depth", type=float, default=4.0)
    parser.add_argument("--pass-radius", type=float, default=2.4)
    parser.add_argument("--borderline-radius", type=float, default=2.1)
    parser.add_argument("--report-tag", default="zn_seed")
    args = parser.parse_args()

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    JOB_ROOT.mkdir(parents=True, exist_ok=True)

    ids = read_json(DATA_DIR / "workflow_table_2_ids.json")["ids"]
    wf = read_json(DATA_DIR / "workflow_table_2_master_full.json")
    table_rows = {r.get("Representative") or r.get("PDB_ID") or r.get("_id"): r for r in wf["rows"]}

    channel_rows = []
    protein_rows = []
    for pid in ids:
        meta = run_caver_for_pid(
            pid,
            seed_mode=args.seed_mode,
            probe_radius=args.probe_radius,
            shell_radius=args.shell_radius,
            shell_depth=args.shell_depth,
            report_tag=args.report_tag,
        )
        out_dir = Path(meta["job_dir"]) / "output"
        tcsv = out_dir / "analysis" / "tunnel_characteristics.csv"
        scsv = out_dir / "summary_precise_numbers.csv"
        _, rows = read_csv_rows(tcsv)
        _, summary_rows = read_csv_rows(scsv)
        analysis_rows = rows if rows else summary_rows
        table_row = table_rows.get(pid, {})
        if not analysis_rows:
            protein_rows.append(
                {
                    "PDB_ID": pid,
                    "Channel_count": 0,
                    "Best_bottleneck_radius_A": "",
                    "Best_throughput": "",
                    "Acetophenone_size_pass_by_caver": "no",
                    "Caver_seed_mode": seed_mode_label(args.seed_mode),
                    "Existing_acetophenone_judgement": table_row.get("苯乙酮_综合判断", ""),
                    "Existing_size_enterable": table_row.get("苯乙酮_尺寸可进入", ""),
                    "Existing_pose_O_to_Zn_A": table_row.get("苯乙酮_O到Zn最短距离_A", ""),
                    "Existing_pose_any_atom_to_Zn_A": table_row.get("苯乙酮_任意原子到Zn最短距离_A", ""),
                    "Existing_note": "CAVER 从默认 Zn seed 未找到可聚类主通道；需要人工重选起点或放宽参数。" + (
                        (" 原始备注：" + str(table_row.get("苯乙酮_判断说明", "") or "").strip())
                        if str(table_row.get("苯乙酮_判断说明", "") or "").strip()
                        else ""
                    ),
                    "Seed_xyz": f"{meta['seed_x']},{meta['seed_y']},{meta['seed_z']}",
                    "Job_dir": meta["job_dir"],
                }
            )
            continue
        best_radius = None
        best_throughput = None
        for row in analysis_rows:
            try:
                bottleneck = float(str(row.get("Bottleneck radius", "") or "").strip())
            except Exception:
                bottleneck = None
            try:
                throughput = float(str(row.get("Throughput", row.get("Average throughput", "")) or "").strip())
            except Exception:
                throughput = None
            channel_pass = classify_radius(bottleneck, args.pass_radius, args.borderline_radius) if bottleneck is not None else "unknown"
            channel_rows.append(
                {
                    "PDB_ID": pid,
                    "Tunnel cluster": row.get("Tunnel cluster", row.get("Tunnel cluster ID", "")),
                    "Tunnel": row.get("Tunnel", ""),
                    "Throughput": row.get("Throughput", row.get("Average throughput", "")),
                    "Cost": row.get("Cost", ""),
                    "Bottleneck radius": row.get("Bottleneck radius", ""),
                    "Length": row.get("Length", ""),
                    "Curvature": row.get("Curvature", ""),
                    "Acetophenone_size_pass_by_caver": channel_pass,
                }
            )
            if bottleneck is not None and (best_radius is None or bottleneck > best_radius):
                best_radius = bottleneck
            if throughput is not None and (best_throughput is None or throughput > best_throughput):
                best_throughput = throughput
        protein_rows.append(
            {
                "PDB_ID": pid,
                "Channel_count": len(analysis_rows),
                "Best_bottleneck_radius_A": "" if best_radius is None else f"{best_radius:.3f}",
                "Best_throughput": "" if best_throughput is None else f"{best_throughput:.6f}",
                "Acetophenone_size_pass_by_caver": classify_radius(best_radius, args.pass_radius, args.borderline_radius) if best_radius is not None else "unknown",
                "Caver_seed_mode": seed_mode_label(args.seed_mode),
                "Existing_acetophenone_judgement": table_row.get("苯乙酮_综合判断", ""),
                "Existing_size_enterable": table_row.get("苯乙酮_尺寸可进入", ""),
                "Existing_pose_O_to_Zn_A": table_row.get("苯乙酮_O到Zn最短距离_A", ""),
                "Existing_pose_any_atom_to_Zn_A": table_row.get("苯乙酮_任意原子到Zn最短距离_A", ""),
                "Existing_note": table_row.get("苯乙酮_判断说明", ""),
                "Seed_xyz": f"{meta['seed_x']},{meta['seed_y']},{meta['seed_z']}",
                "Job_dir": meta["job_dir"],
            }
        )

    channel_csv = REPORT_DIR / f"workflow_table_2_caver_{args.report_tag}_channels.csv"
    protein_csv = REPORT_DIR / f"workflow_table_2_caver_{args.report_tag}_acetophenone_summary.csv"
    with channel_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(channel_rows[0].keys()) if channel_rows else ["PDB_ID"])
        writer.writeheader()
        for row in channel_rows:
            writer.writerow(row)
    with protein_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(protein_rows[0].keys()) if protein_rows else ["PDB_ID"])
        writer.writeheader()
        for row in protein_rows:
            writer.writerow(row)
    print("CHANNEL_CSV", channel_csv)
    print("PROTEIN_CSV", protein_csv)
    print("COUNT", len(protein_rows))


def seed_mode_label(mode: str):
    return "best_pocket_center" if mode == "pocket" else "zn_atom"


if __name__ == "__main__":
    main()
