#!/usr/bin/env python3
"""
Generate a non-breaking organized view of project data.

The application still reads from backend/data and backend/runtime.
This script creates backend/data_catalog as a clearer, category-based
layout that points back to the original files and directories.
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
BACKEND_DIR = PROJECT_DIR / "backend"
DATA_DIR = BACKEND_DIR / "data"
RUNTIME_DIR = BACKEND_DIR / "runtime"
CATALOG_DIR = BACKEND_DIR / "data_catalog"
REGISTRY_PATH = DATA_DIR / "table_registry.json"


def reset_path(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.exists():
        shutil.rmtree(path)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def link_path(link_path: Path, target_path: Path) -> None:
    reset_path(link_path)
    ensure_dir(link_path.parent)
    rel_target = os.path.relpath(target_path, link_path.parent)
    os.symlink(rel_target, link_path)


def maybe_link(link_path: Path, target_path: Path) -> bool:
    if not target_path.exists():
        return False
    link_path.parent.mkdir(parents=True, exist_ok=True)
    link_path_rel_target = os.path.relpath(target_path, link_path.parent)
    reset_path(link_path)
    os.symlink(link_path_rel_target, link_path)
    return True


def load_registry() -> dict:
    if not REGISTRY_PATH.exists():
        return {"tables": []}
    return json.loads(REGISTRY_PATH.read_text())


def prepare_catalog_root() -> None:
    ensure_dir(CATALOG_DIR)
    for name in [
        "manifests",
        "tables",
        "structure_assets",
        "docking_assets",
        "pocket_assets",
        "runtime_assets",
    ]:
        reset_path(CATALOG_DIR / name)
        ensure_dir(CATALOG_DIR / name)


def build_tables_view(registry: dict) -> list[dict]:
    tables_root = CATALOG_DIR / "tables"
    ensure_dir(tables_root / "default_views")
    ensure_dir(tables_root / "bundles")

    maybe_link(tables_root / "table_registry.json", DATA_DIR / "table_registry.json")
    maybe_link(tables_root / "default_views" / "data.json", DATA_DIR / "data.json")
    maybe_link(tables_root / "default_views" / "master_full.json", DATA_DIR / "master_full.json")
    maybe_link(tables_root / "default_views" / "master_table.csv", DATA_DIR / "master_table.csv")

    bundled = []
    for entry in registry.get("tables", []):
        table_id = str(entry.get("id", "")).strip()
        if not table_id:
            continue
        bundle_dir = tables_root / "bundles" / table_id
        ensure_dir(bundle_dir)

        mapping = {
            "data.json": DATA_DIR / str(entry.get("data", "")).strip(),
            "master_full.json": DATA_DIR / str(entry.get("full", "")).strip(),
            "master_table.csv": DATA_DIR / str(entry.get("csv", "")).strip(),
            "ids.json": DATA_DIR / str(entry.get("ids", "")).strip(),
        }

        resolved = {}
        for alias, target in mapping.items():
            if target.name:
                linked = maybe_link(bundle_dir / alias, target)
                if linked:
                    resolved[alias] = str(target.relative_to(PROJECT_DIR))

        meta = dict(entry)
        meta["organized_bundle"] = str(bundle_dir.relative_to(PROJECT_DIR))
        meta["files"] = resolved
        bundled.append(meta)

    return bundled


def build_asset_views() -> dict:
    structure_root = CATALOG_DIR / "structure_assets"
    docking_root = CATALOG_DIR / "docking_assets"
    pocket_root = CATALOG_DIR / "pocket_assets"
    runtime_root = CATALOG_DIR / "runtime_assets"

    maybe_link(structure_root / "structures", DATA_DIR / "structures")
    maybe_link(structure_root / "ligands", DATA_DIR / "ligands")

    maybe_link(docking_root / "diffdock_index.json", DATA_DIR / "diffdock_index.json")
    maybe_link(docking_root / "diffdock_results", DATA_DIR / "diffdock")

    maybe_link(pocket_root / "pocket.csv", DATA_DIR / "pocket.csv")
    maybe_link(pocket_root / "pocket_results", DATA_DIR / "pockets")

    maybe_link(runtime_root / "pdbzn.sqlite", RUNTIME_DIR / "pdbzn.sqlite")
    maybe_link(runtime_root / "exports", RUNTIME_DIR / "exports")
    maybe_link(runtime_root / "diffdock_jobs", RUNTIME_DIR / "diffdock_jobs")
    maybe_link(runtime_root / "fpocket_jobs", RUNTIME_DIR / "fpocket_jobs")
    maybe_link(runtime_root / "tmalign_jobs", RUNTIME_DIR / "tmalign_jobs")
    maybe_link(runtime_root / "workflow_jobs", RUNTIME_DIR / "workflow_jobs")
    maybe_link(runtime_root / "import_structures", RUNTIME_DIR / "import_structures")
    maybe_link(runtime_root / "logs", RUNTIME_DIR / "logs")
    maybe_link(runtime_root / "legacy_pdbzn_workflow.db", BACKEND_DIR / "pdbzn_workflow.db")

    return {
        "structures_dir": "backend/data/structures",
        "ligands_dir": "backend/data/ligands",
        "diffdock_dir": "backend/data/diffdock",
        "pockets_dir": "backend/data/pockets",
        "runtime_dir": "backend/runtime",
    }


def count_items() -> dict:
    return {
        "registered_tables": len(load_registry().get("tables", [])),
        "default_master_rows_file": "backend/data/master_full.json",
        "structures_pdb": len(list((DATA_DIR / "structures").glob("*.pdb"))),
        "structures_cif": len(list((DATA_DIR / "structures").glob("*.cif"))),
        "ligands_sdf": len(list((DATA_DIR / "ligands").glob("*.sdf"))),
        "diffdock_result_dirs": len([p for p in (DATA_DIR / "diffdock").iterdir() if p.is_dir()]),
        "pocket_result_dirs": len([p for p in (DATA_DIR / "pockets").iterdir() if p.is_dir()]),
    }


def write_summary_manifest(tables: list[dict], assets: dict) -> None:
    summary = {
        "generated_from": "scripts/organize_data_catalog.py",
        "catalog_root": str(CATALOG_DIR.relative_to(PROJECT_DIR)),
        "counts": count_items(),
        "tables": tables,
        "asset_dirs": assets,
    }
    manifest_path = CATALOG_DIR / "manifests" / "summary.json"
    manifest_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")


def main() -> None:
    registry = load_registry()
    prepare_catalog_root()
    tables = build_tables_view(registry)
    assets = build_asset_views()
    write_summary_manifest(tables, assets)
    print(f"OK {CATALOG_DIR}")


if __name__ == "__main__":
    main()
