# Data Catalog

This directory is an organized view of project data.

The application still reads from:

- `backend/data/`
- `backend/runtime/`

This catalog is only for human browsing and management. It groups the same
files into clearer categories without changing the paths used by the code.

Some runtime links appear only when the corresponding runtime files already
exist. For example, `pdbzn.sqlite` may be absent before the workflow creates it.

## Layout

- `tables/`
  - `table_registry.json`
  - `default_views/`
  - `bundles/<table_id>/`
- `structure_assets/`
  - `structures`
  - `ligands`
- `docking_assets/`
  - `diffdock_index.json`
  - `diffdock_results`
- `pocket_assets/`
  - `pocket.csv`
  - `pocket_results`
- `runtime_assets/`
  - `pdbzn.sqlite`
  - `exports`
  - `*_jobs`
  - `logs`
- `manifests/summary.json`

## Refresh

Run:

```bash
python3 scripts/organize_data_catalog.py
```

The script recreates the organized view and refreshes the manifest.
