#!/usr/bin/env bash
set -euo pipefail

ID_FILE="${1:-ids.txt}"
OUT_DIR="${2:-data/pdb}"
FAILED_FILE="${3:-}"

if [ ! -f "$ID_FILE" ]; then
  echo "ID list not found: $ID_FILE" >&2
  exit 1
fi

if [ -z "$FAILED_FILE" ]; then
  if [ "$ID_FILE" = "ids.txt" ]; then
    FAILED_FILE="data/failed_ids.txt"
  else
    stem=$(basename "$ID_FILE")
    stem="${stem%.*}"
    FAILED_FILE="data/failed_${stem}.txt"
  fi
fi

mkdir -p "$OUT_DIR"
mkdir -p "$(dirname "$FAILED_FILE")"
: > "$FAILED_FILE"

while read -r id; do
  id=$(echo "$id" | tr '[:lower:]' '[:upper:]' | xargs)
  [ -z "$id" ] && continue

  out="${OUT_DIR}/${id}.cif.gz"

  if [ -f "$out" ]; then
    echo "[SKIP] $id already exists"
    continue
  fi

  echo "[DOWNLOADING] $id"
  if ! curl -fL "https://files.rcsb.org/download/${id}.cif.gz" -o "$out"; then
    echo "[FAILED] $id"
    rm -f "$out"
    echo "$id" >> "$FAILED_FILE"
  fi
done < "$ID_FILE"

echo "Done."
echo "Failed IDs saved to $FAILED_FILE"
