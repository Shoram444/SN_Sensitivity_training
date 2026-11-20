#!/usr/bin/env bash
set -euo pipefail

if (( $# != 0 )); then
    echo "Usage: bash scripts/list_root_events.sh > root_events.csv" >&2
    exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
project_dir="$(cd "$script_dir/.." && pwd -P)"
metadata_csv="$project_dir/data/mathis_data_info_final.csv"
if [[ ! -f "$metadata_csv" ]]; then
    echo "Metadata CSV not found: $metadata_csv" >&2
    exit 1
fi
if ! command -v root >/dev/null; then
    echo "ROOT is not on PATH; load the ROOT module first." >&2
    exit 1
fi
if ! command -v python3 >/dev/null; then
    echo "Python 3 is required to read the metadata CSV." >&2
    exit 1
fi

temp_dir="$(mktemp -d)"
cleanup() {
    rm -f "$temp_dir/files" "$temp_dir/events.csv"
    rmdir "$temp_dir"
}
trap cleanup EXIT

cd "$project_dir"
python3 - "$metadata_csv" > "$temp_dir/files" <<'PY'
import csv
import os
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as source:
    reader = csv.DictReader(source)
    if not reader.fieldnames or "file_path" not in reader.fieldnames:
        raise ValueError("Metadata CSV has no file_path column")
    seen = set()
    for row in reader:
        path = row["file_path"]
        if not path or "\0" in path:
            raise ValueError(f"Invalid file_path on CSV line {reader.line_num}")
        if path not in seen:
            seen.add(path)
            sys.stdout.buffer.write(os.fsencode(path) + b"\0")
    if not seen:
        raise ValueError("Metadata CSV contains no file paths")
PY
export ROOT_EVENTS_LIST="$temp_dir/files"
export ROOT_EVENTS_OUTPUT="$temp_dir/events.csv"

if root -l -b -q "$script_dir/list_root_events.C" >&2; then
    status=0
else
    status=$?
fi
if [[ -f "$ROOT_EVENTS_OUTPUT" ]]; then
    cat "$ROOT_EVENTS_OUTPUT"
fi
exit "$status"
