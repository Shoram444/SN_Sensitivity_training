#!/usr/bin/env bash
set -euo pipefail

default_dir=/sps/nemo/scratch/granjon/double_beta/cimrman_curve/ptd_pol2/root_with_max_threshold_working
if (( $# > 2 )); then
    echo "Usage: bash scripts/list_run_times.sh [DIRECTORY [NAME_PATTERN]] > run_times.csv" >&2
    exit 2
fi
data_dir="${1:-$default_dir}"
pattern="${2:-*.root}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ ! -d "$data_dir" ]]; then
    echo "Directory not found: $data_dir" >&2
    exit 1
fi
if ! command -v root >/dev/null; then
    echo "ROOT is not on PATH; load the ROOT module first." >&2
    exit 1
fi

temp_dir="$(mktemp -d)"
cleanup() {
    rm -f "$temp_dir/files" "$temp_dir/run_times.csv"
    rmdir "$temp_dir"
}
trap cleanup EXIT

find -L "$data_dir" -type f -name "$pattern" -print0 | sort -z > "$temp_dir/files"
if [[ ! -s "$temp_dir/files" ]]; then
    echo "No files matching '$pattern' in $data_dir" >&2
    exit 1
fi
export RUN_TIMES_LIST="$temp_dir/files"
export RUN_TIMES_OUTPUT="$temp_dir/run_times.csv"

if root -l -b -q "$script_dir/list_run_times.C" >&2; then
    status=0
else
    status=$?
fi
if [[ -f "$RUN_TIMES_OUTPUT" ]]; then
    cat "$RUN_TIMES_OUTPUT"
fi
exit "$status"
