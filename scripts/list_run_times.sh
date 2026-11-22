#!/usr/bin/env bash
set -euo pipefail

default_dir=/sps/nemo/scratch/granjon/double_beta/cimrman_curve/ptd_pol2/root_with_max_threshold_working
if (( $# > 1 )); then
    echo "Usage: bash scripts/list_run_times.sh [DIRECTORY] > run_times.csv" >&2
    exit 2
fi
data_dir="${1:-$default_dir}"
data_dir="${data_dir%/}"

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

# Only DIRECTORY/run_XXX/extracted_data.root is read.
find -L "$data_dir" -mindepth 2 -maxdepth 2 -type f -path "$data_dir/run_*/extracted_data.root" -print0 \
    | sort -z > "$temp_dir/files"
if [[ ! -s "$temp_dir/files" ]]; then
    echo "No run_*/extracted_data.root files in $data_dir" >&2
    exit 1
fi
for run_dir in "$data_dir"/run_*/; do
    [[ -e "$run_dir/extracted_data.root" ]] || echo "Skipping $run_dir: no extracted_data.root" >&2
done
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
