#!/usr/bin/env bash

# Test generation and submission without a SLURM installation or ROOT data.
set -euo pipefail
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
tmp=$(mktemp -d "${TMPDIR:-/tmp}/nd-run2-test.XXXXXX")
tmp=$(cd "$tmp" && pwd -P)
run1=""
run2=""
cleanup() {
    local dir
    for dir in "$run1" "$run2"; do
        if [[ -n "$dir" && "$dir" == "$tmp/runs/"* ]]; then
            rm -f "$dir/config.conf" "$dir/settings.jl" "$dir/analysis.jl" "$dir/job.sh" "$dir/job_id.txt"
            rmdir "$dir"
        fi
    done
    rm -f "$tmp/example.conf" "$tmp/sbatch" "$tmp/sbatch-args" "$tmp/bad.conf"
    rmdir "$tmp/runs" "$tmp"
}
trap cleanup EXIT

sed "s|^# RUN_OUTPUT_ROOT=.*|# RUN_OUTPUT_ROOT=$tmp/runs|" \
    "$project_dir/configs/nd_sensitivity.example.conf" > "$tmp/example.conf"

cat > "$tmp/sbatch" <<'BASH'
#!/usr/bin/env bash
printf '%s\n' "$@" > "$SBATCH_ARGS_FILE"
echo 12345
BASH
chmod +x "$tmp/sbatch"
export PATH="$tmp:$PATH"
export SBATCH_ARGS_FILE="$tmp/sbatch-args"

output=$(bash "$project_dir/scripts/run2.sh" "$tmp/example.conf" --dry-run)
run1=${output%%$'\n'*}
run1=${run1#Run directory: }
[[ -f "$run1/config.conf" && -f "$run1/settings.jl" && -f "$run1/analysis.jl" && -f "$run1/job.sh" ]]
[[ ! -e "$SBATCH_ARGS_FILE" ]]
grep -q 'include(joinpath(@__DIR__, "settings.jl"))' "$run1/analysis.jl"
grep -q 'use_om_cut .& (cut_om_size .| cut_om_range)' "$run1/analysis.jl"
grep -q 'algorithm == "ECA"' "$run1/analysis.jl"
bash -n "$run1/job.sh"

output=$(bash "$project_dir/scripts/run2.sh" "$tmp/example.conf")
run2=${output%%$'\n'*}
run2=${run2#Run directory: }
[[ "$run1" != "$run2" ]]
[[ $(<"$run2/job_id.txt") == 12345 ]]
grep -Fx -- "--chdir=$run2" "$SBATCH_ARGS_FILE" >/dev/null
grep -Fx -- "--output=$run2/slurm-%j.out" "$SBATCH_ARGS_FILE" >/dev/null

sed 's/^# RUN_CPUS=.*/# RUN_CPUS=0/' "$tmp/example.conf" > "$tmp/bad.conf"
if bash "$project_dir/scripts/run2.sh" "$tmp/bad.conf" --dry-run >/dev/null 2>&1; then
    echo "Invalid CPU count was accepted" >&2
    exit 1
fi

echo "run2 generation and mocked submission tests passed"
