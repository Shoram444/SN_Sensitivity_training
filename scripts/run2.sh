#!/usr/bin/env bash

# Prepare one reproducible run and submit it to SLURM.
# Usage: bash scripts/run2.sh CONFIG [--dry-run]
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 || (${2:-} != "" && ${2:-} != "--dry-run") ]]; then
    echo "Usage: bash $0 CONFIG [--dry-run]" >&2
    exit 2
fi

config=$1
dry_run=${2:-}
[[ -f "$config" ]] || { echo "Configuration not found: $config" >&2; exit 1; }
config_dir=$(cd "$(dirname "$config")" && pwd -P)
config="$config_dir/$(basename "$config")"
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
template="$project_dir/analysis/run_nd_sensitivity_on_mathis_data_v2.jl"
[[ -f "$template" ]] || { echo "Julia template not found: $template" >&2; exit 1; }

# RUN_ settings are comments, not shell code. The Julia section is executable
# Julia code, so only use configuration files you trust.
read_setting() {
    local key=$1 count
    count=$(grep -c "^# ${key}=" "$config" || true)
    [[ "$count" == 1 ]] || { echo "Expected exactly one # ${key}= line" >&2; return 1; }
    sed -n "s/^# ${key}=//p" "$config"
}

output_root=$(read_setting RUN_OUTPUT_ROOT)
job_name=$(read_setting RUN_JOB_NAME)
partition=$(read_setting RUN_PARTITION)
licenses=$(read_setting RUN_LICENSES)
walltime=$(read_setting RUN_WALLTIME)
memory=$(read_setting RUN_MEMORY)
cpus=$(read_setting RUN_CPUS)
julia_module=$(read_setting RUN_MODULE)
sysimage=$(read_setting RUN_SYSIMAGE)
mail_user=$(read_setting RUN_MAIL_USER)
save_auxiliary_data=false
if grep -q '^# RUN_SAVE_AUXILIARY_DATA=' "$config"; then
    save_auxiliary_data=$(read_setting RUN_SAVE_AUXILIARY_DATA)
fi
[[ "$save_auxiliary_data" == true || "$save_auxiliary_data" == false ]] ||
    { echo "RUN_SAVE_AUXILIARY_DATA must be true or false" >&2; exit 1; }

[[ -n "$output_root" && -n "$partition" && -n "$walltime" && -n "$memory" ]] ||
    { echo "Output root, partition, walltime and memory must not be empty" >&2; exit 1; }
[[ "$job_name" =~ ^[a-zA-Z0-9_.-]+$ ]] ||
    { echo "RUN_JOB_NAME must contain only letters, digits, _, - or ." >&2; exit 1; }
[[ "$cpus" =~ ^[1-9][0-9]*$ ]] ||
    { echo "RUN_CPUS must be a positive integer" >&2; exit 1; }

awk '
    /^# BEGIN JULIA CONFIG$/ { starts++; active = 1; next }
    /^# END JULIA CONFIG$/ { ends++; if (!active) exit 1; active = 0 }
    END { if (starts != 1 || ends != 1 || active) exit 1 }
' "$config" || { echo "Expected one ordered pair of Julia configuration markers" >&2; exit 1; }

# Refuse to generate a partial script if the template has changed.
check_anchor() {
    local pattern=$1 count
    count=$(grep -Ec "$pattern" "$template" || true)
    [[ "$count" == 1 ]] || { echo "Expected one template anchor: $pattern" >&2; return 1; }
}
check_anchor '^global PROJECT_DIR = '
check_anchor '^global signal_name = '
check_anchor '^global α = '
check_anchor '^signal_rows = '
check_anchor '^selected_tables_nd = '
check_anchor '^        if !cut_om_size\[j\] && !cut_om_range\[j\] && '
check_anchor '^    keep = \.!\(cut_same_calo '
check_anchor '^var_bounds_nd = \($'
check_anchor '^var_steps_nd = \($'
check_anchor '^nd_numeric_cuts = \[$'
check_anchor '^cut_steps\(cut\) = '
check_anchor '^algo = PSO\(; options\)$'
check_anchor '^initial_roi = \($'
check_anchor '^x0 = float\.\(encode_grid_roi\(initial_roi, grid\)\)$'

case "$output_root" in /*) ;; *) output_root="$project_dir/$output_root" ;; esac
case "$sysimage" in ""|/*) ;; *) sysimage="$project_dir/$sysimage" ;; esac
mkdir -p "$output_root"
run_dir=$(mktemp -d "$output_root/${job_name}_$(date +%Y%m%d_%H%M%S)_XXXXXX")
run_dir=$(cd "$run_dir" && pwd -P)
cp "$config" "$run_dir/config.conf"
sed -n '/^# BEGIN JULIA CONFIG$/,/^# END JULIA CONFIG$/{ /^# BEGIN JULIA CONFIG$/d; /^# END JULIA CONFIG$/d; p; }' \
    "$run_dir/config.conf" > "$run_dir/settings.jl"

# Insert background selection immediately after the signal lookup. This keeps
# the selected signal and rejects unknown background names before ROOT loading.
cat > "$run_dir/select_backgrounds.jl" <<'JULIA'
if background_names !== :all
    available = Set(String.(selected_files.name[.!selected_files.is_signal]))
    requested = Set(String.(background_names))
    unknown = setdiff(requested, available)
    isempty(unknown) || error("Unknown backgrounds for $phase: $unknown")
    selected_files = selected_files[
        selected_files.is_signal .| in.(selected_files.name, Ref(requested)), :
    ]
end
JULIA

cat > "$run_dir/validate_initial.jl" <<'JULIA'
for name_string in var_names_nd
    name = Symbol(name_string)
    minimum_edge, maximum_edge = var_bounds_nd[name]
    step = var_steps_nd[name]
    step > 0 && minimum_edge < maximum_edge || error("Invalid ROI grid for $name")
    grid_size = (maximum_edge - minimum_edge) / step
    isapprox(grid_size, round(grid_size); atol=1e-8) ||
        error("ROI bounds must align with the step for $name")
    lo, hi = initial_roi[name]
    minimum_edge <= lo < hi <= maximum_edge ||
        error("Initial ROI for $name is outside its bounds or empty")
    for edge in (lo, hi)
        index = (edge - minimum_edge) / step
        isapprox(index, round(index); atol=1e-8) ||
            error("Initial ROI edge $edge for $name is not on the grid")
    end
end
for cut in optimized_cuts
    lo, hi = cut.bounds
    lo <= cut.fixed <= hi || error("Initial threshold for $(cut.name) is outside its bounds")
    index = (cut.fixed - lo) / cut.step
    isapprox(index, round(index); atol=1e-8) ||
        error("Initial threshold for $(cut.name) is not on the grid")
end
JULIA

# Run sed in the new directory so its r commands use simple relative paths.
(
    cd "$run_dir"
    sed \
        -e 's|^global PROJECT_DIR = .*|global PROJECT_DIR = ENV["SN_PROJECT_DIR"]|' \
        -e '/^global signal_name = /,/^global α = /c\
include(joinpath(@__DIR__, "settings.jl"))
' \
        -e '/^@info "Running sensitivity optimization/i\
algorithm in ("PSO", "ECA") || error("Unsupported algorithm: $algorithm")\
optimization_time_limit > 0 || error("Optimization time limit must be positive")
' \
        -e '/^signal_rows = /r select_backgrounds.jl' \
        -e '/^var_bounds_nd = ($/,/^)/d' \
        -e '/^var_steps_nd = ($/,/^)/d' \
        -e '/^nd_numeric_cuts = \[$/,/^]/d' \
        -e '/^cut_steps(cut) = /r validate_initial.jl' \
        -e 's/^        if !cut_om_size\[j\] \&\& !cut_om_range\[j\] \&\& /        if /' \
        -e '/^    keep = \.!(cut_same_calo /c\
    keep = .!((use_same_calo_cut .& cut_same_calo) .| (use_kinks_cut .& cut_kinks) .| (use_om_cut .& (cut_om_size .| cut_om_range)) .| (use_vertex_region_cut .& cut_vertex_region))
' \
        -e '/^selected_tables_nd = /i\
append!(var_names, extra_var_names)\
unique!(var_names)
' \
        -e 's/^algo = PSO(; options)$/algo = algorithm == "PSO" ? PSO(; options) : algorithm == "ECA" ? ECA(; options) : error("Unsupported algorithm: $algorithm")/' \
        -e 's/Starting optimization with PSO algorithm/Starting optimization with $algorithm algorithm/' \
        -e '/^initial_roi = ($/,/^)/d' \
        "$template" > analysis.jl
)
grep -q '^CSV.write(output_dir, df)$' "$run_dir/analysis.jl" ||
    { echo "Generated Julia script is incomplete; check template changes" >&2; exit 1; }
rm "$run_dir/select_backgrounds.jl" "$run_dir/validate_initial.jl"

if [[ "$save_auxiliary_data" == true ]]; then
    cp "$project_dir/analysis/nd_background_report.jl" "$run_dir/auxiliary.jl"
fi

# SLURM options are passed as arguments, so paths with spaces remain intact.
cat > "$run_dir/job.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
BASH
{
    printf 'export SN_PROJECT_DIR=%q\n' "$project_dir"
    printf 'export SN_OUTPUT_PATH=%q\n' "$run_dir/results.csv"
    printf 'export SN_SAVE_AUXILIARY_DATA=%q\n' "$save_auxiliary_data"
    printf 'export JULIA_NUM_THREADS=%q\n' "$cpus"
    printf '%s\n' 'export JULIA_CPU_TARGET=generic'
    printf 'cd %q\n' "$project_dir"
    if [[ -n "$julia_module" ]]; then
        printf 'module load %q\n' "$julia_module"
    fi
    julia_cmd=(julia "--project=$project_dir" -t "$cpus" --pkgimages=no)
    if [[ -n "$sysimage" ]]; then
        julia_cmd+=(-J "$sysimage")
    fi
    printf 'exec '
    printf '%q ' "${julia_cmd[@]}"
    printf '%q\n' "$run_dir/analysis.jl"
} >> "$run_dir/job.sh"

sbatch_args=(
    "--job-name=$job_name"
    "--partition=$partition"
    "--time=$walltime"
    "--mem=$memory"
    "--cpus-per-task=$cpus"
    "--chdir=$run_dir"
    "--output=$run_dir/slurm-%j.out"
    "--error=$run_dir/slurm-%j.err"
)
[[ -z "$licenses" ]] || sbatch_args+=("--licenses=$licenses")
if [[ -n "$mail_user" ]]; then
    sbatch_args+=("--mail-type=END,FAIL" "--mail-user=$mail_user")
fi

echo "Run directory: $run_dir"
if [[ "$dry_run" == "--dry-run" ]]; then
    printf 'Would submit: sbatch'
    printf ' %q' "${sbatch_args[@]}" "$run_dir/job.sh"
    printf '\n'
else
    job_id=$(sbatch --parsable "${sbatch_args[@]}" "$run_dir/job.sh") ||
        { echo "SLURM submission failed; prepared files remain in $run_dir" >&2; exit 1; }
    printf '%s\n' "$job_id" > "$run_dir/job_id.txt"
    echo "Submitted SLURM job: $job_id"
fi
