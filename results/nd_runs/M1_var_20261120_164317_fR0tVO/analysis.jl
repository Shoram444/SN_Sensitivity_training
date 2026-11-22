####
#### Example script for running in batch mode - calculates sensitivity for one phase and signal
####

global PROJECT_DIR = ENV["SN_PROJECT_DIR"]
const save_auxiliary_data = get(ENV, "SN_SAVE_AUXILIARY_DATA", "false") == "true"

# Settings can also be exported from run.sh before starting Julia.
include(joinpath(@__DIR__, "settings.jl"))

algorithm in ("PSO", "ECA") || error("Unsupported algorithm: $algorithm")
optimization_time_limit > 0 || error("Optimization time limit must be positive")
@info "Running sensitivity optimization for phase: $phase, with signal: $signal_name, output will be saved to: $output_dir"
flush(stdout)

using SNSensitivityEstimate, UnROOT, FHist, CSV, DataFramesMeta
using Distributions, Metaheuristics

@info "Packages loaded, loading data..."
flush(stdout)

# STEP 1: Load the data
data_info = CSV.read(joinpath(PROJECT_DIR, "data", "mathis_data_info_final_2.csv"), DataFrame)

# The list contains bad OMs; the first line is a header.
om_list_files = Dict(
    "phase_0" => "list_om_phase_0.txt",
    "phase_1" => "list_om_phase_1.txt",
    "phase_2" => "list_om_phase_2a.txt",
    "phase_3" => "list_om_phase_2b.txt",
    "phase_4" => "list_om_phase_3.txt",
)

function read_bad_oms(path)
    lines = readlines(path)
    isempty(lines) && error("OM list is empty: $path")
    bad = Set{Int}()
    for line in lines[2:end], token in split(line)
        om = tryparse(Int, token)
        om === nothing && continue
        0 <= om <= 711 || error("Invalid OM identifier $om in $path")
        push!(bad, om)
    end
    return bad
end

haskey(om_list_files, phase) || error("No OM list configured for $phase")
om_list_path = joinpath(PROJECT_DIR, "data", "om_lists", om_list_files[phase])
bad_oms = read_bad_oms(om_list_path)
@info "Loaded $(length(bad_oms)) bad OMs for $phase from $om_list_path"

# Only load the selected signal and backgrounds for this phase.
selected_files = @chain data_info begin
    @subset :phase .== phase
    @subset .!occursin.("reduces", :file_path)
    @subset .! :is_signal .| (:name .== signal_name)
end
signal_rows = findall(selected_files.is_signal .& (selected_files.name .== signal_name))
if background_names !== :all
    available = Set(String.(selected_files.name[.!selected_files.is_signal]))
    requested = Set(String.(background_names))
    missing = setdiff(requested, available)
    isempty(missing) || error("Unknown backgrounds for $phase: $missing")
    selected_files = selected_files[
        selected_files.is_signal .| in.(selected_files.name, Ref(requested)), :
    ]
end
length(signal_rows) == 1 || error("Expected exactly one $signal_name signal in $phase")
any(.!selected_files.is_signal) || error("No backgrounds found for $phase")

data_files = [UnROOT.ROOTFile(joinpath(PROJECT_DIR, file_path)) for file_path in selected_files.file_path]
tree_name = "Event"

# We don't need all the variables, only those used by the analysis and cuts.
var_names = [
    "energy_elec_sum",
    "angle_3D_between_ep_em",
    "delta_y_elec",
    "delta_z_elec",
    "diff_time_elec",
    "energy_elec_1",
    "energy_elec_2",
    "hit_the_same_calo_hit",
    "closest_gamma",
    "closest_elec",
    "closest_track",
    "closest_time_track",
    "has_kinks",
    "num_om_elec_f",
    "vertex_3D_start_y",
    "vertex_3D_start_z",
]

# Sensitivity ND
var_names_nd = [
    "energy_elec_sum",
    "angle_3D_between_ep_em",
    "delta_y_elec",
    "delta_z_elec",
    "diff_time_elec",
    "energy_elec_1",
    "energy_elec_2",
]

# The ND package applies these bounds to abs(value), including diff_time_elec.

# Default search resolution (adjust if needed).

# Fixed cuts use comparison and fixed. To optimize a cut, change mode to
# :optimize and set bounds and step for its threshold.
# The ND package applies optimized cuts to abs(value): < max or >= min.
# Include any new variables used by value in var_names above.

for cut in nd_numeric_cuts
    cut.mode in (:fixed, :optimize) || error("Unknown mode for $(cut.name)")
    cut.comparison in ((<), (<=), (>), (>=)) || error("Unknown comparison for $(cut.name)")
    if cut.mode == :optimize
        lo, hi = cut.bounds
        cut.step > 0 && 0 <= lo <= hi || error("Invalid threshold grid for $(cut.name)")
        isapprox((hi - lo) / cut.step, round((hi - lo) / cut.step)) ||
            error("Threshold bounds must align with the step for $(cut.name)")
    end
end
length(unique(cut.name for cut in nd_numeric_cuts)) == length(nd_numeric_cuts) ||
    error("ND numeric cut names must be unique")
optimized_cuts = [cut for cut in nd_numeric_cuts if cut.mode == :optimize]
cut_column(cut) = Symbol("__nd_cut_", cut.name)
cut_steps(cut) = round(Int, (cut.bounds[2] - cut.bounds[1]) / cut.step)
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

append!(var_names, extra_var_names)
unique!(var_names)
selected_tables_nd = [UnROOT.LazyTree(data_file, tree_name, var_names) for data_file in data_files]

# Match the package's ND ROI convention: min <= abs(value) < max.
function filter_nd_bounds(table, var_names, bounds)
    keep = trues(length(table))
    for name in var_names
        lo, hi = bounds[Symbol(name)]
        values = getproperty(table, Symbol(name))
        keep .&= (lo .<= abs.(values)) .& (abs.(values) .< hi)
    end
    return table[keep]
end

const background_spectrum_edges = 0.0:0.05:4.0
energy_bounds = var_bounds_nd.energy_elec_sum
report_bounds = merge(var_bounds_nd, (
    energy_elec_sum = (min(first(energy_bounds), first(background_spectrum_edges)),
                       max(last(energy_bounds), last(background_spectrum_edges))),
))
filtered_data_tables_nd = [
    filter_nd_bounds(table, var_names_nd,
                     save_auxiliary_data && !is_signal ? report_bounds : var_bounds_nd)
    for (table, is_signal) in zip(selected_tables_nd, selected_files.is_signal)
]

@info "Data loaded and filtered, applying Mathis's extra cuts..."
flush(stdout)

# Apply fixed topological cuts before calculating numeric cut variables.
for i in eachindex(filtered_data_tables_nd)
    table = filtered_data_tables_nd[i]
    cut_same_calo = table.hit_the_same_calo_hit .== 1
    cut_kinks = table.has_kinks .== 1

    n_evt = length(table)
    cut_om_size = [length(table.num_om_elec_f[j]) != 2 for j in 1:n_evt]
    cut_om_range = falses(n_evt)
    cut_vertex_region = falses(n_evt)

    for j in 1:n_evt
        if !cut_om_size[j]
            om0, om1 = table.num_om_elec_f[j]
            cut_om_range[j] = !(0 <= om0 < 520 && 0 <= om1 < 520) ||
                              om0 in bad_oms || om1 in bad_oms
        end

        if length(table.vertex_3D_start_y[j]) >= 2 && length(table.vertex_3D_start_z[j]) >= 2
            y_mid = 0.5 * (table.vertex_3D_start_y[j][1] + table.vertex_3D_start_y[j][2])
            z_mid = 0.5 * (table.vertex_3D_start_z[j][1] + table.vertex_3D_start_z[j][2])

            cut_vertex_region[j] = (
                ((y_mid - 1180) / 60)^2 + ((z_mid + 850) / 140)^2 <= 1 ||
                ((y_mid - 1144) / 25)^2 + ((z_mid - 196) / 40)^2 <= 1 ||
                ((y_mid - 1114) / 25)^2 + ((z_mid + 257) / 40)^2 <= 1 ||
                abs(y_mid) > 2360 || abs(z_mid) > 1340
            )
        end
    end

    keep = .!((use_same_calo_cut .& cut_same_calo) .| (use_kinks_cut .& cut_kinks) .| (use_om_cut .& (cut_om_size .| cut_om_range)) .| (use_vertex_region_cut .& cut_vertex_region))
    filtered_data_tables_nd[i] = table[keep]
end

# Optimized numeric cuts are left unfiltered until the ROI is evaluated.
function apply_fixed_numeric_cuts(table, cuts)
    keep = trues(length(table))
    for cut in cuts
        if cut.mode == :fixed
            values = map(cut.value, table)
            keep .&= cut.comparison.(values, cut.fixed)
        end
    end
    return table[keep]
end

filtered_data_tables_nd = [
    apply_fixed_numeric_cuts(table, nd_numeric_cuts)
    for table in filtered_data_tables_nd
]

# Cache only the ROI columns, including derived cuts, before narrowing E_sum
# back to the optimization range. Both the report and PSO reuse these values.
function cache_background_table(table, names, cuts)
    isempty(table) && return table
    columns = NamedTuple(Symbol(name) => collect(getproperty(table, Symbol(name))) for name in names)
    derived = NamedTuple(cut_column(cut) => map(cut.value, table) for cut in cuts)
    return UnROOT.LazyTree(merge(columns, derived))
end

if save_auxiliary_data
    background_rows = findall(!, selected_files.is_signal)
    report_files = selected_files[background_rows, :]
    report_tables = [
        cache_background_table(filtered_data_tables_nd[i], var_names_nd, optimized_cuts)
        for i in background_rows
    ]
    report_by_row = Dict(zip(background_rows, report_tables))
    filtered_data_tables_nd = [
        haskey(report_by_row, i) ?
            filter_nd_bounds(report_by_row[i], ("energy_elec_sum",), var_bounds_nd) :
            table
        for (i, table) in enumerate(filtered_data_tables_nd)
    ]
end

# Empty backgrounds contribute no events. with_variable requires nonempty data.
nd_nonempty = findall(i -> !isempty(filtered_data_tables_nd[i]), eachindex(filtered_data_tables_nd))
for i in eachindex(filtered_data_tables_nd)
    if isempty(filtered_data_tables_nd[i])
        selected_files.is_signal[i] && error("No signal events remain after ND cuts")
        @warn "Skipping empty ND background" name=selected_files.name[i]
    end
end
selected_files = selected_files[nd_nonempty, :]
filtered_data_tables_nd = filtered_data_tables_nd[nd_nonempty]

@info "Mathis's extra cuts applied, creating ND processes..."
flush(stdout)

# STEP 2: Create DataProcessND objects.
# Keep the original batch script's normalization to the full exposure.
cached_names = [var_names_nd; String.(cut_column.(optimized_cuts))]
cached_bounds = merge(var_bounds_nd, NamedTuple(cut_column(cut) => (0.0, Inf) for cut in optimized_cuts))
processes_nd = [
    SNSensitivityEstimate.DataProcessND(
        filtered_data_tables_nd[i],
        String(selected_files.name[i]),
        selected_files.is_signal[i],
        selected_files.activity_bq[i],
        SNparams["t"],
        selected_files.nb_events[i],
        save_auxiliary_data && !selected_files.is_signal[i] ? cached_bounds : var_bounds_nd,
        selected_files.amount[i],
        save_auxiliary_data && !selected_files.is_signal[i] ? cached_names : var_names_nd,
    ) for i in eachindex(filtered_data_tables_nd)
]

processes_nd = let processes = processes_nd
    for cut in optimized_cuts
        if cut.comparison === (>) || cut.comparison === (<=)
            @warn "Optimized cut uses the package's >= or < interval boundary" name=cut.name
        end
        processes = [
            save_auxiliary_data && !process.signal ? process :
                with_variable(process, cut_column(cut), cut.value, (0.0, Inf))
            for process in processes
        ]
    end
    processes
end

signal_matches = get_process(signal_name, processes_nd)
length(signal_matches) == 1 && only(signal_matches).signal ||
    error("Expected exactly one $signal_name signal in $phase")
count(p -> p.signal, processes_nd) == 1 || error("ND optimization requires exactly one signal")
any(p -> !p.signal, processes_nd) || error("No background events remain after ND cuts")

signal_process = only(signal_matches)
const PROCESSES_ND = processes_nd

# STEP 3: Define the optimization problem.
@info "Defining the optimization problem..."
flush(stdout)

# Base ROI edges have two indices each; optimized cuts add one threshold each.
for name in var_names_nd
    lo, hi = var_bounds_nd[Symbol(name)]
    step = var_steps_nd[Symbol(name)]
    step > 0 && lo < hi || error("Invalid ROI grid for $name")
    isapprox((hi - lo) / step, round((hi - lo) / step)) ||
        error("ROI bounds must align with the step for $name")
end
grid_mins = [Float64(var_bounds_nd[Symbol(name)][1]) for name in var_names_nd]
grid_steps = [Float64(var_steps_nd[Symbol(name)]) for name in var_names_nd]
grid_sizes = [
    round(Int, (var_bounds_nd[Symbol(name)][2] - var_bounds_nd[Symbol(name)][1]) /
               var_steps_nd[Symbol(name)])
    for name in var_names_nd
]
grid = GridSpec(Symbol.(var_names_nd), grid_mins, grid_steps, grid_sizes)

searchRange = make_stepRange(grid)
for cut in optimized_cuts
    push!(searchRange, (0, cut_steps(cut)))
end

function nd_thresholds(indices)
    offset = 2 * length(var_names_nd)
    return NamedTuple(
        cut.name => min(cut.bounds[2], cut.bounds[1] + round(Int, indices[offset + i]) * cut.step)
        for (i, cut) in enumerate(optimized_cuts)
    )
end

function nd_roi(indices)
    edges = decode_grid_roi(indices[1:2 * length(var_names_nd)], grid)
    roi = NamedTuple(
        Symbol(name) => (edges[2i-1], edges[2i])
        for (i, name) in enumerate(var_names_nd)
    )
    thresholds = nd_thresholds(indices)
    cut_intervals = NamedTuple(
        cut_column(cut) => (
            cut.comparison in ((<), (<=)) ?
            (0.0, thresholds[cut.name]) :
            (thresholds[cut.name], Inf)
        )
        for cut in optimized_cuts
    )
    return merge(roi, cut_intervals)
end

# The formula approximation is faster; the final result uses the table method.
function prob(proposed_indices)
    roi = nd_roi(proposed_indices)
    edges = Float64[edge for interval in values(roi) for edge in interval]
    return -get_s_to_b(SNparams, α, PROCESSES_ND, edges; approximate="formula")
end

lower_bound = [x[1] for x in searchRange] .|> float
upper_bound = [x[2] for x in searchRange] .|> float
options = Options(;
    time_limit = float(optimization_time_limit),
    verbose = true,
)
bounds = boxconstraints(lb = lower_bound, ub = upper_bound)
algo = algorithm == "PSO" ? PSO(; options) : algorithm == "ECA" ? ECA(; options) : error("Unsupported algorithm: $algorithm")

# Initial ROI from the original batch script, with a symmetric time cut.
x0 = float.(encode_grid_roi(initial_roi, grid))
for cut in optimized_cuts
    initial = clamp(round(Int, (cut.fixed - cut.bounds[1]) / cut.step), 0, cut_steps(cut))
    push!(x0, initial)
end
set_user_solutions!(algo, x0, prob)

@info "Starting optimization with $algorithm algorithm..."
flush(stdout)

result = optimize(prob, bounds, algo)

@info "Optimization finished, processing results..."
flush(stdout)

@show minimum(result)
@show res = minimizer(result)

best = nd_roi(res)
best_thresholds = nd_thresholds(res)
best_sens = get_sensitivityND(SNparams, α, PROCESSES_ND, best; approximate="table")
best_sens.signalEff > 0 || error("ND optimization found no signal events")

println("\n===== BEST RESULT =====")
println("Signal process = ", signal_name)
println("Optimized numeric cuts: ", best_thresholds)
print(best_sens)

# Save results to csv.
mkpath(dirname(output_dir))
df = DataFrame(
    phase = phase,
    best_sensitivity = [best_sens.tHalf],
    best_ROI = [best_sens.roi],
    best_bkg = [best_sens.bkgCounts],
    best_eff = [best_sens.signalEff],
)
CSV.write(output_dir, df)

if save_auxiliary_data
    include(joinpath(@__DIR__, "auxiliary.jl"))
    save_background_report(output_dir, phase, signal_process, report_files, report_tables, best,
                           energy_bounds, background_spectrum_edges, best_sens.bkgCounts)
end
