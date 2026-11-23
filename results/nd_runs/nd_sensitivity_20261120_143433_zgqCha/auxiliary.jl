####
#### Example script for running in batch mode - calculates sensitivity for one phase and signal
####

global PROJECT_DIR = ENV["SN_PROJECT_DIR"]

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

selected_files = selected_files[.!selected_files.is_signal, :]
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

filtered_data_tables_nd = [
    filter_nd_bounds(table, filter(!=("energy_elec_sum"), var_names_nd), var_bounds_nd)
    for table in selected_tables_nd
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


# Appended to the generated auxiliary script after the shared event selection.
import CairoMakie

roi_edge(x::Real) = Float64(x)
roi_edge(x) = x === :Inf ? Inf : error("Invalid ROI edge: $x")

# CSV stores the printed NamedTuple; read only its literal keys and numbers.
function read_best_roi(text, expected_keys)
    parsed = Meta.parse(text)
    parsed isa Expr && parsed.head === :tuple || error("Expected a tuple in best_ROI")
    roi = Dict{Symbol,Tuple{Float64,Float64}}()
    for entry in parsed.args
        entry isa Expr && entry.head === :(=) && length(entry.args) == 2 ||
            error("Invalid entry in best_ROI: $entry")
        key, interval = entry.args
        key isa Symbol && !haskey(roi, key) &&
            interval isa Expr && interval.head === :tuple && length(interval.args) == 2 ||
            error("Invalid ROI interval: $entry")
        lo, hi = roi_edge.(interval.args)
        0 <= lo < hi || error("Invalid ROI bounds for $key: ($lo, $hi)")
        roi[key] = (lo, hi)
    end
    Set(keys(roi)) == Set(expected_keys) || error("best_ROI variables do not match this run")
    return roi
end

function background_counts(table, roi, base_keys, cuts, energy_prefilter, edges)
    energy_lo, energy_hi = roi[:energy_elec_sum]
    prefilter_lo, prefilter_hi = energy_prefilter
    windows = [(key, roi[key]) for key in base_keys if key !== :energy_elec_sum]
    cut_windows = [(cut.value, roi[cut_column(cut)]) for cut in cuts]
    n_roi = 0
    bins = zeros(Int, length(edges) - 1)

    for event in table
        all(lo <= abs(getproperty(event, key)) < hi for (key, (lo, hi)) in windows) || continue
        all(lo <= abs(value(event)) < hi for (value, (lo, hi)) in cut_windows) || continue

        energy = abs(event.energy_elec_sum)
        if energy_lo <= energy < energy_hi && prefilter_lo <= energy < prefilter_hi
            n_roi += 1
        end
        if first(edges) <= energy < last(edges)
            bins[searchsortedlast(edges, energy)] += 1
        end
    end
    return n_roi, bins
end
function save_background_spectrum(path, edges, spectra, energy_roi, phase)
    fig = CairoMakie.Figure(size = (1100, 600))
    ax = CairoMakie.Axis(
        fig[1, 1];
        title = "Expected background for $phase (full exposure)",
        xlabel = "Electron energy sum (MeV)",
        ylabel = "Expected counts / 50 keV",
        yscale = log10,
    )
    CairoMakie.vspan!(ax, energy_roi...; color = (:gray, 0.25), label = "E_sum ROI")

    names = sort!(collect(keys(spectra)))
    colors = CairoMakie.cgrad(:jet1, length(names); categorical = true)
    total = sum(values(spectra))
    function spectrum_line!(bins; kwargs...)
        # A zero-count bin has no position on a logarithmic axis.
        heights = [n > 0 ? n : NaN for n in bins]
        CairoMakie.stairs!(ax, edges, [heights; last(heights)]; step = :post, kwargs...)
    end
    for (i, name) in enumerate(names)
        spectrum_line!(spectra[name]; color = colors[i], linewidth = 2, label = name)
    end
    spectrum_line!(total; color = :black, linewidth = 3, label = "Total")
    CairoMakie.xlims!(ax, first(edges), last(edges))
    positive = filter(>(0), vcat(collect(values(spectra))...))
    if isempty(positive)
        CairoMakie.ylims!(ax, 1e-5, 1.0)
    else
        CairoMakie.ylims!(ax, minimum(positive) / 2, 2 * maximum(total))
    end
    fig[1, 2] = CairoMakie.Legend(fig, ax, "Processes")
    CairoMakie.save(path, fig; px_per_unit = 2)
end

function save_background_report(output_dir, phase, files, tables, var_names, bounds, cuts)
    result = only(CSV.File(output_dir))
    result.phase == phase || error("results.csv phase does not match the run configuration")
    base_keys = Symbol.(var_names)
    roi = read_best_roi(result.best_ROI, [base_keys; cut_column.(cuts)])
    edges = 0.0:0.05:4.0
    spectra = Dict{String,Vector{Float64}}()
    summary = DataFrame(
        process = String[], mc_events_roi = Int[], mc_events_0_4_mev = Int[],
        expected_roi = Float64[], expected_0_4_mev = Float64[],
    )
    for (row, table) in zip(eachrow(files), tables)
        n_roi, bins = background_counts(table, roi, base_keys, cuts, bounds.energy_elec_sum, edges)
        weight = row.amount * row.activity_bq * SNparams["t"] / row.nb_events
        name = String(row.name)
        push!(summary, (name, n_roi, sum(bins), weight * n_roi, weight * sum(bins)))
        # As in ex5, combine PMTs in the plot but keep each process in the CSV.
        group = occursin("pmt", name) ? "pmt" : name
        get!(() -> zeros(length(bins)), spectra, group) .+= weight .* bins
    end

    expected_roi = sum(summary.expected_roi)
    isapprox(expected_roi, result.best_bkg; atol = 1e-10) ||
        error("Recomputed ROI background ($expected_roi) does not match results.csv ($(result.best_bkg))")
    push!(summary, ("Total", sum(summary.mc_events_roi), sum(summary.mc_events_0_4_mev),
                    expected_roi, sum(summary.expected_0_4_mev)))
    println("\n===== EXPECTED BACKGROUND (FULL EXPOSURE) =====")
    show(stdout, MIME"text/plain"(), summary; allrows = true, allcols = true)
    println()

    run_dir = dirname(output_dir)
    CSV.write(joinpath(run_dir, "background_components.csv"), summary)
    save_background_spectrum(joinpath(run_dir, "background_spectrum.png"), edges,
                             spectra, roi[:energy_elec_sum], phase)
end


save_background_report(output_dir, phase, selected_files, filtered_data_tables_nd,
                       var_names_nd, var_bounds_nd, optimized_cuts)
