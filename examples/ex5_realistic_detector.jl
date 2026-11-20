#### 
#### WORK IN PROGRESS
####
import Pkg;
# Pkg.activate("/sps/nemo/scratch/mpetro/Sensitivity_training/")
using SNSensitivityEstimate, UnROOT, FHist, CairoMakie, CSV, DataFramesMeta
using Random, LinearAlgebra, Statistics, Distributions, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, SpecialFunctions, ValueShapes
using Distributions, ColorSchemes
using Metaheuristics



# STEP 1: Load the data

# Start is the same as in ex1_background.jl, we load data and create DataProcess objects.
data_info = CSV.read("data/mathis_data_info_final_2.csv", DataFrame)

# Select which phase we process
phase = "phase_4"
signal_name = "bb0nuM1"

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
om_list_path = joinpath(@__DIR__, "..", "data", "om_lists", om_list_files[phase])
bad_oms = read_bad_oms(om_list_path)
@info "Loaded $(length(bad_oms)) bad OMs for $phase from $om_list_path"

selected_files = @chain data_info begin
    @subset :phase .== phase
    @subset .!occursin.("reduces", :file_path) 
end
# manually add nu0bb path because that's not in reduced files
# append!(selected_files, data_info[data_info.file_path .== "data/mathis_files/bb0nu/phase_0/simu_0/rec_bb0nu_0_0.root", :])


data_files = [UnROOT.ROOTFile(file_path) for file_path in selected_files.file_path] 
tree_name = "Event" 

# we don't need all the variables, we can select a few to use for the analysis.
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
    "vertex_3D_start_z"
    ]

data_tables = [UnROOT.LazyTree(data_file, tree_name, var_names) for data_file in data_files] 


# define the boundaries of the data-cuts, we use simple Esum >300keV, dy < 100mm, dz < 100mm, Pint > 0.01, Pext < 0.01
var_bounds = (
    energy_elec_sum = (0.5, 3.0),
    angle_3D_between_ep_em = (0, 180),
    delta_y_elec = (0, 60),
    delta_z_elec = (0, 65),
    diff_time_elec = (-1.75, 1.75),
    energy_elec_1 = (0.35, 2.4),
    energy_elec_2 = (0.35, 2.4),
)

# Apply bounds only for keys present in `var_bounds`.
filtered_data_tables = [SNSensitivityEstimate.filter_data(table, var_names, var_bounds) for table in data_tables]


# Apply additional topological cuts and store filtered tables back.
for i in eachindex(filtered_data_tables)
    table = filtered_data_tables[i]
    cut_vertex_dist = sqrt.(table.delta_y_elec .* table.delta_y_elec .+ table.delta_z_elec .* table.delta_z_elec) .>= 100
    cut_same_calo   = table.hit_the_same_calo_hit .== 1
    cut_gam         = abs.(table.closest_gamma) .<= 50
    cut_elec        = abs.(table.closest_elec) .<= 25
    cut_track       = abs.(table.closest_track) .< 200
    cut_ttrack      = abs.(table.closest_time_track) .< 10
    cut_kinks       = table.has_kinks .== 1

    n_evt = length(table)
    cut_om_size = [length(table.num_om_elec_f[j]) != 2 for j in 1:n_evt]

    cut_om_range = falses(n_evt)
    for j in 1:n_evt
        if !cut_om_size[j]
            om0 = table.num_om_elec_f[j][1]
            om1 = table.num_om_elec_f[j][2]
            cut_om_range[j] = !(0 <= om0 < 520 && 0 <= om1 < 520) ||
                              om0 in bad_oms || om1 in bad_oms
        end
    end

    cut_vertex_region = falses(n_evt)
    for j in 1:n_evt
        if !cut_om_size[j] && !cut_om_range[j] && length(table.vertex_3D_start_y[j]) >= 2 && length(table.vertex_3D_start_z[j]) >= 2
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

    keep = .!(cut_vertex_dist .| cut_same_calo .| cut_gam .| cut_elec .| cut_track .| cut_ttrack .| cut_kinks .| cut_om_size .| cut_om_range .| cut_vertex_region)
    filtered_data_tables[i] = table[keep]
end


processes = [
    SNSensitivityEstimate.DataProcess(
        collect(filtered_data_tables[i].energy_elec_sum), # Here we specify a single variable to be used. In later stages we will use N-dimensional data vectors.
        String(selected_files.name[i]), 
        selected_files.is_signal[i], 
        selected_files.activity_bq[i], 
        selected_files.time_s[i], 
        selected_files.nb_events[i], 
        0:0.1:4, 
        selected_files.amount[i], 
    ) for i in 1:length(filtered_data_tables)
]



mutable struct HistogramData
    hist::FHist.Hist1D
    name::String
    is_signal::Bool
end

function HistogramData(process::SNSensitivityEstimate.DataProcess)
    hist = get_bkg_counts_1D(process)
    HistogramData(hist, process.isotopeName, process.signal)
end

all_bkg_histograms = [HistogramData(process) for process in processes if !process.signal]

function merge_hist_by_keyword(histograms::Vector{HistogramData}, keyword::String, merge_name::String)
    merged_hist = FHist.Hist1D(;binedges = histograms[1].hist.binedges)
    is_signal = findfirst(hist -> occursin(keyword, hist.name) && hist.is_signal, histograms) !== nothing
    for hist in histograms
        if occursin(keyword, hist.name)
            # if signal changes within keyword error out
            if hist.is_signal != is_signal
                error("Can only merge histograms with the same signal flag for keyword $keyword")
            end
            merge!(merged_hist, hist.hist)
        end

    end
    return HistogramData(merged_hist, merge_name, is_signal)
end

pmt_hists = merge_hist_by_keyword(all_bkg_histograms, "pmt", "pmt")
other_histos = filter(hist -> !occursin("pmt", hist.name), all_bkg_histograms)
all_histos = vcat(other_histos, [pmt_hists])

let
    histos = [hist.hist for hist in all_histos]
    colors = cgrad(:jet1, length(histos); categorical=true)
    for i in 1:length(histos)
        println("Process ", all_histos[i].name)
        println("Integral: ", integral(histos[i]))
    end

    f = CairoMakie.Figure(size = (1000, 600))
    a = CairoMakie.Axis(
        f[1, 1], 
        title = "Histogram of Expected Background Counts, phase $phase data\nExposure of 17.5 kg.yr",
        xlabel = "sum of energies", 
        ylabel = "expected counts",
        yscale = log10
        )
    for i in 1:length(histos)
        CairoMakie.stephist!(a, histos[i], color = colors[i], linewidth = 5)
    end
    CairoMakie.stephist!(a, sum(histos), color = :black, linewidth = 5)

    labels = [hist.name for hist in all_histos]
    push!(labels, "total")
    elements = [CairoMakie.LineElement(color = colors[i], linewidth = 5) for i in 1:length(labels)-1]
    push!(elements, CairoMakie.LineElement(color = :black, linewidth = 5))
    title = "Processes"
    
    ylims!(a, 1e-5, nothing)
    CairoMakie.Legend(f[1,2], elements, labels, title)
    save("test_M1_cuts_background_stacked_background_histogram.png", f, px_per_unit = 2)
    # save("data/out/test_M1_cuts_background_stacked_background_histogram.png", f)
    f
end

let
    adjusted_per_day_histos = [hist.hist .* (1/2.86/365) for hist in all_histos]
    colors = cgrad(:jet1, length(adjusted_per_day_histos); categorical=true)
    f = CairoMakie.Figure(size = (1000, 600))
    a = CairoMakie.Axis(
        f[1, 1], 
        title = "Histogram of Expected Background Counts, phase $phase data", 
        xlabel = "sum of energies", 
        # ylabel = "expected counts per day", 
        ylabel = "log of expected counts", 
        yscale = log10
        )
    p = CairoMakie.stephist!(a, adjusted_per_day_histos[1], color = colors[1], linewidth = 5)
    for i in 2:length(adjusted_per_day_histos)
        CairoMakie.stephist!(a, adjusted_per_day_histos[i], color = colors[i], linewidth = 5)
    end
    CairoMakie.stephist!(a, sum(adjusted_per_day_histos), color = :black, linewidth = 5)



    labels = [hist.name for hist in all_histos]
    push!(labels, "total")
    elements = [CairoMakie.PolyElement(polycolor = colors[i]) for i in 1:length(labels)-1]
    push!(elements, CairoMakie.PolyElement(polycolor = :black))
    title = "Processes"
    
    ylims!(a, 1e-5, nothing)
    CairoMakie.Legend(f[1,2], elements, labels, title)
    save("stacked_background_histogram_per_day_phase_$phase.png", f, px_per_unit = 2)
    # save("data/out/ex5/stacked_background_histogram_per_day_phase_$phase.png", f, px_per_unit = 2)
    f
end


a,b = 2.7, 3.0

function get_restricted_bkg_counts(processes, a, b; mode = "bkg_only")
    # mode = bkg_only -> only consider background processes
    # mode = all -> consider all processes (default behavior)
    # mode = signal_only -> only consider signal processes
    if mode == "bkg_only"
        processes = [p for p in processes if !p.signal]
    elseif mode == "signal_only"
        processes = [p for p in processes if p.signal]
    end
    return sum(restrict.(get_bkg_counts_1D.(processes), a, b))
end

for i=1:length(processes)
    if processes[i].signal
        continue
    end
    println("i: $i, process: ", processes[i].isotopeName, ", n = ", integral(restrict(get_bkg_counts_1D(processes[i]), a, b)), ",  activity = ", processes[i].activity, " Bq", " n_pass = ", length(processes[i].dataVector))
end

println("Total expected counts: ", integral(get_restricted_bkg_counts(processes, 2.7, 3.0; mode = "bkg_only")))
println("\n\n\n")

for i=1:length(processes)
    processes[i].signal ? continue : nothing
    n_pass = length(processes[i].dataVector)
    eff = n_pass/processes[i].nTotalSim
    a = processes[i].activity
    t = 24*3600
    n_exp_day = round(eff*a*t, digits = 2)
    println("i: $i, process: ", processes[i].isotopeName, ", n_exp_day = ", n_exp_day, ", n_pass = ", n_pass)
end
println("Total expected counts per day: ", round(integral(sum(get_bkg_counts_1D.(processes)))*(1/2.86/365), digits = 2))


# Sensitivity 1D
# Pick signal
signal_matches = get_process(signal_name, processes)
length(signal_matches) == 1 || error("Expected exactly one $signal_name signal in $phase")
signal = only(signal_matches)
background = [p for p in processes if !p.signal]

α = 1.65

t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)
best_sens = get_tHalf(SNparams, effbb, expBkgESum, α; approximate="table")
ThalfbbESum = round(best_sens, sigdigits=3)

# Sensitivity ND
# define variables to use
var_names_nd = [ 
    "energy_elec_sum", 
    "angle_3D_between_ep_em", 
    "delta_y_elec",
    "delta_z_elec",
    "diff_time_elec",
    "energy_elec_1",
    "energy_elec_2",
    ]

# Next we define the boundaries of the variables in the search space. These boundaries will be used to define the search space for the optimization algorithm.
var_bounds_nd = (
    energy_elec_sum = (0.300, 4.00),
    angle_3D_between_ep_em = (0, 180),
    delta_y_elec = (0, 200),
    delta_z_elec = (0, 200),
    diff_time_elec = (0, 10),
    energy_elec_1 = (0, 3.500),
    energy_elec_2 = (0, 3.500),
)


# Default search resolution (adjust if needed)
var_steps_nd = (
    energy_elec_sum = 0.05, 
    angle_3D_between_ep_em = 5.0,
    delta_y_elec = 5.0,
    delta_z_elec = 5.0,
    diff_time_elec = 0.01,
    energy_elec_1 = 0.05,
    energy_elec_2 = 0.05,
)

# Fixed cuts use the given comparison. Optimized cuts use the package's
# one-sided intervals on abs(value): value < max or value >= min.
# Bounds and step define the grid of possible thresholds.
nd_numeric_cuts = [
    (
        name = :vertex_distance,
        value = event -> sqrt(event.delta_y_elec^2 + event.delta_z_elec^2),
        comparison = (<),
        mode = :optimize,
        fixed = 300.0,
        bounds = (50.0, 300.0),
        step = 10.0,
    ),
    (
        name = :closest_gamma,
        value = event -> abs(event.closest_gamma),
        comparison = (>),
        mode = :fixed,
        fixed = 50.0,
        bounds = (50.0, 150.0),
        step = 10.0,
    ),
    (
        name = :closest_elec,
        value = event -> abs(event.closest_elec),
        comparison = (>),
        mode = :fixed,
        fixed = 25.0,
        bounds = (25.0, 125.0),
        step = 10.0,
    ),
    (
        name = :closest_track,
        value = event -> abs(event.closest_track),
        comparison = (>=),
        mode = :fixed,
        fixed = 200.0,
        bounds = (200.0, 400.0),
        step = 10.0,
    ),
    (
        name = :closest_time_track,
        value = event -> abs(event.closest_time_track),
        comparison = (>=),
        mode = :fixed,
        fixed = 10.0,
        bounds = (10.0, 50.0),
        step = 5.0,
    ),
]

for cut in nd_numeric_cuts
    cut.mode in (:fixed, :optimize) || error("Unknown mode for $(cut.name)")
    cut.comparison in ((<), (<=), (>), (>=)) || error("Unknown comparison for $(cut.name)")
end
length(unique(cut.name for cut in nd_numeric_cuts)) == length(nd_numeric_cuts) ||
    error("ND numeric cut names must be unique")
optimized_cuts = [cut for cut in nd_numeric_cuts if cut.mode == :optimize]
cut_column(cut) = Symbol("__nd_cut_", cut.name)
cut_steps(cut) = floor(Int, (cut.bounds[2] - cut.bounds[1]) / cut.step)

nd_indices = findall(
    i -> !selected_files.is_signal[i] || selected_files.name[i] == signal_name,
    1:nrow(selected_files),
)
selected_files_nd = selected_files[nd_indices, :]
selected_tables_nd = data_tables[nd_indices]

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
    filter_nd_bounds(table, var_names_nd, var_bounds_nd)
    for table in selected_tables_nd
]


# Apply fixed topological cuts before calculating numeric cut variables.
for i in eachindex(filtered_data_tables_nd)
    table = filtered_data_tables_nd[i]
    # cut_vertex_dist = sqrt.(table.delta_y_elec .* table.delta_y_elec .+ table.delta_z_elec .* table.delta_z_elec) .>= vertex_prefilter_max
    cut_same_calo   = table.hit_the_same_calo_hit .== 1
    # cut_gam         = abs.(table.closest_gamma) .<= 50
    # cut_elec        = abs.(table.closest_elec) .<= 25
    # cut_track       = abs.(table.closest_track) .< 200
    # cut_ttrack      = abs.(table.closest_time_track) .< 10
    cut_kinks       = table.has_kinks .== 1

    n_evt = length(table)
    cut_om_size = [length(table.num_om_elec_f[j]) != 2 for j in 1:n_evt]

    cut_om_range = falses(n_evt)
    for j in 1:n_evt
        if !cut_om_size[j]
            om0 = table.num_om_elec_f[j][1]
            om1 = table.num_om_elec_f[j][2]
            cut_om_range[j] = !(0 <= om0 < 520 && 0 <= om1 < 520) ||
                              om0 in bad_oms || om1 in bad_oms
        end
    end

    cut_vertex_region = falses(n_evt)
    for j in 1:n_evt
        if !cut_om_size[j] && !cut_om_range[j] && length(table.vertex_3D_start_y[j]) >= 2 && length(table.vertex_3D_start_z[j]) >= 2
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

    keep = .!(cut_same_calo .| cut_kinks .| cut_om_size .| cut_om_range .| cut_vertex_region)
    filtered_data_tables_nd[i] = table[keep]
end

# Optimized cuts are left unfiltered until the ROI is evaluated.
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

# with_variable cannot add a column to an empty process; empty backgrounds
# contribute no events and can be omitted from this optimization.
nd_nonempty = findall(i -> !isempty(filtered_data_tables_nd[i]), eachindex(filtered_data_tables_nd))
for i in eachindex(filtered_data_tables_nd)
    if isempty(filtered_data_tables_nd[i])
        selected_files_nd.is_signal[i] && error("No signal events remain after ND cuts")
        @warn "Skipping empty ND background" name=selected_files_nd.name[i]
    end
end
selected_files_nd = selected_files_nd[nd_nonempty, :]
filtered_data_tables_nd = filtered_data_tables_nd[nd_nonempty]

processes_nd = [
    SNSensitivityEstimate.DataProcessND(
        filtered_data_tables_nd[i],
        String(selected_files_nd.name[i]),
        selected_files_nd.is_signal[i],
        selected_files_nd.activity_bq[i],
        selected_files_nd.time_s[i],
        selected_files_nd.nb_events[i],
        var_bounds_nd,
        selected_files_nd.amount[i],
        var_names_nd,
    ) for i in eachindex(filtered_data_tables_nd)
]

processes_nd = let processes = processes_nd
    for cut in optimized_cuts
        if cut.comparison === (>) || cut.comparison === (<=)
            @warn "Optimized cut uses the package's >= or < interval boundary" name=cut.name
        end
        processes = [
            with_variable(process, cut_column(cut), cut.value, (0.0, Inf))
            for process in processes
        ]
    end
    processes
end

nd_signals = get_process(signal_name, processes_nd)
length(nd_signals) == 1 && only(nd_signals).signal ||
    error("Expected exactly one $signal_name signal in $phase")
count(p -> p.signal, processes_nd) == 1 || error("ND optimization requires exactly one signal")
any(p -> !p.signal, processes_nd) || error("ND optimization requires background processes")
isempty(only(nd_signals).data) && error("No signal events remain after ND cuts")

grid_mins = [Float64(var_bounds_nd[Symbol(name)][1]) for name in var_names_nd]
grid_steps = [Float64(var_steps_nd[Symbol(name)]) for name in var_names_nd]
grid_sizes = [
    round(Int, (var_bounds_nd[Symbol(name)][2] - var_bounds_nd[Symbol(name)][1]) /
               var_steps_nd[Symbol(name)])
    for name in var_names_nd
]
nd_grid = GridSpec(Symbol.(var_names_nd), grid_mins, grid_steps, grid_sizes)

nd_ranges = make_stepRange(nd_grid)
for cut in optimized_cuts
    push!(nd_ranges, (0, cut_steps(cut)))
end
nd_bounds = boxconstraints(
    lb = Float64[first(range) for range in nd_ranges],
    ub = Float64[last(range) for range in nd_ranges],
)

function nd_thresholds(indices)
    offset = 2 * length(var_names_nd)
    return NamedTuple(
        cut.name => min(cut.bounds[2], cut.bounds[1] + round(Int, indices[offset + i]) * cut.step)
        for (i, cut) in enumerate(optimized_cuts)
    )
end

function nd_roi(indices)
    edges = decode_grid_roi(indices[1:2 * length(var_names_nd)], nd_grid)
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

function nd_objective(indices)
    roi = nd_roi(indices)
    edges = Float64[edge for interval in values(roi) for edge in interval]
    return -get_s_to_b(SNparams, α, processes_nd, edges; approximate="formula")
end

nd_algorithm = PSO(; options=Options(time_limit=10.0, verbose=true))
nd_initial = Float64[x for n in nd_grid.nsteps for x in (0, n)]
for cut in optimized_cuts
    initial = clamp(
        round(Int, (cut.fixed - cut.bounds[1]) / cut.step),
        0,
        cut_steps(cut),
    )
    push!(nd_initial, initial)
end
set_user_solutions!(nd_algorithm, nd_initial, nd_objective)

nd_result = optimize(nd_objective, nd_bounds, nd_algorithm)

nd_best_indices = minimizer(nd_result)

nd_best_roi = nd_roi(nd_best_indices)

nd_best_thresholds = nd_thresholds(nd_best_indices)

nd_best_sensitivity = get_sensitivityND(
    SNparams,
    α,
    processes_nd,
    nd_best_roi;
    approximate = "table",
)

nd_best_sensitivity.signalEff > 0 || error("ND optimization found no signal events")

println("\n===== BEST ND RESULT =====")
println("Optimized numeric cuts: ", nd_best_thresholds)
println(nd_best_sensitivity)


test_nd_roit = (
    energy_elec_sum = (2.4, 2.9), 
    angle_3D_between_ep_em = (0.0, 180.0), 
    delta_y_elec = (0.0, 65.0), 
    delta_z_elec = (0.0, 65.0), 
    diff_time_elec = (0.0, 1.75), 
    energy_elec_1 = (0.35, 2.5), 
    energy_elec_2 = (0.35, 2.5), 
    __nd_cut_vertex_distance = (0.0, 200.0)
)
test_nd_sensitivity = get_sensitivityND(
    SNparams,
    α,
    processes_nd,
    test_nd_roit;
    approximate = "table",
)
println("Test ND sensitivity: ", test_nd_sensitivity)