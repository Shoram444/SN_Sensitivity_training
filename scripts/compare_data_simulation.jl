# Compare the simulated background spectrum with measured data after the cuts
# of one ND sensitivity run (best ROI from results.csv + fixed cuts from
# settings.jl). Both are normalized to counts per day.
#
#   julia --project=. scripts/compare_data_simulation.jl \
#       --results=results/nd_runs/<run>/results.csv [--name=value ...]
#
# See scripts/README.md ("Data vs simulation") for all options.

using CSV, DataFramesMeta, UnROOT, Distributions
using SNSensitivityEstimate: SNparams

const PROJECT_DIR = dirname(@__DIR__)
include(joinpath(PROJECT_DIR, "analysis", "nd_background_report.jl"))
include(joinpath(PROJECT_DIR, "analysis", "dead_om_map.jl"))

const SECONDS_PER_DAY = 86400.0
const SIM_PREFIX = "data/mathis_files/"

const OPTIONS = Dict(
    "results" => "results.csv of the ND run (required)",
    "settings" => "settings.jl with the fixed cuts [settings.jl next to results]",
    "metadata" => "simulation metadata CSV [data/mathis_data_info_final_2.csv]",
    "data" => "measured data ROOT file [data/real_data/run_3470-4427.root]",
    "run-times" => "CSV from scripts/list_run_times.sh [data/real_data/run_times.csv]",
    "runs" => "'data' (runs present in the data file) or FIRST:LAST [data]",
    "sim-dir" => "directory replacing data/mathis_files in metadata paths, e.g. data/mathis_files_local",
    "output-dir" => "output directory [directory of results]",
    "edges" => "energy-sum binning MIN:STEP:MAX in MeV [0:0.05:4]",
    "chi2-range" => "energy-sum range MIN:MAX in MeV used for chi2 [whole binning]",
)

const BASE_VAR_NAMES = [
    "energy_elec_sum", "angle_3D_between_ep_em", "delta_y_elec", "delta_z_elec",
    "diff_time_elec", "energy_elec_1", "energy_elec_2", "hit_the_same_calo_hit",
    "closest_gamma", "closest_elec", "closest_track", "closest_time_track", "has_kinks",
    "num_om_elec_f", "vertex_3D_start_y", "vertex_3D_start_z",
]

const OM_LIST_FILES = Dict(
    "phase_0" => "list_om_phase_0.txt",
    "phase_1" => "list_om_phase_1.txt",
    "phase_2" => "list_om_phase_2a.txt",
    "phase_3" => "list_om_phase_2b.txt",
    "phase_4" => "list_om_phase_3.txt",
)

function usage()
    lines = ["Usage: julia --project=. scripts/compare_data_simulation.jl --results=PATH [--name=value ...]", ""]
    append!(lines, ["  --$name  $help" for (name, help) in sort!(collect(OPTIONS))])
    return join(lines, '\n')
end

function parse_options(args)
    options = Dict{String,String}()
    for arg in args
        m = match(r"^--([a-z0-9-]+)=(.*)$", arg)
        m === nothing && error("Expected --name=value, got '$arg'\n\n$(usage())")
        haskey(OPTIONS, m[1]) || error("Unknown option --$(m[1])\n\n$(usage())")
        options[m[1]] = m[2]
    end
    haskey(options, "results") || error("Missing --results\n\n$(usage())")
    project_path(p) = isabspath(p) ? p : joinpath(PROJECT_DIR, p)
    results = abspath(project_path(options["results"]))
    get_path(name, default) = abspath(project_path(get(options, name, default)))
    edges = parse_edges(get(options, "edges", "0:0.05:4"))
    chi2_range = haskey(options, "chi2-range") ?
        parse_interval(options["chi2-range"]) : (first(edges), last(edges))
    return (
        results = results,
        settings = haskey(options, "settings") ? get_path("settings", "") :
                   joinpath(dirname(results), "settings.jl"),
        metadata = get_path("metadata", "data/mathis_data_info_final_2.csv"),
        data = get_path("data", "data/real_data/run_3470-4427.root"),
        run_times = get_path("run-times", "data/real_data/run_times.csv"),
        runs = parse_runs(get(options, "runs", "data")),
        sim_dir = haskey(options, "sim-dir") ? get_path("sim-dir", "") : nothing,
        output_dir = haskey(options, "output-dir") ? get_path("output-dir", "") : dirname(results),
        edges = edges,
        chi2_range = chi2_range,
    )
end

function parse_edges(text)
    parts = parse.(Float64, split(text, ':'))
    length(parts) == 3 || error("Edges must be MIN:STEP:MAX, got '$text'")
    lo, step, hi = parts
    n = (hi - lo) / step
    step > 0 && hi > lo && isapprox(n, round(n); atol = 1e-8) ||
        error("Edges '$text' must contain a positive whole number of bins")
    return range(lo, hi; length = round(Int, n) + 1)
end

function parse_interval(text)
    parts = parse.(Float64, split(text, ':'))
    length(parts) == 2 && parts[1] < parts[2] || error("Expected MIN:MAX, got '$text'")
    return (parts[1], parts[2])
end

function parse_runs(text)
    text == "data" && return :data
    parts = parse.(Int, split(text, ':'))
    length(parts) == 2 && parts[1] <= parts[2] || error("Expected --runs=data or FIRST:LAST, got '$text'")
    return parts[1]:parts[2]
end

# results.csv stores the ROI as a printed NamedTuple; read it without eval.
function parse_roi(text)
    expr = Meta.parse(text)
    entries = expr isa Expr && expr.head === :(=) ? [expr] :
              expr isa Expr && expr.head === :tuple ? expr.args :
              error("Cannot parse ROI: $text")
    return NamedTuple(map(entries) do entry
        entry isa Expr && entry.head === :(=) && entry.args[1] isa Symbol ||
            error("Cannot parse ROI entry: $entry")
        interval = entry.args[2]
        interval isa Expr && interval.head === :tuple && length(interval.args) == 2 ||
            error("Cannot parse ROI interval: $interval")
        lo, hi = roi_number.(interval.args)
        lo < hi || error("Empty ROI interval for $(entry.args[1])")
        entry.args[1] => (lo, hi)
    end)
end

roi_number(x::Real) = Float64(x)
roi_number(x::Symbol) = x === :Inf ? Inf : error("Unexpected ROI value: $x")
function roi_number(x::Expr)
    x.head === :call && x.args[1] === :- && length(x.args) == 2 ||
        error("Unexpected ROI value: $x")
    return -roi_number(x.args[2])
end

function read_results(path)
    df = CSV.read(path, DataFrame; types = Dict(:best_ROI => String, :phase => String))
    nrow(df) == 1 || error("Expected one row in $path, found $(nrow(df))")
    row = only(eachrow(df))
    return (phase = row.phase, roi = parse_roi(row.best_ROI), sensitivity = row.best_sensitivity,
            bkg = row.best_bkg, eff = row.best_eff)
end

function load_settings(path)
    isfile(path) || error("Settings not found: $path")
    m = Module(:RunSettings)
    # settings.jl reads the output path from the batch environment.
    withenv("SN_OUTPUT_PATH" => get(ENV, "SN_OUTPUT_PATH", "")) do
        Base.include(m, path)
    end
    get_setting(name) = isdefined(m, name) ? getfield(m, name) : error("$name is not set in $path")
    return (
        phase = get_setting(:phase),
        signal_name = get_setting(:signal_name),
        background_names = get_setting(:background_names),
        extra_var_names = String.(get_setting(:extra_var_names)),
        use_same_calo_cut = get_setting(:use_same_calo_cut),
        use_kinks_cut = get_setting(:use_kinks_cut),
        use_om_cut = get_setting(:use_om_cut),
        use_vertex_region_cut = get_setting(:use_vertex_region_cut),
        use_bad_om_list = isdefined(m, :use_bad_om_list) ? getfield(m, :use_bad_om_list) : true,
        use_dead_om_map = isdefined(m, :use_dead_om_map) ? getfield(m, :use_dead_om_map) : false,
        nd_numeric_cuts = get_setting(:nd_numeric_cuts),
    )
end

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

const CUT_PREFIX = "__nd_cut_"

# Same cuts as the ND analysis: ROI windows on abs(value) (energy sum excluded,
# it is histogrammed), topological cuts, fixed numeric cuts and the optimized
# numeric thresholds from the best ROI.
function build_selection(settings, roi, bad_oms; dead_om = nothing)
    haskey(roi, :energy_elec_sum) || error("The best ROI has no energy_elec_sum window")
    cuts = Dict(cut.name => cut for cut in settings.nd_numeric_cuts)
    windows = Pair{Symbol,Tuple{Float64,Float64}}[]
    optimized = Tuple{Any,Tuple{Float64,Float64}}[]
    for (key, interval) in pairs(roi)
        name = String(key)
        if startswith(name, CUT_PREFIX)
            cut_name = Symbol(name[length(CUT_PREFIX)+1:end])
            haskey(cuts, cut_name) || error("ROI cut $cut_name is not defined in settings.jl")
            cuts[cut_name].mode === :optimize || error("ROI cut $cut_name is not optimized in settings.jl")
            push!(optimized, (cuts[cut_name], interval))
        elseif key !== :energy_elec_sum
            push!(windows, key => interval)
        end
    end
    for cut in settings.nd_numeric_cuts
        cut.mode === :optimize && !any(c === cut for (c, _) in optimized) &&
            error("Optimized cut $(cut.name) is missing from the best ROI")
    end
    var_names = unique!([BASE_VAR_NAMES; settings.extra_var_names; String.(first.(windows))])
    return (
        windows = windows,
        fixed = [cut for cut in settings.nd_numeric_cuts if cut.mode === :fixed],
        optimized = optimized,
        use_same_calo_cut = settings.use_same_calo_cut,
        use_kinks_cut = settings.use_kinks_cut,
        use_om_cut = settings.use_om_cut,
        use_vertex_region_cut = settings.use_vertex_region_cut,
        bad_oms = bad_oms,
        dead_om = dead_om,
        var_names = var_names,
    )
end

function topological_veto(table, sel)
    n_evt = length(table)
    cut_same_calo = table.hit_the_same_calo_hit .== 1
    cut_kinks = table.has_kinks .== 1
    cut_om_size = [length(table.num_om_elec_f[j]) != 2 for j in 1:n_evt]
    cut_om_range = falses(n_evt)
    cut_vertex_region = falses(n_evt)
    for j in 1:n_evt
        if !cut_om_size[j]
            om0, om1 = table.num_om_elec_f[j]
            cut_om_range[j] = !(0 <= om0 < 520 && 0 <= om1 < 520) ||
                              om0 in sel.bad_oms || om1 in sel.bad_oms
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
    return (sel.use_same_calo_cut .& cut_same_calo) .| (sel.use_kinks_cut .& cut_kinks) .|
           (sel.use_om_cut .& (cut_om_size .| cut_om_range)) .|
           (sel.use_vertex_region_cut .& cut_vertex_region)
end

# abs(energy_elec_sum) of the events passing every cut except the energy window.
# data = true takes dead OM veto times from current_timestamp.
function select_energies(table, sel; data = false)
    keep = trues(length(table))
    for (name, (lo, hi)) in sel.windows
        values = abs.(getproperty(table, name))
        keep .&= (lo .<= values) .& (values .< hi)
    end
    any(keep) || return Float64[]
    table = table[keep]
    keep = .!topological_veto(table, sel)
    if sel.dead_om !== nothing
        intervals, live = sel.dead_om
        keep .&= .!dead_om_veto(table, intervals, live; times = data ? collect(table.current_timestamp) : nothing)
    end
    for cut in sel.fixed
        keep .&= cut.comparison.(map(cut.value, table), cut.fixed)
    end
    for (cut, (lo, hi)) in sel.optimized
        values = abs.(map(cut.value, table))
        keep .&= (lo .<= values) .& (values .< hi)
    end
    return Float64.(abs.(collect(table.energy_elec_sum)[keep]))
end

function bin_counts(energies, edges)
    counts = zeros(Int, length(edges) - 1)
    for e in energies
        first(edges) <= e < last(edges) && (counts[searchsortedlast(edges, e)] += 1)
    end
    return counts
end

in_window(energies, (lo, hi)) = count(e -> lo <= e < hi, energies)

function simulation_path(file_path, sim_dir)
    if sim_dir !== nothing && startswith(file_path, SIM_PREFIX)
        return joinpath(sim_dir, file_path[length(SIM_PREFIX)+1:end])
    end
    return isabspath(file_path) ? file_path : joinpath(PROJECT_DIR, file_path)
end

function background_files(metadata_path, settings)
    info = CSV.read(metadata_path, DataFrame; types = Dict(:file_path => String, :name => String))
    rows = info[(info.phase .== settings.phase) .& .!info.is_signal .&
                .!occursin.("reduces", info.file_path), :]
    if settings.background_names !== :all
        requested = Set(String.(settings.background_names))
        missing_names = setdiff(requested, Set(rows.name))
        isempty(missing_names) || error("Unknown backgrounds for $(settings.phase): $missing_names")
        rows = rows[in.(rows.name, Ref(requested)), :]
    end
    nrow(rows) > 0 || error("No backgrounds for $(settings.phase) in $metadata_path")
    allunique(rows.name) || error("Duplicate background names for $(settings.phase) in $metadata_path")
    return rows
end

# Live time of the selected runs. With runs === :data, the runs present in the data.
function live_time(run_times_path, data_runs, runs)
    isfile(run_times_path) || error("Run time CSV not found: $run_times_path (create it with scripts/list_run_times.sh)")
    df = CSV.read(run_times_path, DataFrame; types = Dict(:file_path => String))
    unnamed = df.file_path[ismissing.(df.run_number)]
    isempty(unnamed) || error("Missing run number in $run_times_path for: $(join(unnamed, ", "))")
    duplicated = [r for r in unique(df.run_number) if count(==(r), df.run_number) > 1]
    isempty(duplicated) || error(
        "Runs listed more than once in $run_times_path (merged or duplicate files?): ",
        join(["$r ($(join(df.file_path[df.run_number .== r], ", ")))" for r in duplicated], "; "))
    # Rows without run_time are files that could not be read; fatal only if selected.
    failed = Set(Int.(df.run_number[ismissing.(df.run_time)]))
    times = Dict(Int(r) => Float64(t) for (r, t) in zip(df.run_number, df.run_time) if !ismissing(t))
    data_set = Set(data_runs)
    selected = runs === :data ? sort!(collect(data_set)) : [r for r in runs if haskey(times, r)]
    without_time = [r for r in selected if !haskey(times, r)]
    isempty(without_time) || error("No run time for data runs: $(join(without_time, ", "))")
    unused_failed = sort!(collect(setdiff(failed, data_set)))
    runs === :data || filter!(!in(runs), unused_failed)
    isempty(unused_failed) || @warn "Ignoring unreadable run files outside the selection" runs = unused_failed
    failed_in_range = runs === :data ? Int[] : sort!(collect(intersect(failed, Set(runs))))
    isempty(failed_in_range) || error("No run time for runs in the selected range: $(join(failed_in_range, ", "))")
    if runs !== :data
        no_events = [r for r in selected if !(r in data_set)]
        isempty(no_events) || @warn "Runs with a run time but no data events are counted in the live time" runs = no_events seconds = sum(times[r] for r in no_events)
    end
    isempty(selected) && error("No runs selected")
    return sum(times[r] for r in selected), selected
end

# Poisson likelihood-ratio chi2 (Baker-Cousins) of one bin, with the MC
# normalization of the bin profiled within its statistical uncertainty
# (Barlow-Beeston lite, Conway 2011). Valid for small counts.
function bin_chi2(n, mu, variance)
    tau2 = variance / mu^2
    beta = tau2 > 0 ? ((1 - mu * tau2) + sqrt((mu * tau2 - 1)^2 + 4n * tau2)) / 2 : 1.0
    m = beta * mu
    deviance = 2 * (m - n + (n > 0 ? n * log(n / m) : 0.0))
    return deviance + (tau2 > 0 ? (beta - 1)^2 / tau2 : 0.0)
end

# chi2 sums bin_chi2 over bins with mu > 0 (nothing is fitted, so ndf is the
# number of bins). Pearson's chi2 is kept for reference; it is unreliable when
# bins expect much less than one event.
function chi2_test(data_counts, model_counts, model_variance, mask)
    chi2 = 0.0
    pearson = 0.0
    ndf = 0
    unpredicted = 0
    for i in eachindex(data_counts)
        mask[i] || continue
        n, mu = data_counts[i], model_counts[i]
        if mu <= 0
            unpredicted += n > 0
            continue
        end
        ndf += 1
        chi2 += bin_chi2(n, mu, model_variance[i])
        pearson += (n - mu)^2 / (mu + model_variance[i])
    end
    return (chi2 = chi2, ndf = ndf, p_value = ndf > 0 ? ccdf(Chisq(ndf), chi2) : NaN,
            chi2_pearson = pearson, unpredicted_bins = unpredicted)
end

# Signed sqrt of each bin's chi2 contribution; their squares sum to chi2.
function residuals(data_counts, model_counts, model_variance)
    return [mu > 0 ? sign(n - mu) * sqrt(max(bin_chi2(n, mu, v), 0.0)) : NaN
            for (n, mu, v) in zip(data_counts, model_counts, model_variance)]
end

# Central 68.27% Poisson (Garwood) interval for n observed events.
function poisson_interval(n)
    alpha = (1 - 0.6826894921370859) / 2
    lo = n == 0 ? 0.0 : quantile(Chisq(2n), alpha) / 2
    return lo, quantile(Chisq(2n + 2), 1 - alpha) / 2
end

function step_band(edges, lower, upper)
    x = Float64[]
    lo = Float64[]
    hi = Float64[]
    for i in eachindex(lower)
        append!(x, (edges[i], edges[i+1]))
        append!(lo, (lower[i], lower[i]))
        append!(hi, (upper[i], upper[i]))
    end
    return x, lo, hi
end

function format_power(x)
    exponent = floor(Int, log10(x))
    superscript = map(c -> c == '-' ? '⁻' : "⁰¹²³⁴⁵⁶⁷⁸⁹"[nextind("⁰¹²³⁴⁵⁶⁷⁸⁹", 0, c - '0' + 1)], string(exponent))
    return "$(round(x / 10.0^exponent; digits = 2))×10$superscript"
end

const SIMULATION_COLOR = "#E3000F"
# ²¹⁴Bi (foil) is red-orange in the background report; keep red for the total.
const COMPARISON_STYLES = merge(BACKGROUND_STYLES, Dict("bi_foil" => ("²¹⁴Bi (foil)", "#E69F00", :solid)))

# data_low/data_high: lower and upper data errors (distances from data_rate).
function save_comparison_plot(path, edges, spectra, variances, data_rate, data_low, data_high,
                              pulls, energy_roi, chi2_range, stats_labels, title, subtitle)
    ink = "#172333"
    roi_color = "#6686AD"
    fig = CairoMakie.Figure(size = (1180, 820), fontsize = 19, figure_padding = 24)
    ax = CairoMakie.Axis(
        fig[1, 1];
        title, subtitle,
        titlealign = :left, titlesize = 27, subtitlesize = 17,
        titlecolor = ink, subtitlecolor = "#475569",
        ylabel = "Events / day / $(round(Int, 1000 * step(edges))) keV",
        yscale = log10, xticks = 0:0.5:last(edges), xticklabelsvisible = false,
        xgridvisible = false, ygridcolor = "#E1E6EC",
        topspinevisible = false, rightspinevisible = false,
    )
    pull_ax = CairoMakie.Axis(
        fig[2, 1];
        xlabel = "Electron energy sum, E₁ + E₂ (MeV)", ylabel = "Residual (σ)",
        xticks = 0:0.5:last(edges), xgridvisible = false, ygridcolor = "#E1E6EC",
        topspinevisible = false, rightspinevisible = false,
    )
    CairoMakie.linkxaxes!(ax, pull_ax)
    CairoMakie.rowsize!(fig.layout, 2, CairoMakie.Relative(0.27))
    CairoMakie.rowgap!(fig.layout, 1, 8)

    total = sum(values(spectra))
    total_error = sqrt.(sum(values(variances)))
    positive = filter(>(0), [total; data_rate])
    ymin = isempty(positive) ? 1e-5 : minimum(positive) / 5
    ymax = isempty(positive) ? 1.0 : 3 * maximum([total .+ total_error; data_rate .+ data_high])
    baseline = ymin / 2
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2

    roi_band = CairoMakie.vspan!(ax, energy_roi...; color = (roi_color, 0.17))
    CairoMakie.vspan!(pull_ax, energy_roi...; color = (roi_color, 0.17))
    for axis in (ax, pull_ax)
        CairoMakie.vlines!(axis, collect(energy_roi); color = (roi_color, 0.75),
                           linestyle = :dash, linewidth = 1.2)
    end

    names = sort!(collect(keys(spectra)); by = name -> (-sum(spectra[name]), name))
    fallback_colors = CairoMakie.resample_cmap(:Dark2_8, max(2, length(names)))
    components = map(enumerate(names)) do (i, name)
        label, color, style = get(COMPARISON_STYLES, name, (name, fallback_colors[i], :dot))
        x, y = spectrum_steps(edges, spectra[name], baseline)
        line = CairoMakie.stairs!(ax, x, y; step = :post, color, linestyle = style, linewidth = 1.8)
        (line, label)
    end
    band_x, band_lo, band_hi = step_band(edges, max.(baseline, total .- total_error), max.(baseline, total .+ total_error))
    total_band = CairoMakie.band!(ax, band_x, band_lo, band_hi; color = (SIMULATION_COLOR, 0.22))
    x, y = spectrum_steps(edges, total, baseline)
    total_line = CairoMakie.stairs!(ax, x, y; step = :post, color = SIMULATION_COLOR, linewidth = 3)

    shown = data_rate .> 0
    lower = min.(data_low[shown], data_rate[shown] .- baseline)
    CairoMakie.errorbars!(ax, centers[shown], data_rate[shown], lower, data_high[shown];
                          color = :black, whiskerwidth = 0, linewidth = 1.4)
    data_points = CairoMakie.scatter!(ax, centers[shown], data_rate[shown]; color = :black, markersize = 8)
    CairoMakie.xlims!(ax, first(edges), last(edges))
    CairoMakie.ylims!(ax, ymin, ymax)
    # Whole decades only; the automatic log ticks may fall on half decades.
    decades = ceil(Int, log10(ymin)):floor(Int, log10(ymax))
    ax.yticks = (10.0 .^ decades, [CairoMakie.rich("10", CairoMakie.superscript(string(d))) for d in decades])

    finite = isfinite.(pulls)
    limit = clamp(1.15 * maximum(abs, pulls[finite]; init = 0.0), 3.0, 10.0)
    CairoMakie.hspan!(pull_ax, -2, 2; color = ("#F2C14E", 0.25))
    CairoMakie.hspan!(pull_ax, -1, 1; color = ("#3BA55C", 0.25))
    CairoMakie.hlines!(pull_ax, [0.0]; color = ink, linewidth = 1)
    # Bins outside the chi2 range are faded; clipped residuals are drawn as triangles.
    used = chi2_range[1] .<= centers .< chi2_range[2]
    for (in_chi2, color) in ((true, :black), (false, (:black, 0.3)))
        shown = finite .& (used .== in_chi2)
        marker = [abs(p) > limit ? (p > 0 ? :utriangle : :dtriangle) : :circle for p in pulls[shown]]
        CairoMakie.scatter!(pull_ax, centers[shown], clamp.(pulls[shown], -0.94limit, 0.94limit);
                            color, marker, markersize = 8)
    end
    CairoMakie.ylims!(pull_ax, -limit, limit)

    CairoMakie.Legend(
        fig[1:2, 2],
        [data_points, [total_band, total_line], roi_band, first.(components)...,
         (CairoMakie.MarkerElement(marker = :circle, color = :transparent) for _ in stats_labels)...],
        ["Data", "Total simulation (±1σ MC stat.)", "Energy ROI", last.(components)..., stats_labels...];
        framevisible = false, valign = :center, patchsize = (30, 14), labelsize = 16,
    )
    CairoMakie.save(path, fig; px_per_unit = 3)
    CairoMakie.save(splitext(path)[1] * ".pdf", fig)
end

function run_comparison(config, settings)
    results = read_results(config.results)
    results.phase == settings.phase ||
        error("Phase in $(config.results) ($(results.phase)) differs from settings.jl ($(settings.phase))")
    # The bad OM list follows the run's setting, for simulations and data alike.
    bad_oms = if settings.use_bad_om_list
        haskey(OM_LIST_FILES, settings.phase) || error("No OM list configured for $(settings.phase)")
        read_bad_oms(joinpath(PROJECT_DIR, "data", "om_lists", OM_LIST_FILES[settings.phase]))
    else
        Set{Int}()
    end
    @info "Bad OM list" settings.use_bad_om_list length(bad_oms)
    # Same dead OM map treatment as the analysis: data by event time, simulation
    # by a time drawn from the live time of the whole data file.
    dead_om = if settings.use_dead_om_map
        (read_dead_om_map(joinpath(PROJECT_DIR, DEAD_OM_MAP_FILE), joinpath(PROJECT_DIR, OM_GID_MAP_FILE)),
         data_live_time(config.data, config.run_times))
    else
        nothing
    end
    @info "Dead OM map" settings.use_dead_om_map
    sel = build_selection(settings, results.roi, bad_oms; dead_om)
    edges = config.edges
    energy_roi = results.roi.energy_elec_sum
    @info "Best ROI" results.roi
    @info "Fixed numeric cuts" [(cut.name, cut.comparison, cut.fixed) for cut in sel.fixed]

    # Simulation: counts/day = activity * amount * 86400 / generated events per passing event.
    files = background_files(config.metadata, settings)
    spectra = Dict{String,Vector{Float64}}()
    variances = Dict{String,Vector{Float64}}()
    summary = DataFrame(process = String[], events = Int[], events_roi = Int[],
                        rate_per_day = Float64[], rate_error = Float64[],
                        rate_per_day_roi = Float64[], rate_error_roi = Float64[])
    for row in eachrow(files)
        path = simulation_path(row.file_path, config.sim_dir)
        isfile(path) || error("Simulation file not found: $path (use --sim-dir for a local copy)")
        @info "Reading $(row.name)" path
        table = LazyTree(ROOTFile(path), "Event", sel.var_names)
        energies = select_energies(table, sel; data = false)
        counts = bin_counts(energies, edges)
        n_roi = in_window(energies, energy_roi)
        weight = row.amount * row.activity_bq * SECONDS_PER_DAY / row.nb_events
        push!(summary, (row.name, sum(counts), n_roi, weight * sum(counts),
                        weight * sqrt(sum(counts)), weight * n_roi, weight * sqrt(n_roi)))
        group = occursin("pmt", row.name) ? "pmt" : row.name
        get!(() -> zeros(length(counts)), spectra, group) .+= weight .* counts
        get!(() -> zeros(length(counts)), variances, group) .+= weight^2 .* counts
    end
    # The analysis normalized the ROI background to the full exposure SNparams["t"].
    expected_bkg = sum(summary.rate_per_day_roi) * SNparams["t"] / SECONDS_PER_DAY
    isapprox(expected_bkg, results.bkg; rtol = 1e-6) ||
        @warn "Simulated ROI background differs from best_bkg in results.csv (changed metadata or cuts?)" expected_bkg results.bkg
    mc_rate = sum(values(spectra))
    mc_variance = sum(values(variances))
    push!(summary, ("Total simulation", sum(summary.events), sum(summary.events_roi),
                    sum(summary.rate_per_day), sqrt(sum(abs2, summary.rate_error)),
                    sum(summary.rate_per_day_roi), sqrt(sum(abs2, summary.rate_error_roi))))

    # Data, normalized by the summed run_time of the selected runs.
    @info "Reading data" config.data
    data_table = LazyTree(ROOTFile(config.data), "Event", [sel.var_names; "run_number"; "unix_start_time"; "current_timestamp"])
    event_runs = Int.(collect(data_table.run_number))
    data_runs = unique(event_runs)
    if config.runs !== :data
        in_range = in(Set(config.runs)).(event_runs)
        any(in_range) || error("No data events in runs $(config.runs)")
        data_table = data_table[in_range]
        data_runs = unique(event_runs[in_range])
    end
    live_seconds, runs = live_time(config.run_times, data_runs, config.runs)
    live_days = live_seconds / SECONDS_PER_DAY
    # Some runs store unix_start_time = -1; the span is only a sanity check.
    start_times = filter(>(0), collect(data_table.unix_start_time))
    span_days = isempty(start_times) ? NaN : (maximum(start_times) - minimum(start_times)) / SECONDS_PER_DAY
    @info "Data live time" runs = length(runs) first_run = first(runs) last_run = last(runs) live_days span_of_run_start_days = span_days
    data_energies = select_energies(data_table, sel; data = true)
    data_counts = bin_counts(data_energies, edges)
    data_roi = in_window(data_energies, energy_roi)
    push!(summary, ("Data", sum(data_counts), data_roi, sum(data_counts) / live_days,
                    sqrt(sum(data_counts)) / live_days, data_roi / live_days, sqrt(data_roi) / live_days))

    model_counts = mc_rate .* live_days
    model_variance = mc_variance .* live_days^2
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    chi2_mask = config.chi2_range[1] .<= centers .< config.chi2_range[2]
    test = chi2_test(data_counts, model_counts, model_variance, chi2_mask)
    test.unpredicted_bins > 0 &&
        @warn "Bins with data but no simulated events are excluded from chi2" bins = test.unpredicted_bins
    pulls = residuals(data_counts, model_counts, model_variance)

    data_rate = data_counts ./ live_days
    intervals = poisson_interval.(data_counts)
    data_low = (data_counts .- first.(intervals)) ./ live_days
    data_high = (last.(intervals) .- data_counts) ./ live_days
    bins = DataFrame(
        energy_lo = edges[1:end-1], energy_hi = edges[2:end],
        data_events = data_counts, data_rate_per_day = data_rate,
        data_rate_error_low = data_low, data_rate_error_high = data_high,
        sim_rate_per_day = mc_rate, sim_rate_error = sqrt.(mc_variance),
        residual = pulls, in_chi2 = chi2_mask,
    )
    for name in sort!(collect(keys(spectra)))
        bins[!, "sim_$(name)_rate_per_day"] = spectra[name]
    end
    stats = DataFrame(
        quantity = ["phase", "signal", "best_sensitivity_yr", "best_bkg", "best_eff",
                    "runs", "first_run", "last_run", "live_time_s", "live_time_days",
                    "chi2_range_mev", "chi2", "ndf", "chi2_per_ndf", "p_value",
                    "chi2_pearson", "chi2_pearson_per_ndf", "excluded_bins_without_simulation"],
        value = string.([settings.phase, settings.signal_name, results.sensitivity, results.bkg,
                         results.eff, length(runs), first(runs), last(runs), live_seconds, live_days,
                         "$(config.chi2_range[1])-$(config.chi2_range[2])", test.chi2, test.ndf,
                         test.chi2 / test.ndf, test.p_value, test.chi2_pearson,
                         test.chi2_pearson / test.ndf, test.unpredicted_bins]),
    )

    mkpath(config.output_dir)
    CSV.write(joinpath(config.output_dir, "data_vs_simulation_bins.csv"), bins)
    CSV.write(joinpath(config.output_dir, "data_vs_simulation_summary.csv"), summary)
    CSV.write(joinpath(config.output_dir, "data_vs_simulation_stats.csv"), stats)

    phase_label = replace(settings.phase, "phase_" => "Phase ")
    roi_lo, roi_hi = round.(energy_roi; digits = 3)
    stats_labels = ["χ²/ndf = $(round(test.chi2; digits = 1)) / $(test.ndf) = $(round(test.chi2 / test.ndf; digits = 2))"]
    subtitle = "$phase_label · Runs $(first(runs))–$(last(runs)) ($(length(runs)) runs, " *
               "$(round(live_days; digits = 2)) d) · Energy ROI: $(roi_lo)–$roi_hi MeV · " *
               "T½ > $(format_power(results.sensitivity)) yr ($(settings.signal_name))"
    plot_path = joinpath(config.output_dir, "data_vs_simulation.png")
    save_comparison_plot(plot_path, edges, spectra, variances, data_rate, data_low, data_high,
                         pulls, energy_roi, config.chi2_range, stats_labels,
                         "SuperNEMO | Data vs simulation", subtitle)

    println("\n===== DATA VS SIMULATION (COUNTS / DAY) =====")
    show(stdout, MIME"text/plain"(), summary; allrows = true, allcols = true)
    println("\n")
    show(stdout, MIME"text/plain"(), stats; allrows = true, allcols = true)
    println("\nSaved $plot_path and CSV files to $(config.output_dir)")
    return (; summary, bins, stats, test)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if any(in(("-h", "--help")), ARGS)
        println(usage())
    else
        config = parse_options(ARGS)
        settings = load_settings(config.settings)
        # settings.jl defines closures; run in the latest world so they can be called.
        Base.invokelatest(run_comparison, config, settings)
    end
end
