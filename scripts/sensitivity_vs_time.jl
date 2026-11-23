# Sensitivity (90% CL) versus measurement time for the best ROI of ND runs.
# The signal efficiency is constant and the background (best_bkg, counted over
# SNparams["tYear"]) scales linearly with time.
#
#   julia --project=. scripts/sensitivity_vs_time.jl RESULTS.csv [RESULTS.csv ...] \
#       [--max-time=5] [--combined=PATH.png]
#
# Each run directory gets sensitivity_vs_time.{png,pdf,csv}; --combined also
# draws all signals in one figure.

using CSV, DataFramesMeta, CairoMakie
using SNSensitivityEstimate: SNparams, get_tHalf
using CairoMakie.Makie: latexstring

const PROJECT_DIR = dirname(@__DIR__)

# Reference limits: CUPID-0, Phys. Rev. D 107, 032006 (2023).
const SIGNALS = Dict(
    "bb0nu" => (label = L"0\nu\beta\beta", color = "#0072B2", reference = nothing),
    "bb0nuM1" => (label = L"0\nu\beta\beta\chi^0", color = "#0072B2",
                  reference = (1.2e23, L"CUPID-0, $0\nu\beta\beta\chi^0$")),
    "bb0nuM2" => (label = L"0\nu\beta\beta\chi^0\chi^0", color = "#D55E00",
                  reference = (1.4e22, L"CUPID-0, $0\nu\beta\beta\chi^0\chi^0$")),
)
const EXPOSURE_MARKS = (6.25, 17.5) # kg·yr

# Scalar settings (signal, α) are read from settings.jl next to results.csv.
function setting(settings_text, name, path)
    m = match(Regex("^\\s*(?:global\\s+)?$name\\s*=\\s*(\\S+)", "m"), settings_text)
    m === nothing && error("$name is not set in $path")
    return strip(m[1], '"')
end

function read_run(results_path)
    results = CSV.read(results_path, DataFrame)
    nrow(results) == 1 || error("Expected one row in $results_path")
    settings_path = joinpath(dirname(results_path), "settings.jl")
    text = read(settings_path, String)
    row = only(eachrow(results))
    return (dir = dirname(results_path), signal = String(setting(text, "signal_name", settings_path)),
            α = parse(Float64, setting(text, "α", settings_path)), phase = row.phase,
            efficiency = row.best_eff, background = row.best_bkg, sensitivity = row.best_sensitivity)
end

function sensitivity(run, time_yr)
    background = run.background * time_yr / SNparams["tYear"]
    return get_tHalf(SNparams["W"], SNparams["foilMass"], SNparams["Nₐ"], time_yr, SNparams["a"],
                     run.efficiency, background, run.α; approximate = "table")
end

function sensitivity_table(run, times)
    nominal = sensitivity(run, SNparams["tYear"])
    isapprox(nominal, run.sensitivity; rtol = 1e-6) ||
        @warn "Sensitivity at tYear differs from best_sensitivity" run.dir nominal run.sensitivity
    return DataFrame(
        time_yr = times,
        exposure_kg_yr = SNparams["foilMass"] .* times,
        background_counts = run.background .* times ./ SNparams["tYear"],
        sensitivity_yr = [sensitivity(run, t) for t in times],
    )
end

style(signal) = get(() -> (label = latexstring("\\mathrm{", signal, "}"), color = "#333333", reference = nothing),
                    SIGNALS, signal)

function plot_sensitivity(path, runs, tables, max_time)
    mass = SNparams["foilMass"]
    nominal_time = SNparams["tYear"]
    with_theme(theme_latexfonts()) do
        fig = Figure(size = (900, 600), fontsize = 22, figure_padding = (12, 24, 12, 12))
        ax = Axis(fig[1, 1];
            xlabel = L"Exposure (kg$\,\cdot\,$yr)", ylabel = L"$T_{1/2}$ sensitivity at 90% CL (yr)",
            yscale = log10, yminorticksvisible = true, yminorticks = IntervalsBetween(9),
            xminorticksvisible = true, xgridvisible = false, ygridcolor = (:black, 0.08),
            yminorgridvisible = false)
        top = Axis(fig[1, 1]; xaxisposition = :top, xlabel = "Measurement time (yr)",
            xminorticksvisible = true)
        hidespines!(top); hideydecorations!(top); hidexdecorations!(top; label = false, ticklabels = false, ticks = false, minorticks = false)
        xlims!(ax, 0, mass * max_time); xlims!(top, 0, max_time)

        values = vcat([t.sensitivity_yr[t.time_yr .>= 0.1max_time] for t in tables]...)
        references = [style(r.signal).reference for r in runs]
        append!(values, [ref[1] for ref in references if ref !== nothing])
        ylo, yhi = 10.0^floor(log10(minimum(values) / 1.5)), 10.0^ceil(log10(maximum(values) * 1.5))
        yhi / ylo < 100 && (ylo /= 10) # at least two decades, so ticks fall on powers of ten
        ylims!(ax, ylo, yhi)

        for exposure in EXPOSURE_MARKS
            vlines!(ax, exposure; color = (:black, 0.35), linestyle = :dot, linewidth = 1.2)
            text!(ax, exposure, ylo; text = latexstring(exposure, "\\,\\mathrm{kg\\cdot yr}"),
                  rotation = π / 2, align = (:left, :bottom), offset = (-4, 6), fontsize = 15,
                  color = (:black, 0.5))
        end

        elements, labels = Any[], Any[]
        for (run, table) in zip(runs, tables)
            s = style(run.signal)
            line = lines!(ax, table.exposure_kg_yr, table.sensitivity_yr; color = s.color, linewidth = 3)
            point = scatter!(ax, [mass * nominal_time], [sensitivity(run, nominal_time)];
                             color = s.color, markersize = 13, strokecolor = :white, strokewidth = 1.5)
            push!(elements, [line, point]); push!(labels, latexstring("\\mathrm{SuperNEMO},\\ ", s.label.s[2:end-1]))
            if s.reference !== nothing
                value, label = s.reference
                push!(elements, hlines!(ax, value; color = (s.color, 0.8), linestyle = :dash, linewidth = 2))
                push!(labels, label)
            end
        end
        nominal_label = latexstring("\\mathrm{optimized\\ ROI}\\ (", nominal_time, "\\,\\mathrm{yr})")
        push!(elements, MarkerElement(marker = :circle, color = :gray40, markersize = 13))
        push!(labels, nominal_label)
        axislegend(ax, elements, labels; position = :rb, framevisible = true, framecolor = (:black, 0.15),
                   backgroundcolor = (:white, 0.9), labelsize = 18, patchsize = (34, 14), rowgap = 2)

        phase = replace(first(runs).phase, "phase_" => "Phase ")
        Label(fig[0, 1], L"SuperNEMO $^{82}$Se sensitivity"; fontsize = 26, halign = :left, tellwidth = false)
        Label(fig[0, 1], "ND ROI · $phase background model · Feldman–Cousins";
              fontsize = 17, halign = :right, valign = :bottom, color = :gray35, tellwidth = false)
        save(path, fig; px_per_unit = 4)
        save(splitext(path)[1] * ".pdf", fig)
    end
end

function main(args)
    options = Dict{String,String}()
    paths = String[]
    for arg in args
        m = match(r"^--(max-time|combined)=(.*)$", arg)
        m === nothing ? (startswith(arg, "--") ? error("Unknown option $arg") : push!(paths, arg)) :
                        (options[m[1]] = m[2])
    end
    isempty(paths) && error("Usage: julia --project=. scripts/sensitivity_vs_time.jl RESULTS.csv [...] [--max-time=5] [--combined=PATH.png]")
    resolve(p) = abspath(isabspath(p) ? p : joinpath(PROJECT_DIR, p))
    max_time = parse(Float64, get(options, "max-time", "5"))
    max_time > 0 || error("--max-time must be positive")
    times = collect(range(0.0, max_time; length = 201))[2:end]

    runs = [read_run(resolve(p)) for p in paths]
    tables = [sensitivity_table(run, times) for run in runs]
    for (run, table) in zip(runs, tables)
        CSV.write(joinpath(run.dir, "sensitivity_vs_time.csv"), table)
        plot_sensitivity(joinpath(run.dir, "sensitivity_vs_time.png"), [run], [table], max_time)
        @info "Saved sensitivity_vs_time.{png,pdf,csv}" run.dir run.signal sensitivity_at_tYear = sensitivity(run, SNparams["tYear"])
    end
    if haskey(options, "combined")
        path = resolve(options["combined"])
        mkpath(dirname(path))
        plot_sensitivity(path, runs, tables, max_time)
        @info "Saved combined figure" path
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
