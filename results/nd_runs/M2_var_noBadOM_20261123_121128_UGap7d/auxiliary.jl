# Included after optimization; all event columns are already in memory.
import CairoMakie

function background_counts(table, roi, energy_prefilter, edges)
    energy_lo, energy_hi = roi[:energy_elec_sum]
    prefilter_lo, prefilter_hi = energy_prefilter
    windows = [(key, interval) for (key, interval) in pairs(roi) if key !== :energy_elec_sum]
    n_roi = 0
    bins = zeros(Int, length(edges) - 1)

    for event in table
        all(lo <= abs(getproperty(event, key)) < hi for (key, (lo, hi)) in windows) || continue
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

const BACKGROUND_STYLES = Dict(
    "bb"            => ("2νββ",               "#0072B2", :solid),
    "bi_foil"       => ("²¹⁴Bi (foil)",       "#D55E00", :solid),
    "k"             => ("⁴⁰K (foil)",         "#009E73", :solid),
    "pa"            => ("²³⁴ᵐPa (foil)",      "#CC79A7", :solid),
    "tl"            => ("²⁰⁸Tl (foil)",       "#AA7700", :dash),
    "radon"         => ("Radon (field wires)", "#377EB8", :dash),
    "radon_surface" => ("Radon (surface)",    "#984EA3", :dot),
    "pmt"           => ("PMT glass",          "#6C503A", :dashdot),
)

function spectrum_steps(edges, bins, baseline)
    # Close at both outer edges; zero bins fall below the visible log range.
    return [first(edges); edges], [baseline; map(n -> n > 0 ? n : baseline, bins); baseline]
end

const SIGNAL_LATEX = Dict(
    "bb0nu"   => raw"0\nu\beta\beta",
    "bb0nuM1" => raw"0\nu\beta\beta\chi^{0}\;(n=1)",
    "bb0nuM2" => raw"0\nu\beta\beta\chi^{0}\;(n=2)",
)

function spectrum_subtitle(phase, energy_roi, signal_name, exposure)
    phase_label = replace(phase, "phase_" => "Phase ")
    roi_lo, roi_hi = round.(energy_roi; digits = 3)
    exposure_label = "$(round(exposure; digits = 2)) kg·yr"
    details = "$phase_label · $exposure_label · Energy ROI: $(roi_lo)–$roi_hi MeV"
    signal_tex = get(SIGNAL_LATEX, signal_name, nothing)
    signal_tex === nothing && return "$details · Signal: $signal_name"
    return CairoMakie.Makie.latexstring(raw"\text{", details, raw" · Signal: }", signal_tex)
end

function save_background_spectrum(path, edges, spectra, variances, energy_roi, phase,
                                  signal_name, exposure)
    ink = "#172333"
    roi_color = "#6686AD"
    fig = CairoMakie.Figure(size = (1130, 610), fontsize = 19, figure_padding = 24)
    ax = CairoMakie.Axis(
        fig[1, 1];
        title = "SuperNEMO | Expected background spectrum",
        subtitle = spectrum_subtitle(phase, energy_roi, signal_name, exposure),
        titlealign = :left,
        titlesize = 27,
        subtitlesize = 18,
        titlecolor = ink,
        subtitlecolor = "#475569",
        xlabel = "Electron energy sum, E₁ + E₂ (MeV)",
        ylabel = "Expected events per 50 keV",
        yscale = log10,
        xticks = 0:0.5:4,
        xgridvisible = false,
        ygridcolor = "#E1E6EC",
        topspinevisible = false,
        rightspinevisible = false,
    )
    roi_band = CairoMakie.vspan!(ax, energy_roi...; color = (roi_color, 0.17))
    CairoMakie.vlines!(ax, collect(energy_roi); color = (roi_color, 0.75),
                       linestyle = :dash, linewidth = 1.2)

    names = sort!(collect(keys(spectra)); by = name -> (-sum(spectra[name]), name))
    fallback_colors = CairoMakie.resample_cmap(:Dark2_8, max(2, length(names)))
    total = sum(values(spectra))
    total_error = sqrt.(sum(values(variances)))
    positive = filter(>(0), vcat(collect(values(spectra))...))
    ymin = isempty(positive) ? 1e-5 : minimum(positive) / 2
    ymax = isempty(positive) ? 1.0 : 2 * maximum(total + total_error)
    baseline = ymin / 2
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2

    function spectrum_line!(bins, errors; color, linewidth, linestyle = :solid, opacity = 0.30)
        populated = bins .> 0
        heights = bins[populated]
        sigma = errors[populated]
        CairoMakie.rangebars!(ax, centers[populated], max.(baseline, heights .- sigma),
                              heights .+ sigma; color = (color, opacity), linewidth = 1.2)
        x, y = spectrum_steps(edges, bins, baseline)
        return CairoMakie.stairs!(ax, x, y; step = :post, color, linewidth, linestyle)
    end
    components = [
        let (label, color, style) = get(BACKGROUND_STYLES, name, (name, fallback_colors[i], :dot))
            line = spectrum_line!(spectra[name], sqrt.(variances[name]);
                                  color, linestyle = style, linewidth = 2)
            (line, label)
        end for (i, name) in enumerate(names)
    ]
    total_line = spectrum_line!(total, total_error; color = ink, linewidth = 3, opacity = 0.65)
    CairoMakie.xlims!(ax, first(edges), last(edges))
    CairoMakie.ylims!(ax, ymin, ymax)
    CairoMakie.Legend(
        fig[1, 2],
        [total_line, roi_band, first.(components)...],
        ["Total background", "Energy ROI", last.(components)...],
        "Background components";
        framevisible = false, valign = :center, patchsize = (30, 14), labelsize = 17,
    )
    CairoMakie.save(path, fig; px_per_unit = 3)
    CairoMakie.save(splitext(path)[1] * ".pdf", fig)
end

function save_background_report(output_dir, phase, signal_process, files, tables, roi,
                                energy_prefilter, edges, best_bkg)
    live_time = signal_process.timeMeas
    # ⁸²Se exposure: foil mass × enrichment × live time in years.
    exposure = SNparams["foilMass"] * SNparams["a"] * live_time / (SNparams["t"] / SNparams["tYear"])
    spectra = Dict{String,Vector{Float64}}()
    variances = Dict{String,Vector{Float64}}()
    summary = DataFrame(
        process = String[], activity_bq = Union{Float64,Missing}[],
        mc_events_roi = Int[], mc_events_0_4_mev = Int[],
        expected_roi = Float64[], expected_0_4_mev = Float64[],
        mc_stat_error_roi = Float64[], mc_stat_error_0_4_mev = Float64[],
    )
    for (row, table) in zip(eachrow(files), tables)
        n_roi, bins = background_counts(table, roi, energy_prefilter, edges)
        weight = row.amount * row.activity_bq * live_time / row.nb_events
        name = String(row.name)
        n_full = sum(bins)
        push!(summary, (name, row.amount * row.activity_bq, n_roi, n_full, weight * n_roi, weight * n_full,
                        abs(weight) * sqrt(n_roi), abs(weight) * sqrt(n_full)))
        # As in ex5, combine PMTs in the plot but keep each process in the CSV.
        group = occursin("pmt", name) ? "pmt" : name
        get!(() -> zeros(length(bins)), spectra, group) .+= weight .* bins
        get!(() -> zeros(length(bins)), variances, group) .+= weight^2 .* bins
    end

    expected_roi = sum(summary.expected_roi)
    isapprox(expected_roi, best_bkg; atol = 1e-10) ||
        error("Recomputed ROI background ($expected_roi) does not match the optimization ($best_bkg)")
    # Activities of different sources are not summed. Base.missing: older
    # generated analysis scripts define a global `missing`.
    push!(summary, ("Total", Base.missing, sum(summary.mc_events_roi), sum(summary.mc_events_0_4_mev),
                    expected_roi, sum(summary.expected_0_4_mev),
                    sqrt(sum(abs2, summary.mc_stat_error_roi)),
                    sqrt(sum(abs2, summary.mc_stat_error_0_4_mev))))
    println("\n===== EXPECTED BACKGROUND (FULL EXPOSURE) =====")
    show(stdout, MIME"text/plain"(), summary; allrows = true, allcols = true)
    println()

    run_dir = dirname(output_dir)
    CSV.write(joinpath(run_dir, "background_components.csv"), summary)
    save_background_spectrum(joinpath(run_dir, "background_spectrum.png"), edges,
                             spectra, variances, roi[:energy_elec_sum], phase,
                             signal_process.isotopeName, exposure)
end
