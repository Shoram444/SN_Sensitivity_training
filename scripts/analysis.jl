module SensitivityAnalysis

using ..SensitivityJobs
using SNSensitivityEstimate, UnROOT, Metaheuristics, CSV, DataFramesMeta, CairoMakie

export analyze

function normalization(c)
    params = copy(SNparams)
    seconds_per_year = SNparams["t"] / SNparams["tYear"]
    params["t"] = Float64(c["data"]["exposure_seconds"])
    params["tYear"] = params["t"] / seconds_per_year
    return params
end

function in_bounds(event, name, bounds)
    value = getproperty(event, Symbol(name))
    require(value isa Real && isfinite(value), "Nonfinite/nonnumeric value in branch $name")
    return bounds[1] <= abs(value) < bounds[2]
end

in_roi(event, roi; except=nothing) = all(k == except || in_bounds(event, String(k), roi[k]) for k in keys(roi))

function vertex(event)
    y, z = event.vertex_3D_start_y, event.vertex_3D_start_z
    require(length(y) >= 2 && length(z) >= 2, "Selected two-electron event has incomplete vertex coordinates")
    return (0.5*(y[1]+y[2]), 0.5*(z[1]+z[2]))
end

function excluded_oms(f)
    excluded = Set{Int}()
    if f["om_quality_enabled"]
        for line in readlines(f["om_quality_file"])[2:end], token in split(line)
            id = parse(Int, token)
            require(id >= 0, "Negative OM identifier in quality list")
            push!(excluded, id)
        end
    end
    return excluded
end

function selection_stages(c)
    f = cutset(c)["fixed"]; bad_oms = excluded_oms(f)
    stages = Pair{String,Function}[
        "two_valid_oms" => e -> length(e.num_om_elec_f) == f["required_om_count"] &&
            all(id -> id isa Integer && 0 <= id <= f["om_max_id"], e.num_om_elec_f),
        "vertex_distance" => e -> hypot(e.delta_y_elec, e.delta_z_elec) < f["vertex_distance_lt_mm"],
        "gamma_veto" => e -> abs(e.closest_gamma) > f["closest_gamma_abs_gt_mm"],
        "electron_veto" => e -> abs(e.closest_elec) > f["closest_electron_abs_gt_mm"],
        "track_veto" => e -> abs(e.closest_track) >= f["closest_track_abs_ge_mm"],
        "track_time_veto" => e -> abs(e.closest_time_track) >= f["closest_track_time_abs_ge_ns"],
    ]
    f["reject_same_calo"] && push!(stages, "different_calo" => e -> e.hit_the_same_calo_hit != 1)
    f["reject_kinks"] && push!(stages, "no_kinks" => e -> e.has_kinks != 1)
    f["om_quality_enabled"] && push!(stages, "om_quality" => e -> all(id -> !(id in bad_oms), e.num_om_elec_f))
    push!(stages, "vertex_region" => function(e)
        y, z = vertex(e)
        return abs(y) <= f["vertex_abs_y_le_mm"] && abs(z) <= f["vertex_abs_z_le_mm"]
    end)
    if f["reject_vertex_hotspots"]
        push!(stages, "vertex_hotspots" => function(e)
            y, z = vertex(e)
            return !any(((y-h["y_mm"])/h["radius_y_mm"])^2 +
                        ((z-h["z_mm"])/h["radius_z_mm"])^2 <= 1 for h in f["vertex_hotspots"])
        end)
    end
    for name in variables(c)
        spec = cutset(c)["variables"][name]
        push!(stages, "$(spec["role"])_domain:$name" => e -> in_bounds(e, name, spec["bounds"]))
    end
    return stages
end

function branch_names(c)
    names = [variables(c); ["num_om_elec_f", "delta_y_elec", "delta_z_elec", "closest_gamma",
        "closest_elec", "closest_track", "closest_time_track", "vertex_3D_start_y", "vertex_3D_start_z"]]
    f = cutset(c)["fixed"]
    f["reject_same_calo"] && push!(names, "hit_the_same_calo_hit")
    f["reject_kinks"] && push!(names, "has_kinks")
    return unique(names)
end

function preselect(tree, c)
    stages = selection_stages(c)
    counts = zeros(Int, length(stages)+1)
    keep = trues(length(tree))
    for (i, event) in enumerate(tree)
        counts[1] += 1
        for (j, (_, pass)) in enumerate(stages)
            if !pass(event)
                keep[i] = false
                break
            end
            counts[j+1] += 1
        end
    end
    return tree[keep], ["input"; first.(stages)], counts
end

function open_tree(path, c)
    file = UnROOT.ROOTFile(path)
    branches = branch_names(c)
    tree = UnROOT.LazyTree(file, c["data"]["tree"], branches)
    require(all(n -> Symbol(n) in propertynames(tree), branches), "Missing requested branches in $path")
    return file, tree
end

function build_processes(c, metadata; loader=open_tree)
    names = variables(c)
    bins = NamedTuple(Symbol(n) => Tuple(Float64.(cutset(c)["variables"][n]["bounds"])) for n in names)
    processes = DataProcessND[]
    files = Any[]
    cutflow = NamedTuple[]
    for row in eachrow(metadata)
        @info "Loading process" name=row.name file=row.file_path
        file, tree = loader(row.file_path, c)
        push!(files, file)
        selected, stages, counts = preselect(tree, c)
        for (stage, passed) in zip(stages, counts)
            push!(cutflow, (process=String(row.name), is_signal=Bool(row.is_signal), stage=stage,
                mc_passed=passed, generated=Float64(row.nb_events), efficiency=passed/row.nb_events))
        end
        push!(processes, DataProcessND(selected, String(row.name), Bool(row.is_signal),
            Float64(row.activity_bq), Float64(c["data"]["exposure_seconds"]), Float64(row.nb_events),
            bins, Float64(row.amount), names))
    end
    return (; processes, files, cutflow)
end

function optimization_grid(c)
    specs = cutset(c)["variables"]
    names = filter(n -> specs[n]["role"] == "optimize", variables(c))
    # Config validation ensures integral spans; rounding avoids floor losing a node to floating-point error.
    return GridSpec(Symbol.(names), [Float64(specs[n]["bounds"][1]) for n in names],
        [Float64(specs[n]["step"]) for n in names],
        [round(Int, (specs[n]["bounds"][2]-specs[n]["bounds"][1])/specs[n]["step"]) for n in names])
end

function decode_roi(c, grid, indices)
    edges = decode_grid_roi(indices, grid)
    optimized = Dict(k => (edges[2i-1], edges[2i]) for (i, k) in enumerate(grid.keys))
    return NamedTuple(Symbol(n) => get(optimized, Symbol(n), Tuple(Float64.(cutset(c)["variables"][n]["bounds"]))) for n in variables(c))
end

function optimize_roi(c, processes, params)
    grid = optimization_grid(c)
    if isempty(grid.keys)
        return decode_roi(c, grid, Float64[]), nothing
    end
    specs = cutset(c)["variables"]
    x0 = Float64[index for (i, k) in enumerate(grid.keys) for index in
        (value_to_index(v, grid.mins[i], grid.steps[i]; nsteps=grid.nsteps[i]) for v in specs[String(k)]["initial_roi"])]
    o = c["optimization"]
    objective = function(indices)
        roi = decode_roi(c, grid, indices)
        edges = Float64[x for bounds in values(roi) for x in bounds]
        return -get_s_to_b(params, Float64(o["alpha"]), processes, edges; approximate=o["search_method"])
    end
    options = Options(time_limit=Float64(o["time_limit_seconds"]), seed=o["seed"], verbose=true)
    algorithm = getproperty(Metaheuristics, Symbol(o["algorithm"]))(; options)
    set_user_solutions!(algorithm, x0, objective)
    ranges = make_stepRange(grid)
    bounds = boxconstraints(lb=Float64[first(r) for r in ranges], ub=Float64[last(r) for r in ranges])
    result = optimize(objective, bounds, algorithm)
    roi = decode_roi(c, grid, minimizer(result))
    require(all(v[1] < v[2] for v in values(roi)), "Optimizer did not find a valid ROI")
    return roi, result
end

function component_table(c, processes, roi)
    rows = NamedTuple[]
    for p in processes
        p.signal && continue
        n = count(e -> in_roi(e, roi), p.data)
        weight = p.activity*p.amount*p.timeMeas/p.nTotalSim
        group = only(k for (k, members) in c["outputs"]["background_groups"] if p.isotopeName in members)
        push!(rows, (process=p.isotopeName, group=group, activity_bq=Float64(p.activity),
            amount=p.amount, exposure_seconds=p.timeMeas, generated=p.nTotalSim,
            mc_before_roi=length(p.data), mc_in_roi=n, efficiency=n/p.nTotalSim,
            event_weight=weight, expected_counts=n*weight, mc_stat_error=sqrt(n)*weight))
    end
    return DataFrame(rows)
end

function group_table(components)
    return DataFrame([(group=String(g), mc_in_roi=sum(s.mc_in_roi),
        expected_counts=sum(s.expected_counts), mc_stat_error=sqrt(sum(abs2, s.mc_stat_error)))
        for g in sort(unique(components.group)) for s in (components[components.group .== g, :],)])
end

function spectrum(events, roi, variable, edges, weight)
    counts = zeros(length(edges)-1)
    for e in events
        in_roi(e, roi; except=variable) || continue
        value = abs(getproperty(e, variable))
        i = searchsortedlast(edges, value)
        1 <= i < length(edges) && (counts[i] += weight)
    end
    return counts
end

function save_plot(figure, c, directory, name)
    for format in c["outputs"]["plot_formats"]
        save(joinpath(directory, "$name.$format"), figure)
    end
end

function background_spectrum(c, processes, roi, fc_limit, directory; measured=UnROOT.LazyEvent[])
    out = c["outputs"]; variable = Symbol(out["spectrum_variable"])
    lo, hi = out["spectrum_bounds"]
    edges = collect(range(lo, hi; length=round(Int, (hi-lo)/out["spectrum_bin_width"])+1))
    groups = Dict(k => zeros(length(edges)-1) for k in keys(out["background_groups"]))
    variance = zeros(length(edges)-1)
    signal = only(p for p in processes if p.signal)
    signal_weight = fc_limit/count(e -> in_roi(e, roi), signal.data)
    sig = spectrum(signal.data, roi, variable, edges, signal_weight)
    for p in processes
        p.signal && continue
        group = only(k for (k, members) in out["background_groups"] if p.isotopeName in members)
        weight = p.activity*p.amount*p.timeMeas/p.nTotalSim
        groups[group] .+= spectrum(p.data, roi, variable, edges, weight)
        variance .+= spectrum(p.data, roi, variable, edges, weight^2)
    end
    total = sum(values(groups))
    df = DataFrame(lower=edges[1:end-1], upper=edges[2:end], expected_counts=total,
        mc_stat_error=sqrt.(variance), signal_counts_at_limit=sig)
    for group in sort(collect(keys(groups)))
        df[!, Symbol(group)] = groups[group]
    end
    CSV.write(joinpath(directory, "background_spectrum.csv"), df)
    if c["validation"]["enabled"]
        ndata = spectrum(measured, roi, variable, edges, 1.0)
        sim_days = c["data"]["exposure_seconds"]/86400
        data_days = c["validation"]["measured_exposure_seconds"]/86400
        comparison = DataFrame(lower=edges[1:end-1], upper=edges[2:end],
            simulation_counts_per_day=total./sim_days, simulation_error_per_day=sqrt.(variance)./sim_days,
            data_counts=ndata, data_counts_per_day=ndata./data_days, data_error_per_day=sqrt.(ndata)./data_days)
        CSV.write(joinpath(directory, "data_vs_mc.csv"), comparison)
    end
    "background_spectrum" in out["plots"] || return
    fig = Figure(size=(1400, 650))
    for (i, scale) in enumerate((identity, log10))
        ax = Axis(fig[1, i], xlabel="$(out["spectrum_variable"]) [$(cutset(c)["variables"][out["spectrum_variable"]]["unit"])]",
            ylabel="Expected counts / bin", yscale=scale)
        display_counts(x) = scale === log10 ? max.(x, 1e-10) : x
        for group in sort(collect(keys(groups)))
            v = display_counts(groups[group])
            stairs!(ax, edges, [v; last(v)]; step=:post, label=group)
        end
        for (label, values, style) in (("Total background", total, :solid), ("Signal at sensitivity", sig, :dash))
            v = display_counts(values)
            stairs!(ax, edges, [v; last(v)]; step=:post, color=:black, linestyle=style, label=label)
        end
        vspan!(ax, roi[variable]...; color=(:red, 0.08))
        i == 1 && axislegend(ax; position=:rt, labelsize=12)
    end
    save_plot(fig, c, directory, "background_spectrum")
    if c["validation"]["enabled"]
        data = CSV.read(joinpath(directory, "data_vs_mc.csv"), DataFrame)
        f = Figure(); ax = Axis(f[1, 1], xlabel=out["spectrum_variable"], ylabel="Counts / day / bin")
        center = (data.lower .+ data.upper)./2
        scatter!(ax, center, data.data_counts_per_day; label="Measured phase 3")
        errorbars!(ax, center, data.data_counts_per_day, data.data_error_per_day)
        lines!(ax, center, data.simulation_counts_per_day; label="Simulation, phase-3 conditions")
        axislegend(ax)
        save_plot(f, c, directory, "data_vs_mc")
    end
end

function sensitivity_slices(c, processes, params, roi, directory)
    rows = NamedTuple[]
    for name in variables(c)
        spec = cutset(c)["variables"][name]
        spec["role"] == "optimize" || continue
        key = Symbol(name); lower, upper = roi[key]
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel="$name [$(spec["unit"])]", ylabel="Half-life sensitivity [yr]")
        n = round(Int, (spec["bounds"][2]-spec["bounds"][1])/spec["step"])
        for edge in ("lower", "upper")
            x = Float64[]; y = Float64[]
            for index in 0:n
                v = spec["bounds"][1]+index*spec["step"]
                a, b = edge == "lower" ? (v, upper) : (lower, v)
                a < b || continue
                candidate = merge(roi, NamedTuple{(key,)}(((a, b),)))
                result = get_sensitivityND(params, c["optimization"]["alpha"], processes, candidate;
                    approximate=c["optimization"]["final_method"])
                background = get_roi_bkg_counts(processes, candidate)
                push!(rows, (variable=name, edge=edge, value=v, unit=spec["unit"],
                    sensitivity_years=Float64(result.tHalf), background_counts=Float64(background),
                    signal_efficiency=Float64(result.signalEff)))
                push!(x, v); push!(y, result.tHalf)
            end
            lines!(ax, x, y; label="$edge edge; other edges fixed")
        end
        axislegend(ax)
        save_plot(fig, c, directory, "sensitivity_$(name)")
    end
    isempty(rows) || CSV.write(joinpath(directory, "sensitivity_slices.csv"), DataFrame(rows))
end

function analyze(c, metadata, directory; loader=open_tree)
    params = normalization(c)
    context = build_processes(c, metadata; loader)
    processes = context.processes
    @info "Starting optimization" exposure_seconds=params["t"] years=params["tYear"]
    roi, result = optimize_roi(c, processes, params)
    o = c["optimization"]
    best = get_sensitivityND(params, o["alpha"], processes, roi; approximate=o["final_method"])
    require(isfinite(best.tHalf) && best.tHalf > 0 && best.signalEff > 0, "No positive finite sensitivity found")
    components = component_table(c, processes, roi)
    require(isapprox(sum(components.expected_counts), best.bkgCounts; rtol=1e-10, atol=1e-10), "Reported background counts disagree with package ROI counts")
    signal = only(p for p in processes if p.signal)
    fc_limit = Float64(get_FC(best.bkgCounts, o["alpha"]; approximate=o["final_method"]))
    measured = UnROOT.LazyEvent[]
    if c["validation"]["enabled"]
        for path in c["validation"]["measured_data_files"]
            file, tree = loader(path, c)
            push!(context.files, file)
            selected, _, _ = preselect(tree, c)
            append!(measured, selected)
        end
    end
    summary = (signal=c["data"]["signal"], phase=c["data"]["phase"], exposure_seconds=params["t"],
        exposure_years=params["tYear"], isotope_exposure_kg_year=params["foilMass"]*params["a"]*params["tYear"],
        sensitivity_years=Float64(best.tHalf), background_counts=Float64(best.bkgCounts),
        background_mc_stat_error=sqrt(sum(abs2, components.mc_stat_error)), signal_efficiency=Float64(best.signalEff),
        signal_mc_in_roi=count(e -> in_roi(e, roi), signal.data), signal_generated=signal.nTotalSim,
        fc_signal_limit=fc_limit, alpha=o["alpha"], final_method=o["final_method"], search_method=o["search_method"],
        seed=o["seed"], algorithm=o["algorithm"], iterations=result === nothing ? 0 : result.iteration,
        evaluations=result === nothing ? 1 : result.f_calls, optimizer_seconds=result === nothing ? 0.0 : result.overall_time,
        stop_reason=result === nothing ? "fixed_cuts_only" : string(result.termination_status_code),
        validation_status=c["validation"]["enabled"] ? "comparison_produced" : "disabled")
    for p in processes
        n = count(e -> in_roi(e, roi), p.data)
        push!(context.cutflow, (process=p.isotopeName, is_signal=p.signal, stage="final_roi",
            mc_passed=n, generated=p.nTotalSim, efficiency=n/p.nTotalSim))
    end
    tables = c["outputs"]["tables"]
    "sensitivity" in tables && CSV.write(joinpath(directory, "sensitivity.csv"), DataFrame([summary]))
    "roi" in tables && CSV.write(joinpath(directory, "roi.csv"), DataFrame([
        (variable=n, lower=roi[Symbol(n)][1], upper=roi[Symbol(n)][2],
            unit=cutset(c)["variables"][n]["unit"], value_mode="absolute",
            role=cutset(c)["variables"][n]["role"]) for n in variables(c)]))
    "background_components" in tables && CSV.write(joinpath(directory, "background_components.csv"), components)
    "background_groups" in tables && CSV.write(joinpath(directory, "background_groups.csv"), group_table(components))
    "cutflow" in tables && CSV.write(joinpath(directory, "cutflow.csv"), DataFrame(context.cutflow))
    if c["outputs"]["text_summary"]
        open(joinpath(directory, "summary.txt"), "w") do io
            for (k, value) in pairs(summary)
                println(io, "$k = $value")
            end
            println(io, "absolute-value ROI (lower inclusive, upper exclusive): ", roi)
        end
    end
    write_toml(joinpath(directory, "physics_parameters.toml"), params)
    if "background_spectrum" in c["outputs"]["plots"] || c["validation"]["enabled"]
        background_spectrum(c, processes, roi, fc_limit, directory; measured)
    end
    "sensitivity_slices" in c["outputs"]["plots"] && sensitivity_slices(c, processes, params, roi, directory)
    @info "Analysis complete" sensitivity_years=best.tHalf background_counts=best.bkgCounts
    return summary
end

end
