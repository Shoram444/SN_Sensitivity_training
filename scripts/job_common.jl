module SensitivityJobs

using TOML, Dates, CSV, DataFramesMeta

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
export PROJECT_ROOT, load_config, cutset, variables, read_metadata, resolve_path,
    validate_config, write_toml, render_sbatch, allocate_run, require

require(condition, message) = condition || throw(ArgumentError(message))
resolve_path(path, root=PROJECT_ROOT) = normpath(isabspath(path) ? path : joinpath(root, path))
cutset(c) = c["cut_sets"][c["data"]["cut_set"]]
variables(c) = sort!(collect(keys(cutset(c)["variables"])))
finite_number(x) = x isa Real && !(x isa Bool) && isfinite(x)
positive(x) = finite_number(x) && x > 0

function check_keys(table, required, optional=String[]; context="configuration")
    require(table isa AbstractDict, "$context must be a table")
    require(isempty(setdiff(required, keys(table))), "$context is missing required keys: $(setdiff(required, keys(table)))")
    require(isempty(setdiff(keys(table), [required; optional])), "$context has unknown keys: $(setdiff(keys(table), [required; optional]))")
end

function interval(x, context)
    require(x isa AbstractVector && length(x) == 2 && all(finite_number, x), "$context needs two finite numbers")
    require(x[1] < x[2], "$context lower bound must be smaller than upper bound")
end

on_grid(x, lo, step) = isapprox((x-lo)/step, round((x-lo)/step); atol=1e-9, rtol=1e-10)

function slurm_seconds(value)
    require(value isa String && occursin(r"^(\d+-)?\d+:\d{2}:\d{2}$", value), "SLURM time must be [days-]hours:minutes:seconds")
    pieces = split(value, '-')
    days = length(pieces) == 2 ? parse(Int, pieces[1]) : 0
    h, m, s = parse.(Int, split(pieces[end], ':'))
    require(m < 60 && s < 60 && (length(pieces) == 1 || h < 24), "Invalid SLURM time")
    return days*86400 + h*3600 + m*60 + s
end

function validate_config(c)
    check_keys(c, ["slurm", "runtime", "data", "optimization", "cut_sets", "outputs", "validation"])
    s, r, d, o, out, val = (c[k] for k in ("slurm", "runtime", "data", "optimization", "outputs", "validation"))
    check_keys(s, ["job_name", "partition", "licenses", "time", "mem", "cpus_per_task"])
    for k in ("job_name", "partition", "licenses", "mem")
        require(s[k] isa String && occursin(r"^[A-Za-z0-9_.,:/@+=-]+$", s[k]), "Invalid SLURM $k")
    end
    require(s["cpus_per_task"] isa Integer && !(s["cpus_per_task"] isa Bool) && s["cpus_per_task"] > 0, "cpus_per_task must be a positive integer")
    check_keys(r, ["julia_module", "julia_executable", "cpu_target", "sysimage", "overhead_seconds"])
    require(all(r[k] isa String for k in ("julia_module", "julia_executable", "cpu_target", "sysimage")), "Runtime commands/paths must be strings")
    require(positive(r["overhead_seconds"]), "overhead_seconds must be positive")
    require(!isempty(r["julia_executable"]) && !isempty(r["cpu_target"]), "Julia executable and CPU target cannot be empty")
    check_keys(d, ["metadata_csv", "phase", "signal", "exposure_seconds", "tree", "cut_set", "activity_source"])
    require(all(d[k] isa String && !isempty(d[k]) for k in ("metadata_csv", "phase", "signal", "tree", "cut_set")), "Data paths/names must be nonempty strings")
    require(positive(d["exposure_seconds"]), "Simulation exposure must be positive seconds")
    require(d["activity_source"] == "metadata", "Only reviewed metadata activities are supported")
    check_keys(o, ["algorithm", "time_limit_seconds", "seed", "alpha", "search_method", "final_method"])
    require(o["algorithm"] in ("PSO", "ECA", "SA"), "Supported algorithms: PSO, ECA, SA")
    require(positive(o["time_limit_seconds"]) && positive(o["alpha"]), "Optimization time and alpha must be positive")
    require(o["seed"] isa Integer && !(o["seed"] isa Bool) && 0 <= o["seed"] <= typemax(Int32), "Seed must be an integer in 0:$(typemax(Int32))")
    require(all(o[k] in ("formula", "table") for k in ("search_method", "final_method")), "Use package formula/table FC methods")
    require(slurm_seconds(s["time"]) > o["time_limit_seconds"] + r["overhead_seconds"], "SLURM time must exceed optimization time plus overhead")
    require(haskey(c["cut_sets"], d["cut_set"]), "Unknown cut set $(d["cut_set"])")
    for (name, cs) in c["cut_sets"]
        check_keys(cs, ["value_mode", "fixed", "variables"]; context="cut set $name")
        require(cs["value_mode"] == "absolute", "ND variables must use absolute values, consistently with the package")
        f = cs["fixed"]
        check_keys(f, ["vertex_distance_lt_mm", "reject_same_calo", "closest_gamma_abs_gt_mm",
            "closest_electron_abs_gt_mm", "closest_track_abs_ge_mm", "closest_track_time_abs_ge_ns",
            "reject_kinks", "required_om_count", "om_max_id", "reject_vertex_hotspots", "vertex_hotspots",
            "vertex_abs_y_le_mm", "vertex_abs_z_le_mm", "om_quality_enabled", "om_quality_file"])
        for key in ("reject_same_calo", "reject_kinks", "reject_vertex_hotspots", "om_quality_enabled")
            require(f[key] isa Bool, "$name.$key must be boolean")
        end
        for key in ("vertex_distance_lt_mm", "vertex_abs_y_le_mm", "vertex_abs_z_le_mm")
            require(positive(f[key]), "$name.$key must be positive")
        end
        for key in ("closest_gamma_abs_gt_mm", "closest_electron_abs_gt_mm", "closest_track_abs_ge_mm", "closest_track_time_abs_ge_ns")
            require(finite_number(f[key]) && f[key] >= 0, "$name.$key must be nonnegative")
        end
        require(f["required_om_count"] == 2 && !(f["required_om_count"] isa Bool), "Two-electron selection requires exactly two OMs")
        require(f["om_max_id"] isa Integer && !(f["om_max_id"] isa Bool) && f["om_max_id"] >= 0, "om_max_id must be nonnegative")
        require(f["om_quality_file"] isa String, "om_quality_file must be a path")
        require(!f["om_quality_enabled"] || !isempty(f["om_quality_file"]), "OM quality requires a phase-specific list")
        require(f["vertex_hotspots"] isa AbstractVector, "vertex_hotspots must be an array")
        for h in f["vertex_hotspots"]
            check_keys(h, ["y_mm", "z_mm", "radius_y_mm", "radius_z_mm"])
            require(all(finite_number, values(h)) && positive(h["radius_y_mm"]) && positive(h["radius_z_mm"]), "Invalid hotspot ellipse")
        end
        require(!isempty(cs["variables"]), "A cut set must contain variables")
        for (v, spec) in cs["variables"]
            check_keys(spec, ["role", "bounds", "unit"], ["step", "initial_roi"]; context="variable $v")
            require(spec["role"] in ("fixed", "optimize"), "Variable $v role must be fixed or optimize")
            require(spec["unit"] isa String && !isempty(spec["unit"]), "Variable $v needs units")
            interval(spec["bounds"], "$v bounds")
            require(spec["bounds"][1] >= 0, "Absolute variable $v cannot have negative bounds")
            if spec["role"] == "optimize"
                require(haskey(spec, "step") && positive(spec["step"]) && haskey(spec, "initial_roi"), "$v needs a positive step and initial ROI")
                interval(spec["initial_roi"], "$v initial ROI")
                lo, hi = spec["bounds"]; a, b = spec["initial_roi"]
                require(lo <= a < b <= hi, "$v initial ROI is outside bounds")
                require(all(x -> on_grid(x, lo, spec["step"]), [hi, a, b]), "$v bounds/initial ROI must align with its grid")
            end
        end
    end
    check_keys(out, ["root", "tables", "plots", "plot_formats", "spectrum_variable", "spectrum_bounds", "spectrum_bin_width", "text_summary", "background_groups"])
    require(out["root"] isa String && !isempty(out["root"]), "Output root must be a path")
    for (key, allowed) in (("tables", ["sensitivity", "roi", "background_components", "background_groups", "cutflow"]),
                          ("plots", ["background_spectrum", "sensitivity_slices"]), ("plot_formats", ["pdf", "png"]))
        require(out[key] isa AbstractVector && all(x -> x in allowed, out[key]) && allunique(out[key]), "Invalid or duplicate output $key")
    end
    require(!isempty(out["plot_formats"]) || isempty(out["plots"]), "Plots need at least one format")
    require(out["text_summary"] isa Bool, "text_summary must be boolean")
    require(haskey(cutset(c)["variables"], out["spectrum_variable"]), "Spectrum variable must be in the selected cut set")
    interval(out["spectrum_bounds"], "Spectrum bounds")
    lo, hi = out["spectrum_bounds"]; domain = cutset(c)["variables"][out["spectrum_variable"]]["bounds"]
    require(domain[1] <= lo < hi <= domain[2], "Spectrum bounds must be inside the loaded variable domain")
    require(positive(out["spectrum_bin_width"]) && on_grid(hi, lo, out["spectrum_bin_width"]), "Spectrum bounds must align with bin width")
    require(all(v isa AbstractVector && all(x -> x isa String, v) for v in values(out["background_groups"])), "Background groups must list process names")
    members = String[x for list in values(out["background_groups"]) for x in list]
    require(allunique(members), "A background may belong to only one group")
    check_keys(val, ["enabled", "measured_data_files", "measured_exposure_seconds"])
    require(val["enabled"] isa Bool && val["measured_data_files"] isa AbstractVector && all(x -> x isa String, val["measured_data_files"]), "Invalid measured-data configuration")
    require(positive(val["measured_exposure_seconds"]), "Measured exposure must be positive seconds")
    require(!val["enabled"] || !isempty(val["measured_data_files"]), "Data validation requires measured files")
    return c
end

function load_config(path; root=PROJECT_ROOT)
    c = validate_config(TOML.parsefile(path))
    for (section, key) in (("data", "metadata_csv"), ("outputs", "root"))
        c[section][key] = resolve_path(c[section][key], root)
    end
    if !isempty(c["runtime"]["sysimage"])
        c["runtime"]["sysimage"] = resolve_path(c["runtime"]["sysimage"], root)
    end
    for cs in values(c["cut_sets"])
        f = cs["fixed"]
        isempty(f["om_quality_file"]) || (f["om_quality_file"] = resolve_path(f["om_quality_file"], root))
    end
    c["validation"]["measured_data_files"] = resolve_path.(c["validation"]["measured_data_files"], Ref(root))
    return c
end

function read_metadata(c; check_files=true, root=PROJECT_ROOT)
    path = c["data"]["metadata_csv"]
    require(isfile(path), "Analysis CSV missing: $path. Prepare the ROOT-derived CSV first; no fallback to unverified totals.")
    df = CSV.read(path, DataFrame)
    required = ["phase", "name", "is_signal", "file_path", "activity_bq", "amount", "nb_events", "nb_events_source"]
    require(all(k -> k in names(df), required), "Analysis CSV needs columns: $(join(required, ", "))")
    require(!any(ismissing, Matrix(df[:, required])), "Analysis CSV has missing required values")
    require(all(x -> x isa Bool, df.is_signal), "is_signal must contain true/false")
    chosen = [r.phase == c["data"]["phase"] && (!r.is_signal || r.name == c["data"]["signal"]) for r in eachrow(df)]
    df = df[chosen, :]
    require(count(identity, df.is_signal) == 1, "Selected phase must contain exactly one requested signal")
    require(any(!, df.is_signal), "At least one background is required")
    require(allunique(df.name), "Duplicate process rows in the selected phase")
    for r in eachrow(df)
        require(r.name isa AbstractString && r.file_path isa AbstractString, "Process names and ROOT paths must be strings")
        require(positive(r.nb_events) && isinteger(r.nb_events), "Invalid generated event total for $(r.name)")
        require(positive(r.amount) && finite_number(r.activity_bq) && r.activity_bq >= 0, "Invalid normalization for $(r.name)")
        require(r.nb_events_source == "ROOT:last_event_number", "Unverified generated count for $(r.name)")
    end
    df.file_path = resolve_path.(String.(df.file_path), Ref(root))
    require(allunique(df.file_path), "Multiple selected processes reference the same ROOT file")
    members = [x for list in values(c["outputs"]["background_groups"]) for x in list]
    require(Set(members) == Set(String.(df.name[.!df.is_signal])), "Background groups must cover selected backgrounds exactly once")
    if check_files
        for p in [df.file_path; c["validation"]["enabled"] ? c["validation"]["measured_data_files"] : String[]]
            require(isfile(p), "Missing ROOT input: $p")
        end
        fingerprints = ["root_size_bytes", "root_mtime_unix"]
        if any(k -> k in names(df), fingerprints)
            require(all(k -> k in names(df), fingerprints), "Incomplete ROOT file fingerprints in metadata")
            for row in eachrow(df)
                info = stat(row.file_path)
                require(info.size == row.root_size_bytes && info.mtime == row.root_mtime_unix,
                    "ROOT input changed since count extraction: $(row.file_path). Regenerate analysis metadata.")
            end
        end
    end
    f = cutset(c)["fixed"]
    require(!f["om_quality_enabled"] || isfile(f["om_quality_file"]), "Missing OM quality file")
    require(isempty(c["runtime"]["sysimage"]) || isfile(c["runtime"]["sysimage"]), "Missing configured sysimage")
    return df
end

function write_toml(path, table)
    open(path, "w") do io
        TOML.print(io, table; sorted=true)
    end
end

function allocate_run(c)
    mkpath(c["outputs"]["root"])
    return mktempdir(c["outputs"]["root"]; prefix="run-$(Dates.format(now(UTC), "yyyymmddTHHMMSS"))-", cleanup=false)
end

shellquote(s) = "'" * replace(string(s), "'" => "'\"'\"'") * "'"

function render_sbatch(c, config_path, run_dir; root=PROJECT_ROOT)
    s, runtime = c["slurm"], c["runtime"]
    lines = ["#!/bin/bash"]
    for key in ("job_name", "partition", "licenses", "time", "mem", "cpus_per_task")
        push!(lines, "#SBATCH --$(replace(key, "_" => "-"))=$(s[key])")
    end
    push!(lines, "#SBATCH --output=$(repr(joinpath(run_dir, "slurm-%j.log")))")
    append!(lines, ["", "set -euo pipefail", "cd $(shellquote(root))",
        "export JULIA_NUM_THREADS=$(s["cpus_per_task"])",
        "export JULIA_CPU_TARGET=$(shellquote(runtime["cpu_target"]))"])
    isempty(runtime["julia_module"]) || push!(lines, "module load $(shellquote(runtime["julia_module"]))")
    args = [runtime["julia_executable"], "--project=$root", "--startup-file=no",
        "--threads=$(s["cpus_per_task"])", "--pkgimages=no"]
    isempty(runtime["sysimage"]) || push!(args, "--sysimage=$(runtime["sysimage"])")
    append!(args, [joinpath(root, "scripts", "run.jl"), config_path, "--run-dir", joinpath(run_dir, "results")])
    push!(lines, "exec " * join(shellquote.(args), " "))
    return join(lines, "\n") * "\n"
end

end
