include(joinpath(@__DIR__, "job_common.jl"))
using .SensitivityJobs
using Dates, SHA, Logging, CSV

function main(args)
    require(!isempty(args), "Usage: julia --project=/path/to/project scripts/run.jl CONFIG [--check | --run-dir DIR]")
    path = abspath(first(args))
    rest = args[2:end]
    require(isempty(rest) || rest == ["--check"] || (length(rest) == 2 && rest[1] == "--run-dir"), "Unknown run arguments")
    c = load_config(path)
    metadata = read_metadata(c)
    if rest == ["--check"]
        println("Preflight passed: $(size(metadata, 1)) processes; no ROOT events opened.")
        return
    end
    directory = if isempty(rest)
        joinpath(allocate_run(c), "results")
    else
        abspath(rest[2])
    end
    require(!ispath(directory), "Refusing to overwrite run directory: $directory")
    mkpath(directory)
    println("Run results/log: $directory")
    write_toml(joinpath(directory, "config.toml"), c)
    CSV.write(joinpath(directory, "metadata.csv"), metadata)
    provenance = Dict{String,Any}(
        "status" => "running", "started_utc" => string(now(UTC)),
        "julia_version" => string(VERSION), "threads" => Threads.nthreads(),
        "git_commit" => readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`),
        "git_dirty" => !isempty(readchomp(`git -C $PROJECT_ROOT status --porcelain`)),
        "manifest_sha256" => bytes2hex(sha256(read(joinpath(PROJECT_ROOT, "Manifest.toml")))),
        "metadata_sha256" => bytes2hex(sha256(read(c["data"]["metadata_csv"]))),
        "slurm_job_id" => get(ENV, "SLURM_JOB_ID", ""),
        "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", ""))
    write_toml(joinpath(directory, "run.toml"), provenance)
    open(joinpath(directory, "analysis.log"), "w") do io
        redirect_stdio(stdout=io, stderr=io) do
            with_logger(ConsoleLogger(io)) do
                try
                    include(joinpath(@__DIR__, "analysis.jl"))
                    summary = Base.invokelatest(SensitivityAnalysis.analyze, c, metadata, directory)
                    provenance["status"] = "completed"
                    provenance["sensitivity_package_version"] = string(Base.pkgversion(SensitivityAnalysis.SNSensitivityEstimate))
                    provenance["sensitivity_years"] = summary.sensitivity_years
                catch error
                    provenance["status"] = "failed"
                    provenance["error"] = sprint(showerror, error)
                    @error "Analysis failed; partial outputs are not a completed run" exception=(error, catch_backtrace())
                    rethrow()
                finally
                    provenance["finished_utc"] = string(now(UTC))
                    write_toml(joinpath(directory, "run.toml"), provenance)
                end
            end
        end
    end
    println("Completed: $directory")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
