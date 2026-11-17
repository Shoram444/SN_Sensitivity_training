include(joinpath(@__DIR__, "job_common.jl"))
using .SensitivityJobs

function main(args)
    require(length(args) in (1, 2) && (length(args) == 1 || args[2] == "--dry-run"),
        "Usage: julia --project=/path/to/project scripts/submit.jl CONFIG [--dry-run]")
    c = load_config(abspath(args[1]))
    if length(args) == 2
        directory = joinpath(c["outputs"]["root"], "DRY-RUN")
        println("# Configuration-only preview: no files written, metadata/ROOT inputs not checked.")
        print(render_sbatch(c, joinpath(directory, "config.toml"), directory))
        return
    end
    read_metadata(c)
    directory = allocate_run(c)
    config_path = joinpath(directory, "config.toml")
    write_toml(config_path, c)
    script = joinpath(directory, "job.sbatch")
    write(script, render_sbatch(c, config_path, directory))
    println("Submitting: $script")
    run(`sbatch $script`)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
