include(joinpath(@__DIR__, "job_common.jl"))
using .SensitivityJobs
include(joinpath(@__DIR__, "metadata.jl"))
using .SensitivityMetadata

function main(args)
    require(!isempty(args), "Usage: julia --project=/path/to/project scripts/prepare_metadata.jl CONFIG [--source CSV] [--dry-run] [--overwrite]")
    c = load_config(abspath(first(args)))
    source = joinpath(PROJECT_ROOT, "data", "mathis_data_info_final.csv")
    dry_run = false; overwrite = false; i = 2
    while i <= length(args)
        if args[i] == "--source"
            require(i < length(args), "--source requires a CSV path")
            source = abspath(args[i+1]); i += 2
        elseif args[i] == "--dry-run"
            dry_run = true; i += 1
        elseif args[i] == "--overwrite"
            overwrite = true; i += 1
        else
            throw(ArgumentError("Unknown metadata argument: $(args[i])"))
        end
    end
    if dry_run
        data, paths = select_source(c, source)
        println("Metadata-only extraction preview: $(size(data, 1)) processes in $(c["data"]["phase"]).")
        println("Output: ", c["data"]["metadata_csv"])
        println("No files written or ROOT files opened. File availability and counts are not verified.")
        foreach(println, paths)
        return
    end
    data = prepare_metadata(c, source; overwrite)
    println("Wrote $(size(data, 1)) verified generated totals to $(c["data"]["metadata_csv"]).")
    println("Activities unchanged. last_event_number copied exactly; no +1 adjustment or tree-entry substitution.")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
