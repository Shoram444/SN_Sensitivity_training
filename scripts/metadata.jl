module SensitivityMetadata

using ..SensitivityJobs
using CSV, DataFramesMeta, Dates, SHA

export select_source, read_root_counts, prepare_metadata, write_analysis_csv

function select_source(c, source; root=PROJECT_ROOT)
    df = CSV.read(source, DataFrame)
    required = ["phase", "name", "is_signal", "file_path", "activity_bq", "amount", "nb_events"]
    require(all(n -> n in names(df), required), "Source CSV is missing required columns")
    require(!any(ismissing, Matrix(df[:, required])), "Source CSV has missing required values")
    require(all(x -> x isa Bool, df.is_signal), "Source is_signal must contain true/false")
    keep = [r.phase == c["data"]["phase"] &&
        (!r.is_signal || r.name == c["data"]["signal"]) &&
        !("reduces_files" in splitpath(r.file_path)) for r in eachrow(df)]
    df = df[keep, :]
    require(count(identity, df.is_signal) == 1, "Source must contain exactly one selected signal")
    require(allunique(df.name), "Duplicate selected source processes after excluding reduces_files")
    members = [x for list in values(c["outputs"]["background_groups"]) for x in list]
    require(Set(members) == Set(String.(df.name[.!df.is_signal])), "Configured backgrounds do not match selected source rows")
    paths = resolve_path.(String.(df.file_path), Ref(root))
    require(allunique(paths), "Duplicate selected ROOT paths")
    require(all(!any(ch -> ch in ('\n', '\r'), path) for path in paths), "ROOT paths must not contain newlines")
    return df, paths
end

function parse_counts(path, count)
    require(isfile(path), "ROOT did not produce its metadata report")
    df = CSV.read(path, DataFrame; delim='\t', types=Dict("file_index" => Int, "last_event_number" => Int64), strict=true)
    require(names(df) == ["file_index", "last_event_number"], "Unexpected ROOT metadata report columns")
    require(!any(ismissing, Matrix(df)), "ROOT metadata report contains missing values")
    require(size(df, 1) == count && allunique(df.file_index) &&
        Set(df.file_index) == Set(0:count-1), "Incomplete or duplicate ROOT metadata report")
    require(all(>(0), df.last_event_number), "Generated counts must be positive")
    sort!(df, :file_index)
    return df.last_event_number
end

cpp_string(s) = "\"" * replace(s, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\r" => "\\r", "\t" => "\\t") * "\""

function read_root_counts(paths; root_executable=Sys.which("root"), executor=(cmd, _) -> run(cmd))
    require(root_executable !== nothing, "ROOT is not on PATH. On the cluster, load the ROOT module before extraction.")
    require(!isempty(paths) && allunique(paths), "Extraction needs unique ROOT paths")
    require(all(!any(ch -> ch in ('\n', '\r'), path) for path in paths), "ROOT paths must not contain newlines")
    return mktempdir() do directory
        list = joinpath(directory, "files.txt")
        report = joinpath(directory, "counts.tsv")
        write(list, join(paths, "\n") * "\n")
        macro_path = joinpath(@__DIR__, "read_generated_counts.C")
        invocation = "$macro_path($(cpp_string(list)),$(cpp_string(report)))"
        executor(Cmd([root_executable, "-l", "-b", "-q", invocation]), report)
        parse_counts(report, length(paths))
    end
end

function write_analysis_csv(output, table; overwrite=false)
    require(!isdir(output), "Analysis output must be a CSV file, not a directory")
    require(overwrite || !ispath(output), "Analysis CSV already exists; use --overwrite only after reviewing it")
    mkpath(dirname(output))
    mktemp(dirname(output)) do path, io
        CSV.write(io, table)
        close(io)
        mv(path, output; force=overwrite)
    end
end

function prepare_metadata(c, source; overwrite=false, reader=read_root_counts, root=PROJECT_ROOT)
    output = c["data"]["metadata_csv"]
    source = abspath(source)
    require(resolve_path(output) != source, "Refusing to overwrite the source CSV")
    require(!ispath(output) || realpath(output) != realpath(source), "Output aliases the source CSV")
    require(overwrite || !ispath(output), "Analysis CSV already exists; review it before using --overwrite")
    source_hash = bytes2hex(sha256(read(source)))
    data, paths = select_source(c, source; root)
    for path in paths
        require(isfile(path), "Missing ROOT input: $path")
    end
    before = stat.(paths)
    counts = reader(paths)
    require(length(counts) == length(paths) && all(x -> x isa Integer && x > 0, counts), "Invalid generated totals returned by ROOT")
    after = stat.(paths)
    require(all(a.size == b.size && a.mtime == b.mtime for (a, b) in zip(before, after)), "A ROOT file changed during extraction; rerun")
    require(bytes2hex(sha256(read(source))) == source_hash, "Source CSV changed during extraction; rerun")
    data.nb_events_original = copy(data.nb_events)
    data.nb_events = counts
    data.nb_events_source = fill("ROOT:last_event_number", size(data, 1))
    data.root_size_bytes = [s.size for s in after]
    data.root_mtime_unix = [s.mtime for s in after]
    data.extraction_utc = fill(string(now(UTC)), size(data, 1))
    data.source_csv = fill(source, size(data, 1))
    data.source_csv_sha256 = fill(source_hash, size(data, 1))
    write_analysis_csv(output, data; overwrite)
    return data
end

end
