# Time-dependent dead OM map, shared by the ND analysis and the data comparison.
#
# Data events are vetoed when an electron OM was dead at the event time
# (current_timestamp). Simulated events have no time, so each gets a
# pseudo-random time distributed like the data live time and the same veto
# applies. The time is derived from a hash of the event's energies, so every
# script vetoes the same simulated events.
#
# Map lines: [GID];[yyyy-mm-dd HH:MM:SS/yyyy-mm-dd HH:MM:SS);dead (times in UTC).
# The GID map lines are "<OM number> [GID]".
using CSV, Dates, UnROOT

const DEAD_OM_MAP_FILE = joinpath("data", "dead_om_map_all_phases.csv")
const OM_GID_MAP_FILE = joinpath("data", "snemo_om_num_gid_map.csv")

function read_om_gids(path)
    gids = Dict{String,Int}()
    for (i, line) in enumerate(eachline(path))
        isempty(strip(line)) && continue
        parts = split(strip(line))
        length(parts) == 2 || error("Cannot parse line $i of $path: $line")
        gids[parts[2]] = parse(Int, parts[1])
    end
    return gids
end

# OM number => sorted, merged dead intervals [start, end) in unix seconds.
function read_dead_om_map(map_path, gid_path)
    gids = read_om_gids(gid_path)
    format = dateformat"yyyy-mm-dd HH:MM:SS"
    intervals = Dict{Int,Vector{Tuple{Float64,Float64}}}()
    for (i, line) in enumerate(eachline(map_path))
        isempty(strip(line)) && continue
        fields = split(strip(line), ';')
        length(fields) == 3 || error("Cannot parse line $i of $map_path: $line")
        gid, period, state = fields
        state == "dead" || continue
        m = match(r"^\[(.+)/(.+)\)$", period)
        m === nothing && error("Cannot parse the period on line $i of $map_path: $period")
        haskey(gids, gid) || error("Unknown OM $gid on line $i of $map_path")
        lo = datetime2unix(DateTime(m[1], format))
        hi = datetime2unix(DateTime(m[2], format))
        lo < hi || error("Empty period on line $i of $map_path")
        push!(get!(() -> Tuple{Float64,Float64}[], intervals, gids[gid]), (lo, hi))
    end
    for (om, list) in intervals
        sort!(list)
        merged = [first(list)]
        for (lo, hi) in list[2:end]
            if lo <= last(merged)[2]
                merged[end] = (first(merged[end]), max(last(merged)[2], hi))
            else
                push!(merged, (lo, hi))
            end
        end
        intervals[om] = merged
    end
    return intervals
end

function is_dead(intervals, om, time)
    list = get(intervals, Int(om), nothing)
    list === nothing && return false
    i = searchsortedlast(list, (time, Inf))
    return i > 0 && time < list[i][2]
end

# Live periods of the data: each run from its first event for run_time seconds.
struct LiveTime
    starts::Vector{Float64}
    cumulative::Vector{Float64} # live time before each run
    total::Float64
end

function live_periods(event_runs, event_times, run_times::AbstractDict)
    first_time = Dict{Int,Float64}()
    for (run, time) in zip(event_runs, event_times)
        first_time[run] = min(get(first_time, run, Inf), time)
    end
    runs = sort!(collect(keys(first_time)))
    without = [r for r in runs if !haskey(run_times, r)]
    isempty(without) || error("No run time for data runs: $(join(without, ", "))")
    runs = [r for r in runs if run_times[r] > 0]
    isempty(runs) && error("The data has no live time")
    durations = [Float64(run_times[r]) for r in runs]
    return LiveTime([first_time[r] for r in runs], [0.0; cumsum(durations)[1:end-1]], sum(durations))
end

function read_run_times(path)
    table = CSV.File(path)
    return Dict(Int(r) => Float64(t) for (r, t) in zip(table.run_number, table.run_time)
                if !ismissing(r) && !ismissing(t))
end

function data_live_time(data_path, run_times_path)
    tree = LazyTree(ROOTFile(data_path), "Event", ["run_number", "current_timestamp"])
    return live_periods(Int.(collect(tree.run_number)), collect(tree.current_timestamp),
                     read_run_times(run_times_path))
end

# Time at a fraction u ∈ [0, 1) of the total live time.
function live_time_at(live, u)
    x = u * live.total
    i = searchsortedlast(live.cumulative, x)
    return live.starts[i] + (x - live.cumulative[i])
end

pseudo_uniform(key) = (hash(key) >>> 11) * 2.0^-53

simulated_time(live, energy_1, energy_2, energy_sum) =
    live_time_at(live, pseudo_uniform((Float64(energy_1), Float64(energy_2), Float64(energy_sum))))

# Veto for events in a table (simulation unless times are given).
function dead_om_veto(table, intervals, live; times = nothing)
    veto = falses(length(table.num_om_elec_f))
    for j in eachindex(veto)
        oms = table.num_om_elec_f[j]
        isempty(oms) && continue
        time = times === nothing ?
            simulated_time(live, table.energy_elec_1[j], table.energy_elec_2[j], table.energy_elec_sum[j]) :
            times[j]
        veto[j] = any(om -> is_dead(intervals, om, time), oms)
    end
    return veto
end
