using Test, Dates

include("../analysis/dead_om_map.jl")

unix(s) = datetime2unix(DateTime(s, dateformat"yyyy-mm-dd HH:MM:SS"))

@testset "Dead OM map" begin
    mktempdir() do dir
        gid_path = joinpath(dir, "gids.csv")
        write(gid_path, "0 [1301:0.0.0.0]\n1 [1301:0.0.0.1]\n520 [1231:0.0.0.0.0]\n")
        map_path = joinpath(dir, "map.csv")
        write(map_path, """
            [1301:0.0.0.0];[2025-01-01 00:00:00/2025-01-02 00:00:00);dead
            [1301:0.0.0.0];[2025-01-02 00:00:00/2025-01-03 00:00:00);dead
            [1301:0.0.0.0];[2025-01-05 00:00:00/2025-01-06 00:00:00);dead
            [1231:0.0.0.0.0];[2025-01-01 00:00:00/2025-01-02 00:00:00);dead
            """)
        intervals = read_dead_om_map(map_path, gid_path)
        # Adjacent periods merge; [start, end) boundaries.
        @test intervals[0] == [(unix("2025-01-01 00:00:00"), unix("2025-01-03 00:00:00")),
                               (unix("2025-01-05 00:00:00"), unix("2025-01-06 00:00:00"))]
        @test is_dead(intervals, 0, unix("2025-01-01 00:00:00"))
        @test is_dead(intervals, 0, unix("2025-01-02 12:00:00"))
        @test !is_dead(intervals, 0, unix("2025-01-03 00:00:00"))
        @test !is_dead(intervals, 0, unix("2025-01-04 00:00:00"))
        @test !is_dead(intervals, 1, unix("2025-01-01 12:00:00"))
        @test is_dead(intervals, Int32(520), unix("2025-01-01 12:00:00"))

        write(map_path, "[1301:9.9.9.9];[2025-01-01 00:00:00/2025-01-02 00:00:00);dead\n")
        @test_throws ErrorException read_dead_om_map(map_path, gid_path)
    end

    # Two runs of 100 s and 300 s; run 3 has no live time.
    live = live_periods([1, 1, 2, 3], [1000.0, 1050.0, 5000.0, 9000.0], Dict(1 => 100.0, 2 => 300.0, 3 => 0.0))
    @test live.total == 400.0
    @test live_time_at(live, 0.0) == 1000.0
    @test live_time_at(live, 0.2) == 1080.0
    @test live_time_at(live, 0.5) == 5100.0
    @test_throws ErrorException live_periods([1, 4], [1.0, 2.0], Dict(1 => 10.0))

    # Simulated times are deterministic and follow the live time.
    @test simulated_time(live, 0.5, 0.7, 1.2) == simulated_time(live, 0.5f0 + 0.0, 0.7, 1.2)
    times = [simulated_time(live, rand(), rand(), rand()) for _ in 1:20_000]
    @test all(t -> 1000 <= t < 1100 || 5000 <= t < 5300, times)
    @test isapprox(count(<(2000), times) / length(times), 0.25; atol = 0.02)

    table = (num_om_elec_f = [[0, 1], [1, 2], [0, 1]], energy_elec_1 = [0.5, 0.5, 0.5],
             energy_elec_2 = [0.7, 0.7, 0.7], energy_elec_sum = [1.2, 1.2, 1.2])
    intervals = Dict(0 => [(1000.0, 1100.0)])
    @test dead_om_veto(table, intervals, live; times = [1050.0, 1050.0, 5000.0]) == [true, false, false]
    # Every simulated time falls in run 1 when only run 1 has live time.
    only_run_1 = live_periods([1], [1000.0], Dict(1 => 100.0))
    @test dead_om_veto(table, intervals, only_run_1) == [true, false, true]
end
