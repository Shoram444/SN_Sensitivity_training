using Test, CSV, DataFramesMeta

const SNparams = Dict("t" => 200.0, "tYear" => 2.0, "foilMass" => 6.25)
include("../analysis/nd_background_report.jl")

@testset "Log spectrum steps close at zero bins and outer edges" begin
    x, y = spectrum_steps(0:3, [2.0, 0.0, 5.0], 0.01)
    @test collect(zip(x, y)) == [(0, 0.01), (0, 2.0), (1, 0.01), (2, 5.0), (3, 0.01)]
end

@testset "Background report" begin
    mktempdir() do dir
        path = joinpath(dir, "results.csv")
        energy_prefilter = (1.0, 3.5)
        best_roi = (
            energy_elec_sum = (2.0, 3.0),
            diff_time_elec = (0.0, 2.0),
            __nd_cut_closest_gamma = (50.0, Inf),
        )
        files = DataFrame(
            name = ["bb", "bi_pmt", "tl_pmt"], amount = [2.0, 1.0, 1.0],
            activity_bq = [1.0, 2.0, 3.0], nb_events = [100, 50, 100],
        )
        event(e; time = -1.0, gamma = -50.0) =
            (energy_elec_sum = e, diff_time_elec = time, __nd_cut_closest_gamma = abs(gamma))
        tables = [
            [event(0.5), event(2.0), event(3.0), event(3.75), event(4.0),
             event(2.5; time = -2.0), event(2.5; gamma = 49.0)],
            [event(2.5)],
            [event(2.5), event(3.75)],
        ]
        signal_process = (isotopeName = "bb0nuM1", timeMeas = 100.0)
        save_background_report(path, "phase_4", signal_process, files, tables, best_roi,
                               energy_prefilter, 0.0:0.05:4.0, 9.0)
        @test !isfile(path) # The report uses the in-memory result, not a results CSV.
        report = CSV.read(joinpath(dir, "background_components.csv"), DataFrame)
        @test report.process == ["bb", "bi_pmt", "tl_pmt", "Total"]
        @test report.mc_events_roi == [1, 1, 1, 3]
        @test report.mc_events_0_4_mev == [4, 1, 2, 7]
        @test report.expected_roi == [2.0, 4.0, 3.0, 9.0]
        @test report.expected_0_4_mev == [8.0, 4.0, 6.0, 18.0]
        @test report.mc_stat_error_roi ≈ [2.0, 4.0, 3.0, sqrt(29.0)]
        @test report.mc_stat_error_0_4_mev ≈ [4.0, 4.0, sqrt(18.0), sqrt(50.0)]
        @test open(io -> read(io, 8), joinpath(dir, "background_spectrum.png")) ==
              UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
        @test open(io -> String(read(io, 5)), joinpath(dir, "background_spectrum.pdf")) == "%PDF-"
    end
end
