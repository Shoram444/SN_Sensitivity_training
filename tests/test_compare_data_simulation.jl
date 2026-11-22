using Test, CSV, DataFramesMeta, UnROOT

include("../scripts/compare_data_simulation.jl")

@testset "ROI parsing" begin
    roi = parse_roi("(energy_elec_sum = (2.25, 2.9000000000000004), diff_time_elec = (0.0, 1.5), " *
                    "__nd_cut_closest_gamma = (50.0, Inf))")
    @test roi == (energy_elec_sum = (2.25, 2.9000000000000004), diff_time_elec = (0.0, 1.5),
                  __nd_cut_closest_gamma = (50.0, Inf))
    @test parse_roi("(energy_elec_sum = (-1.0, 2.0))") == (energy_elec_sum = (-1.0, 2.0),)
    @test_throws ErrorException parse_roi("(energy_elec_sum = (2.0, run(`ls`)))")
end

@testset "Event selection" begin
    settings = (
        extra_var_names = String[],
        use_same_calo_cut = true, use_kinks_cut = true,
        use_om_cut = true, use_vertex_region_cut = true,
        nd_numeric_cuts = [
            (name = :closest_gamma, value = e -> abs(e.closest_gamma), comparison = (>),
             mode = :optimize, fixed = 50.0, bounds = (20.0, 150.0), step = 10.0),
            (name = :closest_track, value = e -> abs(e.closest_track), comparison = (>),
             mode = :fixed, fixed = 200.0, bounds = (100.0, 400.0), step = 10.0),
        ],
    )
    roi = (energy_elec_sum = (2.0, 3.0), diff_time_elec = (0.0, 1.5),
           __nd_cut_closest_gamma = (50.0, Inf))
    sel = build_selection(settings, roi, Set([7]))
    @test sel.windows == [:diff_time_elec => (0.0, 1.5)]
    @test only(sel.fixed).name == :closest_track

    good = (energy_elec_sum = 1.2, diff_time_elec = -1.0, closest_gamma = 60.0, closest_track = 250.0,
            hit_the_same_calo_hit = false, has_kinks = false, num_om_elec_f = Int32[1, 2],
            vertex_3D_start_y = [0.0, 0.0], vertex_3D_start_z = [0.0, 0.0])
    events = [
        good,                                               # passes; outside the energy ROI
        merge(good, (energy_elec_sum = 2.5,)),              # passes
        merge(good, (diff_time_elec = 1.5,)),               # ROI window is [lo, hi)
        merge(good, (closest_gamma = 50.0, energy_elec_sum = 2.6)), # optimized threshold is >=
        merge(good, (closest_track = 200.0,)),              # fixed cut keeps its comparison (>)
        merge(good, (has_kinks = true,)),
        merge(good, (hit_the_same_calo_hit = true,)),
        merge(good, (num_om_elec_f = Int32[1, 7],)),        # bad OM
        merge(good, (num_om_elec_f = Int32[1, 600],)),      # outside the main walls
        merge(good, (num_om_elec_f = Int32[1],)),
        merge(good, (vertex_3D_start_y = [0.0, 0.0], vertex_3D_start_z = [1400.0, 1400.0])),
    ]
    columns = NamedTuple(k => [getproperty(e, k) for e in events] for k in keys(good))
    energies = select_energies(LazyTree(columns), sel)
    @test energies == [1.2, 2.5, 2.6]
    @test bin_counts(energies, 0.0:0.5:4.0) == [0, 0, 1, 0, 0, 2, 0, 0]
    @test in_window(energies, roi.energy_elec_sum) == 2

    @test_throws ErrorException build_selection(settings, (energy_elec_sum = (2.0, 3.0),), Set{Int}())
end

@testset "Simulation paths" begin
    @test simulation_path("data/mathis_files/bb/phase_3/bb.root", "/local/sim") == "/local/sim/bb/phase_3/bb.root"
    @test simulation_path("data/majoron_reduced/x.root", "/local/sim") == joinpath(PROJECT_DIR, "data/majoron_reduced/x.root")
    @test simulation_path("data/mathis_files/bb.root", nothing) == joinpath(PROJECT_DIR, "data/mathis_files/bb.root")
end

@testset "Live time" begin
    mktempdir() do dir
        path = joinpath(dir, "run_times.csv")
        write(path, "file_path,run_number,run_time\n\"a.root\",10,100.0\n\"b.root\",11,50.0\n\"c.root\",12,25.0\n")
        @test live_time(path, [11, 10], :data) == (150.0, [10, 11])
        @test (@test_logs (:warn,) live_time(path, [10], 10:13)) == (175.0, [10, 11, 12])
        @test_throws ErrorException live_time(path, [10, 13], :data)

        write(path, "file_path,run_number,run_time\n\"a.root\",10,100.0\n\"merged.root\",10,175.0\n")
        @test_throws ErrorException live_time(path, [10], :data)
        write(path, "file_path,run_number,run_time\n\"a.root\",10,\n")
        @test_throws ErrorException live_time(path, [10], :data)
    end
end

@testset "Chi2" begin
    data = [4, 9, 0, 3, 5]
    model = [4.0, 4.0, 1.0, 0.0, 2.0]
    variance = [0.0, 1.0, 0.0, 0.0, 1.0]
    test = chi2_test(data, model, variance, [true, true, true, true, false])
    @test test.ndf == 3
    @test test.unpredicted_bins == 1
    @test test.chi2 ≈ 0 + 25 / 5 + 1
    @test test.chi2_poisson ≈ 2 * ((4 - 9 + 9 * log(9 / 4)) + 1)
    @test residuals(data, model, variance)[1:3] ≈ [0.0, 5 / sqrt(5), -1.0]
    @test isnan(residuals(data, model, variance)[4])
end

@testset "Comparison plot" begin
    mktempdir() do dir
        path = joinpath(dir, "plot.png")
        edges = 0.0:1.0:4.0
        save_comparison_plot(path, edges, Dict("bb" => [2.0, 1.0, 0.0, 0.0], "bi_pmt" => [0.5, 0.0, 0.1, 0.0]),
                             Dict("bb" => [0.1, 0.1, 0.0, 0.0], "bi_pmt" => [0.01, 0.0, 0.01, 0.0]),
                             [3.0, 1.0, 0.0, 0.5], [1.0, 1.0, 0.0, 0.5], [0.5, 0.0, NaN, NaN],
                             (2.0, 3.0), (0.0, 3.0), "Title", "Subtitle")
        @test open(io -> read(io, 4), path) == UInt8[0x89, 0x50, 0x4e, 0x47]
        @test open(io -> String(read(io, 5)), joinpath(dir, "plot.pdf")) == "%PDF-"
    end
end
