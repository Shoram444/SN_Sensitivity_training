using Test, TOML, CSV, DataFramesMeta
include(joinpath(@__DIR__, "job_common.jl"))
using .SensitivityJobs
include(joinpath(@__DIR__, "analysis.jl"))
using .SensitivityAnalysis
using SNSensitivityEstimate, UnROOT

const EXAMPLE = joinpath(@__DIR__, "job.example.toml")

function small_config()
    c = load_config(EXAMPLE)
    c["outputs"]["background_groups"] = Dict("source" => ["background"])
    c["outputs"]["plots"] = String[]
    c["outputs"]["plot_formats"] = ["png"]
    c["optimization"]["time_limit_seconds"] = 0.01
    c["data"]["signal"] = "signal"
    for spec in values(cutset(c)["variables"])
        spec["role"] = "fixed"
    end
    cutset(c)["variables"]["diff_time_elec"]["bounds"] = [0.0, 2.5]
    return c
end

function synthetic_tree()
    n = 6
    return UnROOT.LazyTree((
        energy_elec_sum=fill(2.7, n), angle_3D_between_ep_em=fill(60.0, n),
        delta_y_elec=fill(-10.0, n), delta_z_elec=fill(10.0, n),
        diff_time_elec=[-1.0, 1.0, 1.5, 2.5, -2.5, 0.0],
        energy_elec_1=fill(1.3, n), energy_elec_2=fill(1.4, n),
        num_om_elec_f=[[1, 2] for _ in 1:n],
        closest_gamma=fill(100.0, n), closest_elec=fill(100.0, n),
        closest_track=fill(300.0, n), closest_time_track=fill(20.0, n),
        hit_the_same_calo_hit=zeros(Int, n), has_kinks=zeros(Int, n),
        vertex_3D_start_y=[[0.0, 0.0] for _ in 1:n],
        vertex_3D_start_z=[[0.0, 0.0] for _ in 1:n]))
end

function synthetic_metadata()
    return DataFrame(phase=fill("phase_3", 2), name=["signal", "background"],
        is_signal=[true, false], file_path=["signal.root", "background.root"],
        activity_bq=[1.0, 1e-8], amount=ones(2), nb_events=[10, 20],
        nb_events_source=fill("ROOT:last_event_number", 2))
end

@testset "Sensitivity workflow: no ROOT files loaded" begin
    @testset "Configuration and scheduler" begin
        c = load_config(EXAMPLE)
        @test c["data"]["exposure_seconds"] == 90192960.0
        @test cutset(c)["variables"]["diff_time_elec"]["bounds"] == [0.0, 5.0]
        for mutate in (
            x -> (x["slurm"]["time"] = "00:00:10"),
            x -> (x["optimization"]["seed"] = -1),
            x -> (cutset(x)["variables"]["diff_time_elec"]["bounds"] = [-5.0, 5.0]),
            x -> (cutset(x)["variables"]["energy_elec_sum"]["step"] = 0.2),
            x -> (x["outputs"]["background_groups"]["duplicate"] = ["bb"]),
            x -> (x["data"]["typo"] = true))
            changed = deepcopy(c); mutate(changed)
            @test_throws ArgumentError validate_config(changed)
        end
        grid = SensitivityAnalysis.optimization_grid(c)
        i = findfirst(==(:energy_elec_sum), grid.keys)
        @test grid.nsteps[i] == 11
        @test grid.mins[i]+grid.nsteps[i]*grid.steps[i] ≈ 3.5
        c["runtime"]["julia_module"] = ""
        c["runtime"]["julia_executable"] = "/bin/echo"
        script = render_sbatch(c, "/tmp/config with 'quote'.toml", "/tmp/run with spaces")
        @test success(pipeline(`bash -n`; stdin=IOBuffer(script)))
        output = read(pipeline(`bash`; stdin=IOBuffer(script)), String)
        @test occursin("/tmp/config with 'quote'.toml", output)
        @test occursin("--project=$(PROJECT_ROOT)", output)
    end

    @testset "Metadata provenance and uniqueness" begin
        mktempdir() do directory
            c = small_config()
            c["data"]["metadata_csv"] = joinpath(directory, "analysis.csv")
            data = synthetic_metadata()
            CSV.write(c["data"]["metadata_csv"], data)
            @test size(read_metadata(c; check_files=false), 1) == 2
            CSV.write(c["data"]["metadata_csv"], vcat(data, data[2:2, :]))
            @test_throws ArgumentError read_metadata(c; check_files=false)
            for mutate in (
                df -> (df.nb_events[2] = 0),
                df -> (df.nb_events_source[2] = "unverified"),
                df -> (df.activity_bq[2] = -1.0))
                changed = deepcopy(data); mutate(changed)
                CSV.write(c["data"]["metadata_csv"], changed)
                @test_throws ArgumentError read_metadata(c; check_files=false)
            end
            other = deepcopy(data)
            other.phase .= "phase_1"
            CSV.write(c["data"]["metadata_csv"], vcat(data, other))
            @test size(read_metadata(c; check_files=false), 1) == 2
            for name in ("signal.root", "background.root")
                write(joinpath(directory, name), "")
            end
            data.file_path = joinpath.(directory, data.file_path)
            CSV.write(c["data"]["metadata_csv"], data)
            config_path = joinpath(directory, "config.toml")
            write_toml(config_path, c)
            command = `$(Base.julia_cmd()) --project=$PROJECT_ROOT --startup-file=no $(joinpath(@__DIR__, "run.jl")) $config_path --check`
            @test occursin("Preflight passed: 2 processes; no ROOT events opened.", read(command, String))
        end
    end

    @testset "Absolute timing, boundary cuts, and fixed variables" begin
        c = small_config(); tree = synthetic_tree()
        selected, _, counts = SensitivityAnalysis.preselect(tree, c)
        @test length(selected) == 4
        @test counts[end] == 4
        @test all(diff(counts) .<= 0)
        @test SensitivityAnalysis.in_bounds(tree[1], "diff_time_elec", [0.0, 2.5])
        @test SensitivityAnalysis.in_bounds(tree[2], "diff_time_elec", [0.0, 2.5])
        @test !SensitivityAnalysis.in_bounds(tree[4], "diff_time_elec", [0.0, 2.5])
        @test !SensitivityAnalysis.in_bounds(tree[5], "diff_time_elec", [0.0, 2.5])
        @test SensitivityAnalysis.in_bounds(tree[6], "diff_time_elec", [0.0, 2.5])
        @test SensitivityAnalysis.in_bounds(tree[1], "delta_y_elec", [0.0, 200.0])
        params = SensitivityAnalysis.normalization(c)
        @test params !== SNparams
        @test params["tYear"] ≈ 2.86
        @test params["foilMass"]*params["a"]*params["tYear"] ≈ 17.4746
        c["data"]["exposure_seconds"] /= 2
        half = SensitivityAnalysis.normalization(c)
        @test half["tYear"] ≈ 1.43
        @test SNparams["tYear"] ≈ 2.86
    end

    @testset "Synthetic results, uncertainty, plots, and per-day validation" begin
        mktempdir() do directory
            c = small_config()
            c["validation"]["enabled"] = true
            c["validation"]["measured_data_files"] = ["measured.root"]
            c["validation"]["measured_exposure_seconds"] = 10*86400.0
            c["outputs"]["plots"] = ["background_spectrum"]
            loader = (_, _) -> (nothing, synthetic_tree())
            summary = analyze(c, synthetic_metadata(), directory; loader)
            expected = 4/20 * 1e-8 * c["data"]["exposure_seconds"]
            @test summary.signal_efficiency ≈ 4/10
            @test summary.background_counts ≈ expected
            @test summary.background_mc_stat_error ≈ sqrt(4)/20*1e-8*c["data"]["exposure_seconds"]
            @test summary.exposure_years ≈ 2.86
            params = SensitivityAnalysis.normalization(c)
            @test summary.sensitivity_years ≈ get_tHalf(params, 4/10, expected, c["optimization"]["alpha"]; approximate="table")
            @test summary.validation_status == "comparison_produced"
            components = CSV.read(joinpath(directory, "background_components.csv"), DataFrame)
            grouped = CSV.read(joinpath(directory, "background_groups.csv"), DataFrame)
            @test sum(components.expected_counts) ≈ sum(grouped.expected_counts) ≈ expected
            comparison = CSV.read(joinpath(directory, "data_vs_mc.csv"), DataFrame)
            @test sum(comparison.data_counts_per_day) ≈ 4/10
            @test sum(comparison.simulation_counts_per_day) ≈ expected/1043.9
            @test all(isfile(joinpath(directory, file)) for file in (
                "sensitivity.csv", "roi.csv", "cutflow.csv", "summary.txt",
                "physics_parameters.toml", "background_spectrum.png", "data_vs_mc.png"))
        end
    end

    @testset "Small optimizer with other variables fixed" begin
        c = small_config()
        time = cutset(c)["variables"]["diff_time_elec"]
        time["role"] = "optimize"
        context = SensitivityAnalysis.build_processes(c, synthetic_metadata(); loader=(_, _) -> (nothing, synthetic_tree()))
        roi, result = SensitivityAnalysis.optimize_roi(c, context.processes, SensitivityAnalysis.normalization(c))
        @test result.f_calls > 0
        @test roi.diff_time_elec[1] < roi.diff_time_elec[2]
        @test roi.energy_elec_sum == (2.4, 3.5)
        @test count(e -> SensitivityAnalysis.in_roi(e, roi), context.processes[1].data) > 0
        mktempdir() do directory
            SensitivityAnalysis.sensitivity_slices(c, context.processes,
                SensitivityAnalysis.normalization(c), roi, directory)
            slices = CSV.read(joinpath(directory, "sensitivity_slices.csv"), DataFrame)
            @test Set(slices.edge) == Set(["lower", "upper"])
            @test all(isfinite, slices.sensitivity_years)
            @test all(slices.background_counts .>= 0)
            @test isfile(joinpath(directory, "sensitivity_diff_time_elec.png"))
        end
    end
end
