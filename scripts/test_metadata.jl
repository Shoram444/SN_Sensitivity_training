using Test, CSV, DataFramesMeta, SHA
include(joinpath(@__DIR__, "job_common.jl"))
using .SensitivityJobs
include(joinpath(@__DIR__, "metadata.jl"))
using .SensitivityMetadata

function fixture(directory)
    c = load_config(joinpath(@__DIR__, "job.example.toml"))
    c["data"]["signal"] = "signal"
    c["data"]["metadata_csv"] = joinpath(directory, "analysis.csv")
    c["outputs"]["background_groups"] = Dict("source" => ["background"])
    data = DataFrame(phase=["phase_3", "phase_3", "phase_1", "phase_3", "phase_3"],
        name=["signal", "background", "background", "background", "other_signal"],
        is_signal=[true, false, false, false, true],
        file_path=["signal with 'quote'.root", "background.root", "missing_phase1.root",
            "reduces_files/background.root", "missing_other_signal.root"],
        activity_bq=[1.0, 1e-8, 1e-8, 1e-8, 1.0], amount=ones(5),
        nb_events=[10, 20, 20, 20, 10], phase_time_s=fill(6615636.0, 5))
    source = joinpath(directory, "source.csv")
    CSV.write(source, data)
    for path in data.file_path[1:2]
        write(joinpath(directory, path), "")
    end
    return c, source, data
end

@testset "Generated-count metadata" begin
    @testset "ROOT report shape and ordering" begin
        mktempdir() do directory
            report = joinpath(directory, "counts.tsv")
            CSV.write(report, DataFrame(file_index=[1, 0], last_event_number=[202, 101]); delim='\t')
            @test SensitivityMetadata.parse_counts(report, 2) == [101, 202]
            for table in (
                DataFrame(file_index=[0, 0], last_event_number=[101, 202]),
                DataFrame(file_index=[0], last_event_number=[101]),
                DataFrame(file_index=[0, 1], last_event_number=[101, 0]),
                DataFrame(file_index=[0, 1], last_event_number=[101, -1]),
                DataFrame(file_index=[0, 2], last_event_number=[101, 202]),
                DataFrame(file_index=[0, 1], wrong_column=[101, 202]))
                CSV.write(report, table; delim='\t')
                @test_throws ArgumentError SensitivityMetadata.parse_counts(report, 2)
            end
            calls = Ref(0)
            executor = function(cmd, output)
                calls[] += 1
                @test occursin("read_generated_counts.C(", last(cmd.exec))
                CSV.write(output, DataFrame(file_index=[1, 0], last_event_number=[202, 101]); delim='\t')
            end
            counts = read_root_counts(["a with spaces.root", "b'quote.root"];
                root_executable="mock-root-not-executed", executor)
            @test counts == [101, 202]
            @test calls[] == 1
            @test_throws ArgumentError read_root_counts(["duplicate", "duplicate"]; root_executable="unused")
            @test_throws ArgumentError read_root_counts(["a\nb"]; root_executable="unused")
        end
    end

    @testset "Separate CSV, provenance, refusal to overwrite" begin
        mktempdir() do directory
            c, source, original = fixture(directory)
            original_bytes = read(source)
            selected, paths = select_source(c, source; root=directory)
            @test selected.name == ["signal", "background"]
            @test length(paths) == 2
            output = prepare_metadata(c, source; reader=_ -> [101, 202], root=directory)
            @test output.nb_events == [101, 202]
            @test output.nb_events_original == [10, 20]
            @test output.activity_bq == original.activity_bq[1:2]
            @test output.phase_time_s == original.phase_time_s[1:2]
            @test all(==("ROOT:last_event_number"), output.nb_events_source)
            @test all(==(bytes2hex(sha256(original_bytes))), output.source_csv_sha256)
            @test read(source) == original_bytes
            @test size(read_metadata(c; root=directory), 1) == 2
            existing = read(c["data"]["metadata_csv"])
            @test_throws ArgumentError prepare_metadata(c, source; reader=_ -> error("must not run"), root=directory)
            @test read(c["data"]["metadata_csv"]) == existing
            @test_throws ErrorException prepare_metadata(c, source; overwrite=true,
                reader=_ -> error("injected ROOT failure"), root=directory)
            @test read(c["data"]["metadata_csv"]) == existing
            for reader in (_ -> [101], _ -> [101, 0], _ -> [101.0, 202.0])
                @test_throws ArgumentError prepare_metadata(c, source; overwrite=true, reader, root=directory)
                @test read(c["data"]["metadata_csv"]) == existing
            end
            write(paths[2], "changed")
            @test_throws ArgumentError read_metadata(c; root=directory)
            refreshed = prepare_metadata(c, source; overwrite=true, reader=_ -> [111, 222], root=directory)
            @test refreshed.nb_events == [111, 222]
            @test size(read_metadata(c; root=directory), 1) == 2
            alias = deepcopy(c); alias["data"]["metadata_csv"] = source
            @test_throws ArgumentError prepare_metadata(alias, source; overwrite=true, root=directory)
            @test read(source) == original_bytes
            CSV.write(source, vcat(original, original[2:2, :]))
            @test_throws ArgumentError select_source(c, source; root=directory)
        end
    end

    @testset "Changing inputs never publishes partial metadata" begin
        mktempdir() do directory
            c, source, _ = fixture(directory)
            @test_throws ArgumentError prepare_metadata(c, source; root=directory,
                reader=paths -> (write(paths[1], "changed"); [101, 202]))
            @test !ispath(c["data"]["metadata_csv"])
            @test_throws ArgumentError prepare_metadata(c, source; root=directory,
                reader=_ -> (open(io -> write(io, "\n"), source, "a"); [101, 202]))
            @test !ispath(c["data"]["metadata_csv"])
        end
    end

    @testset "ROOT integration on temporary parameter-only files" begin
        root = Sys.which("root")
        if root === nothing
            @info "ROOT integration not run: load the ROOT module on the cluster and rerun test_metadata.jl"
            @test_skip false
        else
            mktempdir() do directory
                macro_path = joinpath(directory, "create_parameters.C")
                write(macro_path, """
                #include <string>
                #include "TFile.h"
                #include "TParameter.h"
                void create_parameters(const char* dir) {
                    std::string base(dir);
                    { TFile f((base+"/int64.root").c_str(), "RECREATE");
                      TParameter<Long64_t> p("last_event_number", 5000000000LL); p.Write(); }
                    { TFile f((base+"/int.root").c_str(), "RECREATE");
                      TParameter<Int_t> p("last_event_number", 123); p.Write(); }
                    { TFile f((base+"/double.root").c_str(), "RECREATE");
                      TParameter<Double_t> p("last_event_number", 456.0); p.Write(); }
                    { TFile f((base+"/fractional.root").c_str(), "RECREATE");
                      TParameter<Double_t> p("last_event_number", 12.5); p.Write(); }
                    { TFile f((base+"/missing.root").c_str(), "RECREATE"); }
                }
                """)
                run(Cmd([root, "-l", "-b", "-q", "$macro_path($(SensitivityMetadata.cpp_string(directory)))"]))
                @test read_root_counts(joinpath.(directory, ["int64.root", "int.root", "double.root"])) == [5000000000, 123, 456]
                @test_throws ProcessFailedException read_root_counts([joinpath(directory, "fractional.root")])
                @test_throws ProcessFailedException read_root_counts([joinpath(directory, "missing.root")])
            end
        end
    end
end
