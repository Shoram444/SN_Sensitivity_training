using Test

include("../scripts/sensitivity_vs_time.jl")

@testset "Sensitivity versus time" begin
    mktempdir() do dir
        write(joinpath(dir, "settings.jl"), "global signal_name = \"bb0nuM1\"\nglobal α = 1.64485362695147\n")
        run = (signal = "bb0nuM1", α = 1.64485362695147, efficiency = 0.02, background = 400.0)
        nominal = sensitivity(run, SNparams["tYear"])
        write(joinpath(dir, "results.csv"), "phase,best_sensitivity,best_ROI,best_bkg,best_eff\nphase_4,$nominal,\"(a = (0, 1))\",400.0,0.02\n")
        parsed = read_run(joinpath(dir, "results.csv"))
        @test parsed.signal == "bb0nuM1" && parsed.α == run.α && parsed.background == 400.0

        table = sensitivity_table(parsed, [0.5, SNparams["tYear"]])
        @test table.background_counts ≈ [400.0 * 0.5 / SNparams["tYear"], 400.0]
        @test table.exposure_kg_yr ≈ SNparams["foilMass"] * SNparams["a"] .* [0.5, SNparams["tYear"]]
        @test table.exposure_kg_yr[2] ≈ 17.475 atol = 0.01
        @test table.sensitivity_yr[2] ≈ nominal
        @test table.sensitivity_yr[1] < nominal
    end
end
