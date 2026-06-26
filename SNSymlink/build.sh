module load julia
export JULIA_DEPOT_PATH="$HOME/sps_mpetro/PROGRAMS/.julia_sym"
export JULIA_CPU_TARGET="generic"

rm -rf "$JULIA_DEPOT_PATH/compiled/v1.10"
rm -rf "$JULIA_DEPOT_PATH/pkgimages/v1.10"


# use JULIA_DEPOT_PATH and JULIA_CPU_TARGET and build a sysimage for the SNSensitivityEstimate package and its dependencies
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile(); using PackageCompiler; create_sysimage([:SNSensitivityEstimate, :UnROOT, :FHist, :CSV, :DataFramesMeta, :BAT, :BinnedModels, :StatsBase, :DensityInterface, :IntervalSets, :SpecialFunctions, :ValueShapes]; sysimage_path=\"SNSymlink/SNPackageSysimage.so\", cpu_target=ENV[\"JULIA_CPU_TARGET\"])"

# test
julia --project=. -J "SNSymlink/SNPackageSysimage.so" -e 'using SNSensitivityEstimate; println("SNSensitivityEstimate loaded successfully with sysimage!")'
