#### 
#### WORK IN PROGRESS
####
import Pkg;
Pkg.activate("/sps/nemo/scratch/mpetro/Sensitivity_training/")
using SNSensitivityEstimate, UnROOT, FHist, CairoMakie, CSV, DataFramesMeta
using Random, LinearAlgebra, Statistics, Distributions, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, SpecialFunctions, ValueShapes
using Distributions, ColorSchemes


# STEP 1: Load the data

# Start is the same as in ex1_background.jl, we load data and create DataProcess objects. 
data_info = CSV.read("data/mathis_data_info_final.csv", DataFrame)

# Select which phase we process
phase = "phase_3"

selected_files = @chain data_info begin
    @subset :phase .== phase
    # @subset .!occursin.("bb0nu", :file_path)
    @subset .!occursin.("reduces", :file_path) 
end
# manually add nu0bb path because that's not in reduced files
append!(selected_files, data_info[data_info.file_path .== "data/mathis_files/bb0nu/phase_0/simu_0/rec_bb0nu_0_0.root", :])


data_files = [UnROOT.ROOTFile(file_path) for file_path in selected_files.file_path] 
tree_name = "Event" 

# we don't need all the variables, we can select a few to use for the analysis.
var_names = [ 
    "energy_elec_sum", 
    "angle_3D_between_ep_em", 
    "delta_y_elec",
    "delta_z_elec",
    "diff_time_elec",
    "energy_elec_1",
    "energy_elec_2",
    "hit_the_same_calo_hit",
    "closest_gamma",
    "closest_elec",
    "closest_track",
    "closest_time_track",
    "has_kinks",
    "num_om_elec_f",
    "vertex_3D_start_y",
    "vertex_3D_start_z"
    ]

data_tables = [UnROOT.LazyTree(data_file, tree_name, var_names) for data_file in data_files] 


# define the boundaries of the data-cuts, we use simple Esum >300keV, dy < 100mm, dz < 100mm, Pint > 0.01, Pext < 0.01
var_bounds = (
    energy_elec_sum = (0.7, 4),
    angle_3D_between_ep_em = (0, 180),
    delta_y_elec = (0, 100),
    delta_z_elec = (0, 100),
    diff_time_elec = (-1.5, 1.5),
    energy_elec_1 = (0.35, 4),
    energy_elec_2 = (0.35, 4),
)

# Apply bounds only for keys present in `var_bounds`.
filtered_data_tables = [SNSensitivityEstimate.filter_data(table, var_names, var_bounds) for table in data_tables]


# Apply additional topological cuts and store filtered tables back.
for i in eachindex(filtered_data_tables)
    table = filtered_data_tables[i]
    cut_vertex_dist = sqrt.(table.delta_y_elec .* table.delta_y_elec .+ table.delta_z_elec .* table.delta_z_elec) .>= 100
    cut_same_calo   = table.hit_the_same_calo_hit .== 1
    cut_gam         = abs.(table.closest_gamma) .<= 50
    cut_elec        = abs.(table.closest_elec) .<= 25
    cut_track       = abs.(table.closest_track) .< 200
    cut_ttrack      = abs.(table.closest_time_track) .< 10
    cut_kinks       = table.has_kinks .== 1

    n_evt = length(table)
    cut_om_size = [length(table.num_om_elec_f[j]) != 2 for j in 1:n_evt]

    cut_om_range = falses(n_evt)
    for j in 1:n_evt
        if !cut_om_size[j]
            om0 = table.num_om_elec_f[j][1]
            om1 = table.num_om_elec_f[j][2]
            cut_om_range[j] = (om0 > 520 || om1 > 520)
        end
    end

    cut_vertex_region = falses(n_evt)
    for j in 1:n_evt
        if !cut_om_size[j] && !cut_om_range[j] && length(table.vertex_3D_start_y[j]) >= 2 && length(table.vertex_3D_start_z[j]) >= 2
            y_mid = 0.5 * (table.vertex_3D_start_y[j][1] + table.vertex_3D_start_y[j][2])
            z_mid = 0.5 * (table.vertex_3D_start_z[j][1] + table.vertex_3D_start_z[j][2])

            cut_vertex_region[j] = (
                ((y_mid - 1180) / 60)^2 + ((z_mid + 850) / 140)^2 <= 1 ||
                ((y_mid - 1144) / 25)^2 + ((z_mid - 196) / 40)^2 <= 1 ||
                ((y_mid - 1114) / 25)^2 + ((z_mid + 257) / 40)^2 <= 1 ||
                abs(y_mid) > 2360 || abs(z_mid) > 1340
            )
        end
    end

    # `cut_om_good` intentionally omitted here: OM quality map/helper is defined externally.
    keep = .!(cut_vertex_dist .| cut_same_calo .| cut_gam .| cut_elec .| cut_track .| cut_ttrack .| cut_kinks .| cut_om_size .| cut_om_range .| cut_vertex_region)
    filtered_data_tables[i] = table[keep]
end


processes = [
    SNSensitivityEstimate.DataProcess(
        collect(filtered_data_tables[i].energy_elec_sum), # Here we specify a single variable to be used. In later stages we will use N-dimensional data vectors.
        String(selected_files.name[i]), 
        selected_files.is_signal[i], 
        selected_files.activity_bq[i], 
        selected_files.time_s[i], 
        selected_files.nb_events[i], 
        0:0.1:4, 
        selected_files.amount[i], 
    ) for i in 1:length(filtered_data_tables)
]



mutable struct HistogramData
    hist::FHist.Hist1D
    name::String
    is_signal::Bool
end

function HistogramData(process::SNSensitivityEstimate.DataProcess)
    hist = get_bkg_counts_1D(process)
    HistogramData(hist, process.isotopeName, process.signal)
end

all_bkg_histograms = [HistogramData(process) for process in processes if !process.signal]

function merge_hist_by_keyword(histograms::Vector{HistogramData}, keyword::String, merge_name::String)
    merged_hist = FHist.Hist1D(;binedges = histograms[1].hist.binedges)
    is_signal = findfirst(hist -> occursin(keyword, hist.name) && hist.is_signal, histograms) !== nothing
    for hist in histograms
        if occursin(keyword, hist.name)
            # if signal changes within keyword error out
            if hist.is_signal != is_signal
                error("Can only merge histograms with the same signal flag for keyword $keyword")
            end
            merge!(merged_hist, hist.hist)
        end

    end
    return HistogramData(merged_hist, merge_name, is_signal)
end

pmt_hists = merge_hist_by_keyword(all_bkg_histograms, "pmt", "pmt")
other_histos = filter(hist -> !occursin("pmt", hist.name), all_bkg_histograms)
all_histos = vcat(other_histos, [pmt_hists])

let
    colors = cgrad(:managua100, length(all_histos); categorical=true)
    histos = [hist.hist for hist in all_histos]

    f = CairoMakie.Figure(size = (1000, 600))
    a = CairoMakie.Axis(
        f[1, 1], 
        title = "Stacked Histogram of Expected Background Counts, phase $phase data\nExposure of 17.5 kg.yr", 
        xlabel = "sum of energies", 
        ylabel = "expected counts", 
        # ylabel = "log of expected counts", 
        # yscale = log10
        )
    p = CairoMakie.hist!(a, sum(histos), color = colors[1])
    for i in length(histos):-1:1
        CairoMakie.hist!(a, sum(histos[1:i]), color = colors[i])
    end

    labels = [String(p.name) for p in all_histos]
    elements = [CairoMakie.PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    title = "Processes"
    
    ylims!(a, 1e-3, nothing)
    CairoMakie.Legend(f[1,2], elements, labels, title)
    # save("data/out/ex1/stacked_background_histogram.png", f)
    f
end

let
    adjusted_per_day_histos = [hist.hist .* (1/2.86/365) for hist in all_histos]
    colors = cgrad(:jet1, length(adjusted_per_day_histos); categorical=true)
    f = CairoMakie.Figure(size = (1000, 600))
    a = CairoMakie.Axis(
        f[1, 1], 
        title = "Histogram of Expected Background Counts, phase $phase data", 
        xlabel = "sum of energies", 
        # ylabel = "expected counts per day", 
        ylabel = "log of expected counts", 
        yscale = log10
        )
    p = CairoMakie.stephist!(a, adjusted_per_day_histos[1], color = colors[1], linewidth = 5)
    for i in 2:length(adjusted_per_day_histos)
        CairoMakie.stephist!(a, adjusted_per_day_histos[i], color = colors[i], linewidth = 5)
    end
    CairoMakie.stephist!(a, sum(adjusted_per_day_histos), color = :black, linewidth = 5)



    labels = [hist.name for hist in all_histos]
    push!(labels, "total")
    elements = [CairoMakie.PolyElement(polycolor = colors[i]) for i in 1:length(labels)-1]
    push!(elements, CairoMakie.PolyElement(polycolor = :black))
    title = "Processes"
    
    ylims!(a, 1e-5, nothing)
    CairoMakie.Legend(f[1,2], elements, labels, title)
    save("data/out/ex5/stacked_background_histogram_per_day_phase_$phase.png", f, px_per_unit = 2)
    f
end



for i=1:length(processes)
    println("i: $i, process: ", processes[i].isotopeName, ", n = ", integral(get_bkg_counts_1D(processes[i])), ",  activity = ", processes[i].activity, " Bq", " n_pass = ", length(processes[i].dataVector))
end

println("Total expected counts: ", integral(sum(get_bkg_counts_1D.(processes))))
println("\n\n\n")

for i=1:length(processes)
    n_pass = length(processes[i].dataVector)
    eff = n_pass/processes[i].nTotalSim
    a = processes[i].activity
    t = 24*3600
    n_exp_day = round(eff*a*t, digits = 2)
    println("i: $i, process: ", processes[i].isotopeName, ", n_exp_day = ", n_exp_day, ", n_pass = ", n_pass)
end
println("Total expected counts per day: ", round(integral(sum(get_bkg_counts_1D.(processes)))*(1/2.86/365), digits = 2))


# Sensitivity
# Pick signal
signal = get_process("bb0nu", processes)[1]
background = [p for p in processes if !p.signal]

α = 1.65

t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)
best_sens = get_tHalf(SNparams, effbb, expBkgESum, α; approximate="table")
ThalfbbESum = round(best_sens, sigdigits=3)
