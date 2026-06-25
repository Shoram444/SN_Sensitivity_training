##############################################################################
##############################################################################
##############################################################################
##
## DESCRIPTION: This example shows how to calculate the frequentist (Feldman-Cousins) sensitivity for a given signal process
## using the SNSensitivityEstimate package. 
##############################################################################
##############################################################################
##############################################################################


import Pkg;
Pkg.activate("/sps/nemo/scratch/mpetro/Sensitivity_training/")
using SNSensitivityEstimate, UnROOT, FHist, CairoMakie, CSV, DataFramesMeta
using Distributions 

# STEP 1: Load the data

# Start is the same as in ex1_background.jl, we load data and create DataProcess objects. 
data_info = CSV.read("data/data_info.csv", DataFrame)
file_paths = joinpath.("data", data_info.file_name) 
data_files = [UnROOT.ROOTFile(file_path) for file_path in file_paths] 
tree_name = "tree" 
variables = keys(data_files[1][tree_name]) 
data_tables = [UnROOT.LazyTree(data_file, tree_name, variables) for data_file in data_files] 


# STEP 1.1: apply data-cuts to the data before creating DataProcess objects.
# In ex1, we did not impose any data-cuts on the data.
# If we want somewhat meaningful results on sensitivity, we must first filter data to some meaningful data-cuts:

# First we define variables we want to cut
var_names = [ 
    "sumE", 
    "dy", 
    "dz",
    "Pint",
    "Pext",
    ]

# Next we define the boundaries of the data-cuts, we use simple Esum >300keV, dy < 100mm, dz < 100mm, Pint > 0.01, Pext < 0.01
var_bounds = (
    sumE = (300, 3500),
    dy = (0, 100),
    dz = (0, 100),
    Pint = (0.01, 1),
    Pext = (0, 0.01),
)

# next we use filter_data! function from SNSensitivityEstimate to filter the data before loading it to DataProcess.
filtered_data_tables = [SNSensitivityEstimate.filter_data(table, var_names, var_bounds) for table in data_tables]


processes = [
    SNSensitivityEstimate.DataProcess(
        collect(filtered_data_tables[i].sumE), # Here we specify a single variable to be used. In later stages we will use N-dimensional data vectors.
        String(data_info.process[i]), 
        data_info.is_signal[i], 
        data_info.activity[i], 
        data_info.time_s[i], 
        data_info.n_sim[i], 
        0:100:3500, 
        data_info.amount[i], 
    ) for i in 1:length(filtered_data_tables)
]

# select signal and background processes

signal_name = "bb0nu_foil_bulk"
signal = SNSensitivityEstimate.get_process(signal_name, processes)[1]

background_names = data_info.process[data_info.process .!= signal_name]
background = [SNSensitivityEstimate.get_process(String(name), processes)[1] for name in background_names]

# STEP 2: The general idea of the 1D frequentist approach is to first find the most optimal ROI in 1 selected variable (in this case sumE) for a given signal process.
# Then for the optimal ROI, we calculate the signal efficiency and expected background counts (run it through FC tables) and then calculate the sensitivity.
# As shown in #6090-v1 the way we do this is split the ESum range into bins and then for each bin, we calculate the sensitivity and then find the bin with the best sensitivity.

# (Optional) STEP 2.1: we can first show a very simple approach - choose a single ROI and calculate the sensitivity for that ROI. This is not optimal, but it is a good starting point to understand the process.
my_ROI = (2500, 3200) # some arbitrary ROI, we will later find the optimal ROI.

my_eff = SNSensitivityEstimate.lookup(signal, my_ROI[1], my_ROI[2]) 
my_expBkg = SNSensitivityEstimate.get_bkg_counts_ROI(my_ROI[1], my_ROI[2], background...)
println("For ROI = $my_ROI, the signal efficiency is $my_eff and the expected background counts is $my_expBkg")

# Now easily, for the sensitivity @90% simply:
α = 1.64485362695147
my_sens = SNSensitivityEstimate.get_tHalf(SNparams, my_eff, my_expBkg, α; approximate="table")
println("For ROI = $my_ROI, the sensitivity is $my_sens")

# We can visualize what this looks like in the background model plot with zoom to ROI

let
    # calculate signal activity for the given sensitivity
    my_0nu_acitvity=SNSensitivityEstimate.halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], my_sens*365*24*3600)
    SNSensitivityEstimate.set_activity!(signal, my_0nu_acitvity) 
    
    bkg_processes = [p for p in processes if !p.signal]
    bkg_histograms = [SNSensitivityEstimate.get_bkg_counts_1D(p) for p in bkg_processes]
    signal_histogram = SNSensitivityEstimate.get_bkg_counts_1D(signal)
    colors = ["#041E42", "#BE4D00", "#951272", "#006630", "#005C8A", "#FFB948", "#605643", "#302D23"] # get our supernemo colors
    f = CairoMakie.Figure(size = (1800, 800), fontsize = 24)
    a1 = CairoMakie.Axis(f[1, 1], title = "Background model with signal overlay", xlabel = "sum of energies", ylabel = "expected counts", yscale = log10)
    a2 = CairoMakie.Axis(f[1, 2], title = "Zoom to ROI", xlabel = "sum of energies", ylabel = "expected counts")

    for axis in (a1, a2)
        for i in length(bkg_histograms):-1:1
            CairoMakie.stephist!(axis, sum(bkg_histograms[1:i]), color = colors[i], linewidth = 5)
        end
        CairoMakie.stephist!(axis, signal_histogram, color = :black, linewidth = 4, linestyle = :dash)
        CairoMakie.vspan!(axis, my_ROI[1], my_ROI[2], color = (:black, 0.08))
    end

    labels = [String(p.isotopeName) for p in bkg_processes]
    push!(labels, "signal")

    elements = [CairoMakie.LineElement(color = colors[i], linewidth = 5) for i in 1:length(bkg_processes)]
    push!(elements, CairoMakie.LineElement(color = :black, linestyle = :dash, linewidth = 5))
    title = "Processes"
    
    ylims!(a1, 1e-3, 1e5)
    ylims!(a2, 0, 5)
    xlims!(a2, my_ROI[1], my_ROI[2])

    CairoMakie.Legend(f[1,3], elements, labels, title, patchsize = (45, 20))
    save("data/out/ex3/my_roi_plot.png", f, px_per_unit = 2)
    f

end

# STEP 3: Now we can find the optimal ROI by scanning through all possible ROIs and calculating the sensitivity for each ROI. We can then find the ROI with the best sensitivity.
# We do so by creating the so called t12Map, which is a map of sensitivities for each ROI. The ROI with the best sensitivity is the one with the lowest t12 value.

t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)
best_sens = get_tHalf(SNparams, effbb, expBkgESum, α; approximate="table")
ThalfbbESum = round(best_sens, sigdigits=3)

with_theme(theme_latexfonts()) do
    f = CairoMakie.Figure(size = (1500, 950), fontsize = 34, figure_padding = 18)
    a = CairoMakie.Axis(
        f[1, 1],
        xlabel = "ROI min (keV)",
        ylabel = "ROI max (keV)",
        title = "1D frequentist sensitivity map in E_sum",
    )

    p = CairoMakie.plot!(a, t12MapESum, colormap = :coolwarm)
    CairoMakie.scatter!(a, [best_t12ESum[:minBinEdge]], [best_t12ESum[:maxBinEdge]], color = :darkgreen, marker = :xcross, markersize = 42, strokewidth = 5)
    CairoMakie.vlines!(a, [best_t12ESum[:minBinEdge]], color = (:darkgreen, 0.9), linestyle = :dash, linewidth = 4)
    CairoMakie.hlines!(a, [best_t12ESum[:maxBinEdge]], color = (:darkgreen, 0.9), linestyle = :dash, linewidth = 4)
    CairoMakie.textlabel!(
        a,
        best_t12ESum[:minBinEdge],
        best_t12ESum[:maxBinEdge],
        text = "best ROI = ($(best_t12ESum[:minBinEdge]), $(best_t12ESum[:maxBinEdge])) keV\nT_1/2 = $(ThalfbbESum) yr",
        text_align = (:left, :bottom),
        offset = (-450, -80),
        text_color = :black,
        fontsize = 28,
        background_color = (:white, 0.92),
        strokecolor = :darkgreen,
        strokewidth = 6,
        cornerradius = 8,
        padding = (12, 10, 12, 10),
    )

    CairoMakie.Colorbar(f[1, 2], p, label = "sensitivity (yr)", height = Relative(0.85))
    save("data/out/ex3/frequentist_sensitivity_map.png", f, px_per_unit = 2)
    f
end


# Now the same background overlay plot but with best ROI

let
    best_ROI = (best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge])

    # calculate signal activity for the optimized sensitivity and ROI
    best_0nu_activity = SNSensitivityEstimate.halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], best_sens*365*24*3600)
    SNSensitivityEstimate.set_activity!(signal, best_0nu_activity)

    bkg_processes = [p for p in processes if !p.signal]
    bkg_histograms = [SNSensitivityEstimate.get_bkg_counts_1D(p) for p in bkg_processes]
    signal_histogram = SNSensitivityEstimate.get_bkg_counts_1D(signal)
    colors = ["#041E42", "#BE4D00", "#951272", "#006630", "#005C8A", "#FFB948", "#605643", "#302D23"]
    f = CairoMakie.Figure(size = (1800, 800), fontsize = 24)
    a1 = CairoMakie.Axis(f[1, 1], title = "Background model with best signal overlay", xlabel = "sum of energies", ylabel = "expected counts", yscale = log10)
    a2 = CairoMakie.Axis(f[1, 2], title = "Zoom to best ROI", xlabel = "sum of energies", ylabel = "expected counts")

    for axis in (a1, a2)
        for i in length(bkg_histograms):-1:1
            CairoMakie.stephist!(axis, sum(bkg_histograms[1:i]), color = colors[i], linewidth = 5)
        end
        CairoMakie.stephist!(axis, signal_histogram, color = :black, linewidth = 4, linestyle = :dash)
        CairoMakie.vspan!(axis, best_ROI[1], best_ROI[2], color = (:darkgreen, 0.10))
    end

    labels = [String(p.isotopeName) for p in bkg_processes]
    push!(labels, "signal")

    elements = [CairoMakie.LineElement(color = colors[i], linewidth = 5) for i in 1:length(bkg_processes)]
    push!(elements, CairoMakie.LineElement(color = :black, linestyle = :dash, linewidth = 5))

    ylims!(a1, 1e-3, 1e5)
    ylims!(a2, 0, 5)
    xlims!(a2, best_ROI[1], best_ROI[2])
    CairoMakie.textlabel!(
        a2,
        best_ROI[1],
        4.6,
        text = "best ROI = ($(best_ROI[1]), $(best_ROI[2])) keV\nT_1/2 = $(ThalfbbESum) yr",
        text_align = (:left, :top),
        offset = (100, 0),
        text_color = :black,
        fontsize = 24,
        background_color = (:white, 0.92),
        strokecolor = :darkgreen,
        strokewidth = 4,
        cornerradius = 8,
        padding = (10, 8, 10, 8),
    )

    CairoMakie.Legend(f[1,3], elements, labels, "Processes", patchsize = (45, 20))
    save("data/out/ex3/best_roi_plot.png", f, px_per_unit = 2)
    f
end

