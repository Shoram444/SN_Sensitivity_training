import Pkg;
Pkg.activate(".")
using SNSensitivityEstimate, UnROOT, FHist, CairoMakie, CSV, DataFramesMeta

# STEP 1: Load the data
# STEP 1.1: We start by loading the data/data_info.csv file which holds the information about the data files we want to analyze. 
# This file contains the paths to the ROOT files and the corresponding necessary information to process them.

data_info = CSV.read("data/data_info.csv", DataFrame)

# STEP 1.2: We then load the data files specified in the data_info DataFrame.

# the . in julia is used to broadcast operations over arrays. 
file_paths = joinpath.("data", data_info.file_name) 

# The [x for x in y] in julia is a list comprehension that creates a new array by applying the expression x to each element in the array y.
data_files = [UnROOT.ROOTFile(file_path) for file_path in file_paths] 

# STEP 1.3: Turn UnROOT.ROOTFile objects to Tabular objects of LazyTree type
tree_name = "tree" # name of the root tree
variables = keys(data_files[1][tree_name]) # get the variable names from the first file

data_tables = [UnROOT.LazyTree(data_file, tree_name, variables) for data_file in data_files] 

# (OPTIONAL) STEP 1.4: We can inspect the loaded data by printing the first few rows of the first data table.
println(first(data_tables[1], 5)) # print the first 5 rows of the first data table

# (OPTIONAL) STEP 1.5: We can also visualize the distributions of selected variables and processes using the CairoMakie package. In this example, we will plot the distributions of the "sumE" variable for each process.:

let
    binning = 0:100:3500 # define the binning for the histogram in format start:step:end

    f = CairoMakie.Figure(size = (800, 600))
    a = CairoMakie.Axis(f[1, 1], title = "sumE Distributions", xlabel = "sum of energies", ylabel = "normalized counts")

    for (i, data_table) in enumerate(data_tables)
        process_name = data_info.process[i]
        sumE_values = data_table.sumE # extract the sumE values from the data table
        hist = FHist.Hist1D(sumE_values; binedges = binning) # create a normalized histogram
        hist = FHist.normalize(hist) # normalize the histogram
        CairoMakie.stephist!(a, hist, label = process_name) # plot the histogram as a line plot
    end

    CairoMakie.axislegend(a, position = :lt) # add a legend to the plot
    f
end

# STEP 2: Create SNSensitivtyEstimate.DataProcess - the basic data structure of the SNSensitivityEstimate package. This structure holds the data and the necessary information to perform the sensitivity analysis.

# The following fields are defined for the SNSensitivityEstimate.DataProcess structure:
#
#    + dataVector::Vector{<:Real} - vector of initial data points
#    + isotopeName::String - isotope name 
#    + signal::Bool - signal or background
#    + activity::Real - activity of the given process in units [Bq]
#    + timeMeas::Real - time of measurement in units [s]
#    + nTotalSim::Real - total number of originally simulated events from Falaise
#    + bins::AbstractRange - binning to be used in format (min:step:max)
#    + amount::Real - mass [kg] or volume [l / m³] of the object where the isotope is present: i.e. source foil
#    + efficiency::Hist2D - the 2D histogram with efficiencies per bin 

processes = [
    SNSensitivityEstimate.DataProcess(
        collect(data_tables[i].sumE), # Here we specify a single variable to be used. In later stages we will use N-dimensional data vectors.
        String(data_info.process[i]), 
        data_info.is_isgnal[i], 
        data_info.activity[i], 
        data_info.time_s[i], 
        data_info.n_sim[i], 
        0:100:3500, 
        data_info.amount[i], 
    ) for i in 1:length(data_tables)
]

# STEP 3: get a histogram of background counts of each process. This is done using the SNSensitivtyEstimate.get_bkg_counts_1D function

# To calculate the expected bkg counts we use the following equation:
# n_i = a_i * m_i * t * e_i
# Where i is the i-th bin for which we calculate the expected background counts, n_i is the expected number of background counts in units [counts],
# a_i is activity in units [Bq/x], 
# m_i is the amount X of the process in units [x], 
# t is the time of measurement in units [s], 
# e_i is the efficiency of the process in units [1/s]. 
# e_i is calculated as e_i = n_pass_i / n_sim, where n_pass_i is the number of events that pass the cuts in the i-th bin and n_sim is the total number of simulated events.

# we can define a simple function that calculates n_i by defining the edges of the i-th bin:
function get_n_i(process::SNSensitivityEstimate.DataProcess, bin_lower::Real, bin_high::Real)
    n_pass_i = 0
    for i in 1:length(process.dataVector)
        if process.dataVector[i] >= bin_lower && process.dataVector[i] < bin_high
            n_pass_i += 1
        end
    end
    e_i = n_pass_i / process.nTotalSim
    n_i = process.activity * process.amount * process.timeMeas * e_i
    return n_i
end

# Let's test n_i for bb_foil_bulk in bin defined by 2000 - 2100 keV
n_i_test = get_n_i(processes[5], 2000, 2100)
println("Expected background counts for bb_foil_bulk in bin 2000-2100 keV: ", round(n_i_test, digits=2), " counts")


# Repeating this process we could build a full histogram by finding out how many events we expect in each bin for each process.
# Luckily, the SNSensitivityEstimate package has a built-in function that does this for us. The function is called get_bkg_counts_1D and it takes a SNSensitivityEstimate.DataProcess object as input and returns a histogram of expected background counts for each bin.
bb_histogram = get_bkg_counts_1D(processes[5]) # get the background counts for the bb_foil_bulk process

# plot the histogram
let
    f = CairoMakie.Figure(size = (800, 600))
    a = CairoMakie.Axis(f[1, 1], title = "Expected Background Counts for bb_foil_bulk", xlabel = "sum of energies", ylabel = "expected counts")
    CairoMakie.stephist!(a, bb_histogram, label = "bb_foil_bulk")
    CairoMakie.axislegend(a, position = :lt)
    f
end

# Or even simpler, we can get all histograms at once by using the get_bkg_counts_1D function on the array of processes:
filter!(p -> !p.signal, processes) # filter out the signal processes
all_histograms = get_bkg_counts_1D.(processes)

# STEP 4: Plot the stacked histogram and play around with obtaining counts per process, saving to tables etc.

# STEP 4.1: Plot
let
    colors = ["#041E42", "#BE4D00", "#951272", "#006630", "#005C8A", "#FFB948", "#605643", "#302D23"] # get our supernemo colors
    f = CairoMakie.Figure(size = (1000, 600))
    a = CairoMakie.Axis(f[1, 1], title = "Stacked Histogram of Expected Background Counts", xlabel = "sum of energies", ylabel = "log of expected counts", yscale = log10)
    p = CairoMakie.hist!(a, sum(all_histograms), color = colors[1])
    for i in length(all_histograms):-1:1
        CairoMakie.hist!(a, sum(all_histograms[1:i]), color = colors[i])
    end

    labels = [String(p.isotopeName) for p in processes]
    elements = [CairoMakie.PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    title = "Processes"
    
    ylims!(a, 1e-3, 1e5)
    CairoMakie.Legend(f[1,2], elements, labels, title)
    f
end

# STEP 4.2: Get counts per process and save to table
counts_per_process = [integral(hist) for hist in all_histograms] # we integrate each histogram to get the total counts per process
df_counts_per_process = DataFrame(isotopeName = [String(p.isotopeName) for p in processes], counts = counts_per_process)
CSV.write("data/out/ex1/counts_per_process.csv", df_counts_per_process)

# STEP 4.3: Get counts in ROI
a,b = 2700, 3000
counts_per_process_roi = [integral(restrict(hist, a, b)) for hist in all_histograms] # we first restrict the histogram to the ROI and then integrate to get the total counts per process in the ROI
df_counts_per_process_roi = DataFrame(isotopeName = [String(p.isotopeName) for p in processes], counts = counts_per_process_roi)
CSV.write("data/out/ex1/counts_per_process_roi.csv", df_counts_per_process_roi)

