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
using SNSensitivityEstimate, UnROOT, FHist, CairoMakie, CSV, DataFramesMeta, Metaheuristics
using Distributions 


# STEP 1: Load the data

# Start is the same as in ex1_background.jl, we load data and create DataProcess objects. 
data_info = CSV.read("data/data_info.csv", DataFrame)
file_paths = joinpath.("data", data_info.file_name) 
data_files = [UnROOT.ROOTFile(file_path) for file_path in file_paths] 
tree_name = "tree" 
variables = keys(data_files[1][tree_name]) 
data_tables = [UnROOT.LazyTree(data_file, tree_name, variables) for data_file in data_files] 


# STEP 1: define variables and their ranges that we are going to work with

# First we define variables we want to cut
var_names = [ 
    "sumE", 
    "dy", 
    "dz",
    "lPint",
    "lPext",
    ]

# Next we define the boundaries of the variables in the search space. These boundaries will be used to define the search space for the optimization algorithm.
var_bounds = (
    sumE = (300, 3500),
    dy = (0, 200),
    dz = (0, 200),
    lPint = (0., 100.0),
    lPext = (0, 100.0),
)

selected_tables = [UnROOT.LazyTree(data_file, tree_name, var_names) for data_file in data_files] 


# STEP 2: Create DataProcessND objects
# We do not filter data here, because that is the role of the optimization process. 
# In the ND example we no longer work with DataProcess object, but rather DataProcessND objects, which are designed to work with N-dimensional data vectors.

processes_nd = [
    SNSensitivityEstimate.DataProcessND(
        selected_tables[i], # Here we pass the whole data-set
        String(data_info.process[i]), 
        data_info.is_signal[i], 
        data_info.activity[i], 
        data_info.time_s[i], 
        data_info.n_sim[i], 
        var_bounds, # bounds for each variable
        data_info.amount[i], 
        var_names, # variable names
    ) for i in 1:length(selected_tables)
]


# STEP 3: define signal and background processes
signal_name = "bb0nu_foil_bulk"
signal_process = get_process(signal_name, processes_nd)[1]

background_names = data_info.process[data_info.process .!= signal_name]
background_process = [SNSensitivityEstimate.get_process(String(name), processes_nd)[1] for name in background_names]



# Test simple roi, you can play around by hand to get the feeling.
# NOTE: the ROI must specify a window for EVERY variable in `var_names`.
my_roi = (
    sumE = (2000, 3100),
    dy = (0, 200),
    dz = (0, 200),
    lPint = (0.0, 50.0),
    lPext = (0.0, 50.0),
)

# there are 2 approximation methods "table" and "formula", 
# table is more precise but takes longer to calculate, 
# formula is fine for quick tests
my_sensitivity = get_sensitivityND(SNparams, α, processes_nd, my_roi; approximate="table") 


# print your result, it shows the sensitivity, efficiency, background and roi
print(my_sensitivity) 


# STEP 4: Now we can set up the optimization problem to find the best ROI that maximizes the sensitivity.

# STEP 4.1: Define the step sizes for each variable, this will define the resolution of the search space.
# Build the grid decoder from the signal process + per-variable step sizes.
# The grid extent comes from `var_bounds`; the resolution from `var_steps`.
grid = GridSpec(signal_process, var_steps)


# STEP 4.2: Define the optimization problem, we want to minimize the negative sensitivity, which is equivalent to maximizing the sensitivity.
prob(proposed_indices) = - get_s_to_b(SNparams, α, processes_nd, proposed_indices, grid; approximate="formula") # we use the "formula" approximation here because it is faster, but you can also use "table" if you want more precise results.


# The search range is in INDEX units: (0, nsteps) for each ROI edge.
# `make_stepRange(grid)` returns these integer bounds automatically.
searchRange = make_stepRange(grid)


# Lower and upper bounds for the optimizer, derived from the index search range.
lower_bound = [x[1] for x in searchRange] .|> float
upper_bound = [x[2] for x in searchRange] .|> float


##############################################
# Now we can set up the optimization problem and run it
##############################################

# first the optimization options, basically the hyperparameters of the optimization algorithm
# you can leave them empty - default values will be used,
# or you can try to tune yourself to see what the results give, basically you want to find a good balance between time and quality of results
# For sure set up at least a time_limit and verbosity to true to see progress
# The full list of parameters can be found in https://jmejia8.github.io/Metaheuristics.jl/stable/api/#Metaheuristics.Options
options = Options(;
    # x_tol = 1.0, # tolerance in position space
    # f_tol = 1e-5, # tolerance in function value
    time_limit = 10*20.0, # time limit in seconds, must include .0 to be a Float64 value, doesn't take integers
    verbose = true, # verbosity on
    # iterations = 15, # number of iterations / generations
)

# Set up boxconstraints bounds - just the lower and upper bounds we defined before
bounds = boxconstraints(lb = lower_bound, ub = upper_bound)

# Initial guess for the ROI, IN INDEX UNITS (not physical units).
# Convert a physical edge to an index with: index = (edge - min) / step.
# Order matches `var_names`: sumE, dy, dz, lPint, lPext.
# e.g. sumE in (2700, 3100) with min=300, step=100 -> indices (24, 28).
x0 = float.([
    value_to_index(2700, grid.mins[1], grid.steps[1]; nsteps=grid.nsteps[1]), value_to_index(3100, grid.mins[1], grid.steps[1]; nsteps=grid.nsteps[1]),  # sumE  -> (2700, 3100)
    value_to_index(0, grid.mins[2], grid.steps[2]; nsteps=grid.nsteps[2]), value_to_index(100, grid.mins[2], grid.steps[2]; nsteps=grid.nsteps[2]),   # dy    -> (0, 100)
    value_to_index(0, grid.mins[3], grid.steps[3]; nsteps=grid.nsteps[3]), value_to_index(140, grid.mins[3], grid.steps[3]; nsteps=grid.nsteps[3]),   # dz    -> (0, 140)
    value_to_index(0.0, grid.mins[4], grid.steps[4]; nsteps=grid.nsteps[4]), value_to_index(4.0, grid.mins[4], grid.steps[4]; nsteps=grid.nsteps[4]),   # lPint -> (0.0, 4.0)
    value_to_index(2.0, grid.mins[5], grid.steps[5]; nsteps=grid.nsteps[5]), value_to_index(110.0, grid.mins[5], grid.steps[5]; nsteps=grid.nsteps[5]),   # lPext -> (2.0, 110.0)
])


# Set up optimizer algorithm, I generally use ECA - supports parallel evaluation, PSO - also can use parallel, SA - single core, but pretty good
algo = PSO(;options)
# algo = ECA(;options)

# Set the initial guess for the optimizer
set_user_solutions!(algo, x0, prob)


# Now we run the optimization, depending on how much time you gave it in the options, this can take a while
result = optimize(prob, bounds, algo)

# We can show the best result found
@show minimum(result)

# and store its (index-space) value
@show res = minimizer(result)


# Decode the best step-indices back to a physical ROI (note the extra `grid`
# argument), then calculate the precise sensitivity with the "table" method.
# `res = minimizer(result)` is the best index-vector found above.
best = get_best_ROI_ND(res, signal_process, grid)
best_sens = get_sensitivityND(SNparams, α, processes_nd, best; approximate="table")

# finally print the result
println("\n===== BEST RESULT =====")
print(best_sens)


# PLOT results

let
    binning = 0:100:3500
    plot_var = :sumE
    best_ROI = best
    min_a, min_b = best_ROI[plot_var]
    best_ROI = merge(best_ROI, (sumE = (binning[1], binning[end]),))
    best_sens = get_sensitivityND(SNparams, α, processes_nd, best; approximate = "table").tHalf



    # calculate signal activity for the optimized sensitivity and ROI
    best_0nu_activity = SNSensitivityEstimate.halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], best_sens*365*24*3600)
    SNSensitivityEstimate.set_activity!(signal_process, best_0nu_activity)

    bkg_histograms = [SNSensitivityEstimate.get_roi_bkg_counts_hist(p, best_ROI, binning, plot_var) for p in background_processes]
    signal_histogram = SNSensitivityEstimate.get_roi_bkg_counts_hist(signal_process, best_ROI, binning, plot_var)
    colors = ["#041E42", "#BE4D00", "#951272", "#006630", "#005C8A", "#FFB948", "#605643", "#302D23"]
    f = CairoMakie.Figure(size = (1800, 800), fontsize = 24)
    a1 = CairoMakie.Axis(f[1, 1], title = "Background model with best signal overlay", xlabel = "sum of energies", ylabel = "expected counts", yscale = log10)
    a2 = CairoMakie.Axis(f[1, 2], title = "Zoom to best ROI", xlabel = "sum of energies", ylabel = "expected counts")

    for axis in (a1, a2)
        for i in length(bkg_histograms):-1:1
            CairoMakie.stephist!(axis, sum(bkg_histograms[1:i]), color = colors[i], linewidth = 5)
        end
        CairoMakie.stephist!(axis, signal_histogram, color = :black, linewidth = 4, linestyle = :dash)
        CairoMakie.vspan!(axis, min_a, min_b, color = (:darkgreen, 0.10))
    end

    labels = [String(p.isotopeName) for p in background_processes]
    push!(labels, "signal")

    elements = [CairoMakie.LineElement(color = colors[i], linewidth = 5) for i in 1:length(background_processes)]
    push!(elements, CairoMakie.LineElement(color = :black, linestyle = :dash, linewidth = 5))

    ylims!(a1, 1e-3, 1e5)
    ylims!(a2, 0, 5)
    xlims!(a2, min_a, min_b)
    CairoMakie.textlabel!(
        a2,
        min_a,
        4.6,
        text = "T_1/2 = $(round(best_sens, digits=3)) yr",
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
    save("/sps/nemo/scratch/mpetro/Sensitivity_training/data/out/ex4/best_roi_plot.png", f, px_per_unit = 2)
    f
end

