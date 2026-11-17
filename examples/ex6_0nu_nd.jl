#### 
#### Example script for running in batch mode - calculates sensitivity for 0nu for phase 3
####

global PROJECT_DIR="/sps/nemo/scratch/mpetro/Sensitivity_training/" # Change this to your project directory

# Select which phase we process
global signal_name = "bb0nu"
global phase = "phase_3"
global optimization_time_limit = 1 * 3600  # make sure your optimization time limit is about 1hr less than job time limit to account for overhead tasks
global output_dir = joinpath(PROJECT_DIR, "data", "out", "ex6_$(phase)_results_1.csv")
global α = 1.64485362695147


@info "Running sensitivity optimization for phase: $phase, with signal: $signal_name, output will be saved to: $output_dir"
flush(stdout)

using SNSensitivityEstimate, UnROOT, FHist, CSV, DataFramesMeta
using Distributions, Metaheuristics


@info "Packages loaded, loading data..."
flush(stdout)

# STEP 1: Load the data
# Start is the same as in ex1_background.jl, we load data and create DataProcess objects. 
data_info = CSV.read(joinpath(PROJECT_DIR, "data", "mathis_data_info_final.csv"), DataFrame)


selected_files = @chain data_info begin
    @subset :phase .== phase
    @subset .!occursin.("reduces", :file_path) 
end

data_files = [UnROOT.ROOTFile(joinpath(PROJECT_DIR, file_path)) for file_path in selected_files.file_path] 
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



# Sensitivity ND
# define variables to use
var_names_nd = [ 
    "energy_elec_sum", 
    "angle_3D_between_ep_em", 
    "delta_y_elec",
    "delta_z_elec",
    "diff_time_elec",
    "energy_elec_1",
    "energy_elec_2",
    ]

# Next we define the boundaries of the variables in the search space. These boundaries will be used to define the search space for the optimization algorithm.
var_bounds_nd = (
    energy_elec_sum = (2.40, 3.50),
    angle_3D_between_ep_em = (0, 180),
    delta_y_elec = (0, 200),
    delta_z_elec = (0, 200),
    diff_time_elec = (-5, 5),
    energy_elec_1 = (0, 3.50),
    energy_elec_2 = (0, 3.50),
)


# Default search resolution (adjust if needed)
var_steps_nd = (
    energy_elec_sum = 0.05, # 50keV
    angle_3D_between_ep_em = 5.0, # 5 degree
    delta_y_elec = 5.0, # 5mm
    delta_z_elec = 5.0, # 5mm
    diff_time_elec = 0.25, # 0.25ns
    energy_elec_1 = 0.05, # 50keV
    energy_elec_2 = 0.05, # 50keV
)

selected_tables_nd = [UnROOT.LazyTree(data_file, tree_name, var_names) for data_file in data_files] 
filtered_data_tables_nd = [SNSensitivityEstimate.filter_data(table, var_names_nd, var_bounds_nd) for table in selected_tables_nd]

@info "Data loaded and filtered, applying Mathis's extra cuts..."
flush(stdout)

# Apply Mathis's extra cuts

# Apply additional topological cuts and store filtered tables back.
for i in eachindex(filtered_data_tables_nd)
    table = filtered_data_tables_nd[i]
    cut_vertex_dist = sqrt.(table.delta_y_elec .* table.delta_y_elec .+ table.delta_z_elec .* table.delta_z_elec) .>= 300
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
    filtered_data_tables_nd[i] = table[keep]
end

@info "Mathis's extra cuts applied, creating ND processes..."
flush(stdout)

# create ND processes

# STEP 2: Create DataProcessND objects
# We do not filter data here, because that is the role of the optimization process. 
# In the ND example we no longer work with DataProcess object, but rather DataProcessND objects, which are designed to work with N-dimensional data vectors.

processes_nd = [
    SNSensitivityEstimate.DataProcessND(
        filtered_data_tables_nd[i], # Here we pass the whole data-set
        String(selected_files.name[i]), 
        selected_files.is_signal[i], 
        selected_files.activity_bq[i], 
        selected_files.time_s[i], 
        selected_files.nb_events[i], 
        var_bounds_nd, # bounds for each variable
        data_info.amount[i], 
        var_names_nd, # variable names
    ) for i in 1:length(selected_tables_nd)
]

# decalre what's signal and background

signal_process = get_process(signal_name, processes_nd)[1]

background_names = data_info.name[data_info.name .!= signal_name]
background_processes = [SNSensitivityEstimate.get_process(String(name), processes_nd)[1] for name in background_names]


# define the optimization problem

@info "Defining the optimization problem..."
flush(stdout)

grid = GridSpec(signal_process, var_steps_nd)


# STEP 4.2: Define the optimization problem, we want to minimize the negative sensitivity, which is equivalent to maximizing the sensitivity.
prob(proposed_indices) = - get_s_to_b(SNparams, α, processes_nd, proposed_indices, grid; approximate="formula") # we use the "formula" approximation here because it is faster, but you can also use "table" if you want more precise results.


# The search range is in INDEX units: (0, nsteps) for each ROI edge.
# `make_stepRange(grid)` returns these integer bounds automatically.
searchRange = make_stepRange(grid)


# Lower and upper bounds for the optimizer, derived from the index search range.
lower_bound = [x[1] for x in searchRange] .|> float
upper_bound = [x[2] for x in searchRange] .|> float

options = Options(;
    # x_tol = 1.0, # tolerance in position space
    # f_tol = 1e-4, # tolerance in function value
    time_limit = float(optimization_time_limit), # time limit in seconds, must include .0 to be a Float64 value, doesn't take integers
    verbose = true, # verbosity on
    # iterations = 1000, # number of iterations / generations
)

# Set up boxconstraints bounds - just the lower and upper bounds we defined before
bounds = boxconstraints(lb = lower_bound, ub = upper_bound)

algo = PSO(;options)

# initial condition, adjust according to your ideas

x0 = float.([
    value_to_index(2.700, grid.mins[1], grid.steps[1]; nsteps=grid.nsteps[1]), value_to_index(3.100, grid.mins[1], grid.steps[1]; nsteps=grid.nsteps[1]),  # sumE  -> (2.700, 3.100)MeV
    value_to_index(0, grid.mins[2], grid.steps[2]; nsteps=grid.nsteps[2]), value_to_index(180, grid.mins[2], grid.steps[2]; nsteps=grid.nsteps[2]),   # phi    -> (0, 180)deg
    value_to_index(0, grid.mins[3], grid.steps[3]; nsteps=grid.nsteps[3]), value_to_index(100, grid.mins[3], grid.steps[3]; nsteps=grid.nsteps[3]),   # dy    -> (0, 100)mm
    value_to_index(0, grid.mins[4], grid.steps[4]; nsteps=grid.nsteps[4]), value_to_index(140, grid.mins[4], grid.steps[4]; nsteps=grid.nsteps[4]),   # dz    -> (0, 140)mm
    value_to_index(-2.5, grid.mins[5], grid.steps[5]; nsteps=grid.nsteps[5]), value_to_index(2.5, grid.mins[5], grid.steps[5]; nsteps=grid.nsteps[5]),   # dT    -> (-2.5, 2.5)ns
    value_to_index(0.300, grid.mins[6], grid.steps[6]; nsteps=grid.nsteps[6]), value_to_index(2.900, grid.mins[6], grid.steps[6]; nsteps=grid.nsteps[6]),  # E1  -> (0.3, 3.100)MeV
    value_to_index(0.300, grid.mins[7], grid.steps[7]; nsteps=grid.nsteps[7]), value_to_index(2.900, grid.mins[7], grid.steps[7]; nsteps=grid.nsteps[7]),  # E2  -> (0.3, 3.100)MeV
])

set_user_solutions!(algo, x0, prob)

# optimize and wait for convergence or time limit

@info "Starting optimization with PSO algorithm..."
flush(stdout)

result = optimize(prob, bounds, algo)

@info "Optimization finished, processing results..."
flush(stdout)

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

# Save results to csv
df = DataFrame(
    phase = phase,
    best_sensitivity = [best_sens.tHalf],
    best_ROI = [best_sens.roi],
    best_bkg = [best_sens.bkgCounts],
    best_eff = [best_sens.signalEff],
)
CSV.write(output_dir, df)