global signal_name = "bb0nuM2"
global phase = "phase_4"
global optimization_time_limit = 48*3600.0 # seconds; leave enough SLURM time for loading data
global output_dir = ENV["SN_OUTPUT_PATH"]
global α = 1.64485362695147
algorithm = "PSO" # PSO or ECA

# :all or a list of background names from the metadata CSV.
background_names = :all
# background_names = ["bb", "bi_foil"]

# Extra ROOT fields referenced by custom numeric cut formulas.
extra_var_names = String[]

# Fixed topological cuts. Turning off the OM cut still requires the phase's
# OM list to be present because the script loads it at startup.
use_same_calo_cut = true
use_kinks_cut = true
use_om_cut = true
use_vertex_region_cut = true

# The ND package applies all windows to abs(value): min <= abs(value) < max.
var_bounds_nd = (
    energy_elec_sum = (1.0, 3.50),
    angle_3D_between_ep_em = (0, 180),
    delta_y_elec = (0, 200),
    delta_z_elec = (0, 200),
    diff_time_elec = (0, 5),
    energy_elec_1 = (0, 3.50),
    energy_elec_2 = (0, 3.50),
)

var_steps_nd = (
    energy_elec_sum = 0.05,
    angle_3D_between_ep_em = 5.0,
    delta_y_elec = 5.0,
    delta_z_elec = 5.0,
    diff_time_elec = 0.25,
    energy_elec_1 = 0.05,
    energy_elec_2 = 0.05,
)

# Each initial edge must be within the bounds and on the corresponding grid.
initial_roi = (
    energy_elec_sum = (2.200, 2.850),
    angle_3D_between_ep_em = (0, 180),
    delta_y_elec = (0, 65),
    delta_z_elec = (0, 65),
    diff_time_elec = (0, 1.75),
    energy_elec_1 = (0.300, 2.900),
    energy_elec_2 = (0.300, 2.900),
)

# :fixed uses comparison(value, fixed). :optimize uses fixed as its initial
# threshold and searches bounds in steps of step. Optimized cuts operate on
# abs(value) with < or >=; > and <= therefore differ at the boundary.
# Add new ROOT fields used by value to extra_var_names above.
nd_numeric_cuts = [
    (
        name = :vertex_distance,
        value = event -> sqrt(event.delta_y_elec^2 + event.delta_z_elec^2),
        comparison = (<),
        mode = :optimize,
        fixed = 300.0,
        bounds = (50.0, 300.0),
        step = 50.0,
    ),
    (
        name = :closest_gamma,
        value = event -> abs(event.closest_gamma),
        comparison = (>),
        mode = :optimize,
        fixed = 50.0,
        bounds = (20.0, 150.0),
        step = 10.0,
    ),
    (
        name = :closest_elec,
        value = event -> abs(event.closest_elec),
        comparison = (>),
        mode = :optimize,
        fixed = 25.0,
        bounds = (5.0, 125.0),
        step = 10.0,
    ),
    (
        name = :closest_track,
        value = event -> abs(event.closest_track),
        comparison = (>),
        mode = :fixed,
        fixed = 200.0,
        bounds = (100.0, 400.0),
        step = 10.0,
    ),
    (
        name = :closest_time_track,
        value = event -> abs(event.closest_time_track),
        comparison = (>),
        mode = :optimize,
        fixed = 10.0,
        bounds = (5.0, 50.0),
        step = 5.0,
    ),
]
