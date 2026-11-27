global signal_name = "bb0nuM1"
global phase = "phase_4"
global optimization_time_limit = 1000 # seconds; leave enough SLURM time for loading data
global output_dir = ENV["SN_OUTPUT_PATH"]
global α = 1.64485362695147
algorithm = "PSO" # PSO or ECA

# :all or a list of background names from the metadata CSV.
background_names = :all
# background_names = ["bb", "bi_foil"]

# Extra ROOT fields referenced by custom numeric cut formulas.
extra_var_names = String[]

# Fixed topological cuts.
use_same_calo_cut = true
use_kinks_cut = true
use_om_cut = true
use_vertex_region_cut = true
# Dead OMs: the time-dependent map (data/dead_om_map_all_phases.csv) vetoes an
# OM only while it was dead, in data and simulation alike. The static bad OM
# list (data/om_lists) removes its OMs for the whole phase; keep it off when
# the map is used. The OM cut then requires two main-wall OMs.
use_bad_om_list = false
use_dead_om_map = true

# The ND package applies all windows to abs(value): min <= abs(value) < max.
var_bounds_nd = (
    energy_elec_sum = (0.5, 3.50),
    angle_3D_between_ep_em = (0, 180),
    delta_y_elec = (0, 500),
    delta_z_elec = (0, 500),
    diff_time_elec = (0, 10),
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
    energy_elec_sum = (1.15, 2.850),
    angle_3D_between_ep_em = (0, 180),
    delta_y_elec = (0, 200),
    delta_z_elec = (0, 200),
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
        bounds = (10.0, 600.0),
        step = 10.0,
    ),
    (
        name = :closest_gamma,
        value = event -> abs(event.closest_gamma),
        comparison = (>),
        mode = :optimize,
        fixed = 50.0,
        bounds = (0.0, 500.0),
        step = 10.0,
    ),
    (
        name = :closest_elec,
        value = event -> abs(event.closest_elec),
        comparison = (>),
        mode = :optimize,
        fixed = 10.0,
        bounds = (0.0, 500.0),
        step = 5.0,
    ),
    (
        name = :closest_track,
        value = event -> abs(event.closest_track),
        comparison = (>),
        mode = :optimize,
        fixed = 200.0,
        bounds = (0.0, 800.0),
        step = 10.0,
    ),
    (
        name = :closest_time_track,
        value = event -> abs(event.closest_time_track),
        comparison = (>),
        mode = :optimize,
        fixed = 10.0,
        bounds = (0.0, 200.0),
        step = 5.0,
    ),
]

# Added by run2.sh --use_dead_om_map=true
global use_dead_om_map = true
