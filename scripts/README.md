# ND sensitivity runs

`run2.sh` prepares and submits one SLURM job for one configuration. It uses
[`run_nd_sensitivity_on_mathis_data_v2.jl`](../analysis/run_nd_sensitivity_on_mathis_data_v2.jl)
as a template; it does not modify the template itself.

## Preparation

- The cluster needs SLURM, Julia, and the dependencies in `Project.toml` and
  `Manifest.toml`. ROOT files are read through UnROOT.
- The project must contain the files referenced by
  [`mathis_data_info_final_2.csv`](../data/mathis_data_info_final_2.csv) and
  the phase-specific lists in [`data/om_lists/`](../data/om_lists/).
- Copy [`nd_sensitivity.example.conf`](../configs/nd_sensitivity.example.conf)
  and edit the copy for your run. The project root is determined from the
  location of `run2.sh`; run the script from your cluster clone.
- Adjust the `RUN_` settings to the cluster. `RUN_MODULE` loads the Julia
  module; leave `RUN_SYSIMAGE` empty to run without a custom sysimage. The
  configured sysimage must be compatible with the project's package versions.

## Submit

From the project root:

```bash
# Create the run directory and preview the sbatch command; do not submit.
bash scripts/run2.sh configs/nd_sensitivity.example.conf --dry-run

# Create a separate run directory and submit it.
bash scripts/run2.sh configs/nd_sensitivity.example.conf
```

Use the path to your edited configuration in place of the example. The script
prints the run directory and, on submission, the SLURM job ID. Use `squeue -j
JOB_ID` to check a running job. Inspect the error log if the job fails.

By default, each invocation creates a unique directory in `results/nd_runs/`.
The dry run creates the directory and files too, but neither submits nor runs
the analysis. Submitted runs contain:

- `config.conf`: copy of the supplied configuration.
- `settings.jl`: extracted Julia configuration.
- `analysis.jl`: generated script, adapted from the template using `sed`.
- `job.sh`: generated batch script; `job_id.txt`: submitted job ID.
- `slurm-JOBID.out` and `slurm-JOBID.err`: SLURM logs.
- `results.csv`: sensitivity and best ROI after a successful run.

The project environment, input ROOT files, metadata and sysimage are not
snapshotted. Keep them unchanged if runs must remain comparable.

### Optional background report

Set `# RUN_SAVE_AUXILIARY_DATA=true` in the configuration (`false` by default).
The runner copies `auxiliary.jl` into the run directory and the analysis
includes it after successful optimization, in the same Julia process.
The report uses the best ROI and the preselected background columns already
in memory. It does not reopen ROOT files or repeat the fixed cuts. With this
option enabled, the analysis retains background events over 0–4 MeV for the
report, so it can use more memory; the optimization keeps its original
energy-sum preselection.

- `background_spectrum.png`: logarithmic step plot of expected counts for the
  full exposure, in 50 keV bins from 0 to 4 MeV. The best energy-sum ROI is
  shaded, not applied as a cut; the configured optimization energy-sum
  prefilter is also lifted. All other cuts remain active. PMTs are grouped
  as in the ex5 plot. Curves drop to the bottom of the log axis at zero-count
  bins. Vertical bars show ±1σ MC statistical uncertainty, clipped at the
  bottom of the axis where necessary. A vector copy is saved as
  `background_spectrum.pdf`. The subtitle shows the signal used to optimize
  the ROI and the exposure in kg·yr (process live time × foil mass).
- `background_components.csv`: passing MC events and expected counts per
  background (plus a total) for the full best ROI and for
  `0 <= E_sum < 4 MeV` with only the energy-sum window replaced. The ROI total
  is checked against `best_bkg` in `results.csv`.

The plot and CSV use the weighted Poisson estimate `sqrt(sum(w^2))` for MC
statistical uncertainty, where `w` is each passing event's normalization
weight. Independent contributions add in quadrature, including PMT groups
and totals. Activity and other systematic uncertainties are not included.

## Configuration

The `# RUN_...=` lines configure SLURM and output paths. The section between
`# BEGIN JULIA CONFIG` and `# END JULIA CONFIG` is **executable Julia code**:
use only trusted configuration files.

- Set `signal_name`, `phase` and `background_names` (`:all` or a list of names
  in the metadata). Exactly one signal is used.
- Choose `algorithm = "PSO"` or `"ECA"`. `optimization_time_limit` is in
  seconds and defaults to **10 seconds** for testing. Set `RUN_WALLTIME`
  longer to allow for data loading and final result calculation.
- `var_bounds_nd` and `var_steps_nd` define the search grid. `initial_roi`
  defines one starting ROI; its edges must be inside the bounds and on the grid.
- Switch the topological cuts with the `use_..._cut` booleans. The OM list is
  still loaded at startup even when `use_om_cut` is false.
- In `nd_numeric_cuts`, `mode = :fixed` applies the configured comparison and
  `fixed` threshold before optimization. `mode = :optimize` searches the
  threshold over `bounds` and `step`; `fixed` becomes its initial value and
  must lie on that grid. Add ROOT fields needed by custom expressions to
  `extra_var_names`.

ND windows use `min <= abs(value) < max`, including for timing. Optimized
numeric cuts also use absolute values and the package's `<` or `>=` boundary:
an optimized `<=` behaves as `<`, and `>` behaves as `>=`. Fixed numeric cuts
use their configured comparison exactly. The script warns about the changed
boundary for optimized `>` and `<=`.

Normalization remains the full exposure from `SNparams["t"]`, as in the
original batch script. Phase 4 uses the paths recorded in the metadata CSV;
the OM list is selected by the configured data phase.

## Tests

```bash
bash tests/test_run2.sh
julia --project=. tests/test_nd_background_report.jl
```

These tests cover generation, mocked submission, and the background report on
synthetic events without SLURM or production ROOT files.

## Listing ROOT event parameters

On the cluster, load ROOT and list the `last_event_number` parameter for each
distinct path in the **original** `data/mathis_data_info_final.csv`:

```bash
module load root
bash scripts/list_root_events.sh > root_events.csv
```

Python 3 is required to read the CSV. Paths are relative to the project root.
Errors go to stderr; affected rows have an empty `nb_events` and the command
exits nonzero. The output is a separate report and does not change either
metadata CSV. The parameter is copied exactly, without adding one or reading
tree entries; its count-versus-last-ID convention still needs confirmation.

## Data vs simulation

### 1. Run times (cluster)

Each run has its own directory `run_XXX` containing `extracted_data.root`,
which stores the run duration in the `run_time` parameter. On the cluster, list
it for every `run_*/extracted_data.root` in the default directory
`/sps/nemo/scratch/granjon/double_beta/cimrman_curve/ptd_pol2/root_with_max_threshold_working`:

```bash
module load root
bash scripts/list_run_times.sh > data/real_data/run_times.csv
# or another directory with the same layout:
bash scripts/list_run_times.sh DIRECTORY > data/real_data/run_times.csv
```

No other files are read. The CSV has `file_path,run_number,run_time`; the run
number comes from the directory name (leading zeros removed). A warning is
printed when the first event's `run_number` differs, and for `run_*`
directories without `extracted_data.root`. Errors go to stderr and the
command exits nonzero; unreadable files keep their run number but have an
empty `run_time`. The comparison ignores such runs with a warning unless they
are selected (a data run, or within `--runs=FIRST:LAST`), and rejects run
numbers listed more than once.

### 2. Comparison

```bash
julia --project=. scripts/compare_data_simulation.jl \
    --results=results/nd_runs/<run>/results.csv
```

The script reads the best ROI from `results.csv` and the fixed cuts
(topological switches and `:fixed` numeric cuts) from `settings.jl` next to
it. It applies the same cuts as the analysis to the phase backgrounds from the
metadata CSV and to the data. The energy-sum window is not applied; it is
shaded, as in the background report. Simulations are normalized to counts/day
with `amount * activity_bq * 86400 / nb_events` per event, data with the summed
`run_time` (assumed in seconds) of the selected runs.

Options (`--name=value`; relative paths are relative to the project root):

- `--settings`: settings file [`settings.jl` next to the results].
- `--metadata`: [`data/mathis_data_info_final_2.csv`].
- `--data`: [`data/real_data/run_3470-4427.root`].
- `--run-times`: [`data/real_data/run_times.csv`].
- `--runs`: `data` uses the runs present in the data file [default];
  `FIRST:LAST` keeps only data events from those runs and sums the run time
  of all listed runs in the range (warning about runs without events).
- `--sim-dir`: replaces the `data/mathis_files` prefix of metadata paths,
  e.g. `--sim-dir=data/mathis_files_local` for a local copy made with
  `./data/download_mathis_files.sh [PHASE ...]`. On the cluster,
  `data/mathis_files` is an untracked symlink to the simulation directory;
  create it once with `ln -s <simulation directory> data/mathis_files`.
- `--output-dir`: [directory of the results file].
- `--edges`: energy-sum binning `MIN:STEP:MAX` in MeV [`0:0.05:4`].
- `--chi2-range`: `MIN:MAX` in MeV used for χ² [the whole binning].

Outputs:

- `data_vs_simulation.png`/`.pdf`: simulated components (PMTs grouped), total
  simulation (red) with ±1σ MC statistical band, and data (black) with central
  68% Poisson (Garwood) intervals. The legend shows χ²/ndf;
  below, each bin's residual `sign(n − μ) · sqrt(χ²_bin)`, whose squares sum
  to χ². Bins outside the χ² range are faded; residuals beyond the axis are
  drawn as triangles.
- `data_vs_simulation_bins.csv`: per-bin rates, errors and residuals.
- `data_vs_simulation_summary.csv`: events and counts/day per background,
  total simulation and data, over the full binning and in the energy ROI.
- `data_vs_simulation_stats.csv`: live time, runs, sensitivity, χ²/ndf and
  p-value. `chi2` is the Poisson likelihood ratio (Baker–Cousins)
  `Σ 2(μ − n + n ln(n/μ))`, with each bin's simulation allowed to vary within
  its MC statistical error (Barlow–Beeston lite), so it remains valid for bins
  with few or no events. Nothing is fitted, so `ndf` is the number of bins with
  `μ > 0`, and `p_value` is the χ²(ndf) probability of exceeding `chi2`.
  `chi2_pearson` (`Σ (n − μ)² / (μ + σ²_MC)`) is kept for reference; it is
  unreliable for bins expecting far below one event. Bins with data but no
  simulation are excluded and counted. Bins expecting almost nothing add to
  `ndf` but hardly to χ²; restrict them with `--chi2-range` if needed.

The script warns if the simulated ROI background times `SNparams["t"]` does not
reproduce `best_bkg` in `results.csv`, e.g. after metadata changes.

```bash
julia --project=. tests/test_compare_data_simulation.jl
```
