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
