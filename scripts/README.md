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
```

This tests generation and mocked submission without SLURM or production ROOT
files. It is not a cluster or full physics-analysis test.

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
