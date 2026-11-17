# Sensitivity jobs

## Status

Runner/submission and metadata-extraction implementation.
The required `mathis_analysis_info.csv` has **not** been generated.
Do not substitute the original CSV: its event totals have not been verified.

## Commands

Set the repository location explicitly; the working directory does not matter:

```sh
PROJECT_DIR=/path/to/SN_Sensitivity_training
```

1. Preview SLURM commands without reading metadata or ROOT events or writing files:

   ```sh
   julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/submit.jl" "$PROJECT_DIR/scripts/job.example.toml" --dry-run
   ```

2. On the cluster, create the second CSV from ROOT parameters only:

   ```sh
   module load root
   julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/prepare_metadata.jl" "$PROJECT_DIR/scripts/job.example.toml"
   ```

   Add `--dry-run` to preview selected paths locally without opening ROOT files.
   The source defaults to the original metadata CSV; use `--source CSV` for a
   reviewed alternative. Only configured backgrounds and the selected signal
   in the selected phase are extracted, excluding `reduces_files` rows.

3. Check metadata and input-file existence only:

   ```sh
   julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/run.jl" "$PROJECT_DIR/scripts/job.example.toml" --check
   ```

4. On the cluster, submit the configured analysis:

   ```sh
   julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/submit.jl" "$PROJECT_DIR/scripts/job.example.toml"
   ```

For a direct cluster run, use `run.jl CONFIG` without `--check`.
No entry point instantiates packages, downloads data, or modifies the source CSV.

## Configuration and normalization

- Paths resolve relative to the repository. SLURM CPUs set Julia threads.
  An empty module or sysimage setting disables that feature.
- Simulation uses **2.86 years / 90192960 seconds**, phase-3 conditions, and
  the package's isotope mass, enrichment, and Feldman-Cousins methods.
  A copied physics dictionary keeps process times and sensitivity normalization
  synchronized without changing package globals.
- Validation is off by default. When enabled, simulations are divided by their
  simulated exposure in days and measured counts by measured phase-3 live days.
  The measured duration remains provisional.
- Named cut sets separate fixed topology from variable windows. Set a variable's
  `role` to `fixed` to use its bounds without optimizing it. Otherwise supply
  `step` and `initial_roi`. Bounds and starting edges must lie on the grid.
- All ND variable windows apply to **absolute values**, including timing, both
  during preselection and optimization: lower inclusive, upper exclusive.
  Topological threshold names explicitly distinguish `gt`, `ge`, `lt`, and `le`.

The analysis CSV needs `phase,name,is_signal,file_path,activity_bq,amount,nb_events`
and `nb_events_source`, with the latter equal to `ROOT:last_event_number`.
It must have unique selected processes and files, exactly the requested signal,
and background groups covering every selected background once. Other phases and
other signals are not loaded. Missing inputs or unverified totals cause an error.
See [activity_comparison.csv](activity_comparison.csv) before changing activities;
its candidates use notebook defaults, not the optional fitted PMT factors.

The extractor copies positive integral `last_event_number` values exactly:
no `+1`, no selected-tree entry counts, no event-tree reads. The producer's
count-versus-last-ID convention still needs confirmation against simulation
generation records. It accepts integer ROOT parameters and exact integer-valued
double parameters; missing, fractional, unsupported, or nonpositive values fail.
Original counts, activities, and measured phase times are retained alongside
source-CSV hash, extraction timestamp, and ROOT file size/mtime. Analysis preflight
rejects changed ROOT files when these fingerprints are present.

Existing output CSVs are protected. `--overwrite` explicitly regenerates the
second CSV from the chosen source and replaces any manual edits to that output.
The source is never overwritten, and extraction errors do not publish partial
metadata. Regenerate metadata when changing the configured phase or signal.

## Results

Each submission gets a unique directory under the configured output root:
resolved configuration, SLURM script/log, and a `results/` directory.
The runner records metadata, hashes, software/job provenance, analysis logs,
and `run.toml` status (`running`, `completed`, or `failed`).
Partial outputs in a failed run are not final results.

Requested CSVs contain sensitivity/FC limit/exposure, ROI edges and units,
per-component/group backgrounds with generated/passing MC counts and
`sqrt(sum(w^2))` errors, and cumulative cutflow. Text summaries and PDF/PNG plots
are configurable. Plot CSVs accompany plots. Spectra apply all optimized cuts
except the plotted variable; shaded ROI shows its final window. Signal is scaled
to the final FC limit, not the placeholder signal activity in metadata.
Sensitivity slices change one ROI edge while holding all others fixed; they
are not re-optimized profiles. Validation output is a comparison, not a pass/fail
claim, and is never marked passed merely because a file was produced.

## Small-scale tests

```sh
julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/test_workflow.jl"
julia --project="$PROJECT_DIR" "$PROJECT_DIR/scripts/test_metadata.jl"
```

Tests use synthetic in-memory events, temporary metadata, and small plots.
They do not load real ROOT files or submit jobs.
When ROOT is available, metadata tests additionally create temporary parameter-only
ROOT files to test the real C++ reader. Otherwise that integration check is explicitly
skipped; run it on the cluster before extracting production metadata.

## Bounded cluster smoke test

Run this before production metadata extraction or submission. Dependencies must
already be available in the cluster project environment. Submit from the repository
root, like the existing batch workflow.
This requests one CPU and 8 GB for at most 10 minutes, excluding queue wait;
adjust partition/license names if your cluster setup differs.

```sh
cd "$PROJECT_DIR"
sbatch scripts/smoke_test.sh
```

The batch script loads Julia and ROOT, uses the project environment without a
custom sysimage, and calls Julia directly; it does not use `srun`.
Submission returns a job ID; output goes to `smoke-JOBID.log` in the repository
root. Inspect it with `tail -n 80 smoke-JOBID.log`, replacing `JOBID` with that ID.
If the log is not present yet, check `squeue -j JOBID` for a pending job.

Expected: **45 metadata checks and 52 workflow checks pass**, with no skipped ROOT
check. Missing/fractional-parameter error messages are intentional negative tests;
the final test summaries determine success. Tests use temporary parameter-only
ROOT files and synthetic events, not production event trees. Outputs are temporary.
This checks ROOT extraction, optimizer/package compatibility, tables, plots,
validation normalization, and CLI preflight, not real detector-data correctness.

If either suite fails or the allocation times out, stop and retain the first error;
do not submit production. Once both pass, follow Commands 2-4 above:
extract actual metadata, run preflight, then submit. The optimizer's 10-second
limit does **not** limit full-file loading, plotting, or total job duration.
