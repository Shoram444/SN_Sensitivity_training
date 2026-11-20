# Standardized sensitivity jobs

Next session: read [HANDOFF.md](HANDOFF.md); inspect the cluster batch smoke-test log first.

Status: Steps 1-5 complete locally; cluster smoke test/extraction pending.
Verified: 52 workflow + 42 metadata checks pass; no real ROOT files loaded.
ROOT/C++ integration is explicitly skipped locally; run it on the cluster.

## Steps

1. **Compare inputs** (~30 min): sensitivity-relevant notebook differences. Done.
2. **Define configuration** (~20 min): [job.example.toml](job.example.toml). Done.
3. **Build the workflow** (~60-90 min): Julia analysis, SLURM submission, strict
   configuration checks, no-data dry-run, unique run directories. Done.
4. **Reconcile metadata and validation** (~45-60 min): read `last_event_number`
   into a second analysis CSV, preserving source metadata and extraction provenance.
   Extractor and optional data/MC checks ready; actual CSV awaits cluster extraction.
   Activities still await colleague review.
5. **Verify cluster use** (~20-30 min): local syntax/configuration and synthetic
   tests done; [README](README.md) has a 10-minute cluster smoke test for you to run.

## Decisions

- Phase-3 detector conditions for 2.86 years = 90192960 seconds (~17.5 kg yr).
  Simulated process exposure and sensitivity normalization must use this same time.
  Data validation compares counts/day: simulation / 1043.9 days, data / 76.569861 days.
  Measured duration is provisional (notebook says 77.04 days); verify ROOT metadata.
- Existing fixes: 171 -> 19 unique backgrounds; five measured signal durations
  corrected for future validation. Working script now uses `SNparams["t"]` for all
  process times. Keep the intentional 10-second optimizer test.
- Extract ROOT `last_event_number` once per unique file into a second analysis CSV.
  Never replace generated totals with selected tree entries or overwrite source CSV.
- [activity_comparison.csv](activity_comparison.csv): 22 differences, five unchanged;
  notebook cell 3 defaults, PMT fit factors disabled. Activities await colleague review.
  Mappings: bb/bb2nu, bi_foil/214Bi, pa/234Pa, k/40K, tl/208Tl.
- Exclude missing extra backgrounds/neutrons and notebook-specific output code.
  Keep angle independently optimized and use existing package Feldman-Cousins methods.

## Configuration contract

- Paths resolve from repository root; SLURM CPUs set Julia threads. Empty sysimage
  disables it. Validate job time against optimizer limit plus overhead before submit.
- A named cut set contains fixed topology and variables with role `fixed` or
  `optimize`. Fixed variables use bounds; optimized ones also need step/initial ROI.
  Reject invalid bounds, off-grid edges, missing branches/files, and duplicate processes.
- Baseline preserves current topology/timing branch. OM quality is explicitly
  off until a phase-3 list exists. Approved energy-sum step is 0.1 instead of 0.2
  so bounds and initial ROI align; all seven current ND variables remain optimized.
- Each run saves resolved TOML, metadata snapshot, seed/version/job provenance,
  logs, and CSVs: sensitivity/exposure/FC limit/efficiency/counts/convergence;
  ROI edges with units; component/group backgrounds with MC counts and sqrt(sum(w^2));
  cumulative cutflow. Also save a readable summary and requested PDF/PNG plots.
- Sensitivity slices vary one ROI edge with others fixed at the optimum, not
  re-optimization. Data validation stays disabled until files exist; never label it passed.

Timing decision: use absolute values consistently in the new workflow.

## Boundaries

- New work stays in `scripts/`; existing files change only with explicit approval.
- No local long jobs, event-data loading, downloads, or notebook execution.
- Stop on bugs/problems and ask who should fix them.
- Approve physics choices and metadata changes; never silently skip missing inputs.
- Estimates exclude approval waits/cluster runs; unavailable validation is not a pass.
