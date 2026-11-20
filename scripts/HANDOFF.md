# Resume here

Updated: 2026-11-17. Implementation is ready locally; cluster verification and
physics-input review are still pending. Do not start production yet.

## First action next session

Ask for the cluster smoke-test log/job ID. The user may already have submitted
the batch job, but neither submission nor its result has been confirmed here.
Inspect an existing result before submitting another job.

If not submitted, run these commands **on the cluster, from the repository root**:

```sh
git pull --ff-only
sbatch scripts/smoke_test.sh
```

Inspect `smoke-JOBID.log` in the repository root using the returned job ID:

```sh
squeue -j JOBID
tail -n 80 smoke-JOBID.log
```

Expected: **45 metadata checks + 52 workflow checks pass**, with no skipped ROOT
check, ending with `Smoke test completed successfully.` Missing/fractional ROOT
parameter errors are expected negative tests; inspect the final summaries.

The batch script requests 1 CPU, 8 GB, and 10 minutes (queue wait additional).
It loads `julia root`, sets a generic CPU target, disables package images, and
uses no custom sysimage. Dependencies must already be available on the cluster.
Tests create temporary parameter-only ROOT files and synthetic events, not
production event trees. Outputs are temporary apart from the batch log.

## Current evidence and limitations

- Local suites: **42 metadata + 52 workflow checks passed**. Local ROOT/C++
  integration was explicitly skipped because ROOT is unavailable.
- Configuration-only extraction preview selected **20 processes**: requested
  phase-3 signal plus 19 backgrounds. No real ROOT files were opened.
- Submission preview and shell syntax checks passed. The new batch script also
  passed mocked success/failure checks, including stopping after a failed suite.
  This does not establish that it has run successfully on the cluster.
- Earlier `srun` obtained an allocation but failed with
  `Communication connection failure` before Julia started. No tests ran.
  The user prefers the established `sbatch run.sh` style; use `sbatch`, not `srun`.
- The actual `scripts/mathis_analysis_info.csv` has **not** been generated here.
  Activities have **not** been changed. No production job/result is confirmed.

## Remaining sequence

1. Inspect/fix the cluster smoke test, stopping to ask who should fix any problem.
2. Review [activity_comparison.csv](activity_comparison.csv) with the colleague.
   There are 22 differences and five unchanged rows; these are candidates, not
   approved corrections. Optional notebook PMT fit multipliers were excluded.
3. Apply only agreed activities to a reviewed metadata copy, supplied through
   `--source`. Extract real generated counts into the separate analysis CSV on
   the cluster; preserve the original source CSV. See [README](README.md).
4. Run metadata/input preflight, then submit production only after the input
   review and checks pass. Measured-data validation stays disabled until ready.

The count extractor itself does not depend on activities, so extraction can
technically happen earlier. However, `--overwrite` regenerates from the chosen
source and replaces manual edits in the derived CSV: do not lose reviewed changes.

## Approved scientific decisions

- Sensitivity projects **2.86 years = 90192960 seconds**, under **phase-3
  conditions**, approximately **17.4746 kg yr** with package mass/enrichment.
  Do not combine/phase-weight detector conditions. Process exposure and the
  sensitivity normalization must use the same time.
- Data/MC validation compares counts/day: simulated counts / **1043.9 days**;
  measured counts / **76.569861 days**. The measured phase-3 time is provisional:
  the notebook has conflicting saved output (77.04 days). Confirm before final
  measured-data claims; do not silently change it.
- Keep all seven current ND variables optimized by default, including angle.
  Do **not** adopt the colleague's fixed 50-degree angle cut.
- New workflow uses absolute ND values consistently, including timing:
  lower-inclusive/upper-exclusive windows. Timing bounds `[0,5]`, initial
  `[0,2.5]`; approved energy-sum step **0.1 MeV**, replacing off-grid 0.2.
- Keep existing package Feldman-Cousins methods and intentional **10-second**
  optimizer test limit. This limit does not bound data loading or total job time.
  Exclude unavailable extra source components/neutrons. OM quality is disabled
  pending a real phase-3 list.

## Implementation map

- [job.example.toml](job.example.toml): SLURM, data, optimization, named cut sets,
  outputs, background groups, and optional measured-data comparison.
- [job_common.jl](job_common.jl), [submit.jl](submit.jl), [run.jl](run.jl):
  strict validation, configuration-only submission preview, metadata-only
  preflight, unique output directories, snapshots, hashes, logs, and run status.
- [analysis.jl](analysis.jl): package-backed optimization, full result/background
  tables with MC errors, cutflow, spectra, sensitivity slices, and optional
  identically cut data/MC comparison. Slices vary one edge, not re-optimized
  profiles. Comparison output is not claimed to be a validation pass.
- [prepare_metadata.jl](prepare_metadata.jl), [metadata.jl](metadata.jl),
  [read_generated_counts.C](read_generated_counts.C): copy exact positive
  `last_event_number` parameters without reading event trees, retaining original
  counts/activities and extraction provenance. No automatic `+1`; confirm the
  producer's count-versus-final-ID convention against generation records.
  File-size/mtime checks reject changed ROOT inputs; writes protect the source
  and existing output unless explicit overwrite is requested.
- [test_metadata.jl](test_metadata.jl), [test_workflow.jl](test_workflow.jl),
  [smoke_test.sh](smoke_test.sh): focused tests and bounded batch entry point.

Package details relevant to future changes:

- The package's ND ROI code uses absolute values, whereas its `filter_data`
  helper uses signed strict inequalities; the new workflow handles selection
  consistently instead of mixing these APIs.
- Sensitivity normalization independently uses `SNparams["tYear"]`; the workflow
  copies the physics dictionary and synchronizes seconds/years.
- Package floating-point grid flooring and ROI histogram-cut behavior required
  explicit grid construction and direct selected-event spectra in this workflow.

## Existing source changes and Git state

Inputs were [the working ND example](../examples/run_nd_sensitivity_on_mathis_data.jl),
[the batch example](../examples/run.sh), and
[the colleague's notebook](<../examples/ex8_sensibilite_finale(2).ipynb>).
Approved existing-file changes:

- Working ND example: deduplicate selected-phase backgrounds (171 -> 19);
  normalize all processes with `SNparams["t"]`. Its timing semantics were not changed.
- [Original metadata CSV](../data/mathis_data_info_final.csv): five signal
  measured-phase durations corrected. Activities/generated totals remain unchanged.

Last pushed code baseline: **2beabaf**, `Use sbatch for bounded cluster smoke tests`.
Previous **c667e08** contains the main workflow. A rejected commit accidentally
included a 186 MB partial ROOT download; it was amended before pushing, without
force-push. Local data was kept. Root ignore rules now cover downloaded data,
ROOT files/partial downloads, and logs; metadata CSVs/scripts remain tracked.

This handoff and its link in the plan are local documentation additions, not
committed/pushed as part of this wrap-up. Do not assume the cluster has them.

## Working agreement

- Keep changes in `scripts/` unless explicit approval covers existing files.
- No local long jobs, real event-data loading, downloads, or notebook execution.
- **On any bug/problem, stop and ask whether the assistant or user should fix it.**
- Keep responses short and actionable; approve substantial new steps individually.
  Ask questions through the structured question tool.
- Use an explicit Julia project, never a global environment. Locally:
  `--project=/Users/maros.petro/Work/Phd/SN_Sensitivity_training`.
- Preserve unrelated/user changes. Do not install dependencies speculatively,
  silently skip missing inputs, adopt candidate activities, or call unrun checks passed.

Suggested opening prompt for the next session:

> Read scripts/HANDOFF.md and scripts/plan.md. Resume with the sbatch smoke-test
> result, not production. Keep responses short and stop to ask on problems.
