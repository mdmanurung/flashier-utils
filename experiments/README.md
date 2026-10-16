# Experiments index

Top-level index tying each benchmark/study folder to the scripts that ran it
and the outputs it produced. Package source (`R/`, `tests/`, `man/`) is
separate; nothing indexed here runs during `R CMD check`. Conventions for
adding new studies are at the bottom.

## Study index

### `benchmarks/` — backend baselines and sync smoke tests

- **`fu_baseline.R`** — testthat baseline tabulation across dependency
  versions (e.g. flashier 1.0.60 vs 1.0.7): loads the package via
  `pkgload::load_all()` and tabulates passed/failed/error/skipped per test
  file. *Outputs*: per-file baseline tables written to the session workspace;
  not committed.
- **`sync3_smoke.R`** — SYNC-3 release-gate smoke test: refits one Stephenson
  cell type (`CD4_naive`, K30 masked pseudobulk) through the migrated
  flashier.utils backend and checks parity with the committed masked fit
  (n factors, ELBO, `L_pm`/`F_pm`/score diffs, fallback flag). *Outputs*:
  console parity receipt; recorded in the cellprograms repo
  (`reports/SYNC_ROADMAP.md`, cross-repo) and the sync commit history. No
  figures.

No committed figures or reports yet in this repo; when a study produces
them, commit them under the study folder (`outputs/` or `figs/`) and add a
row here.

## Conventions for new experiments

1. **One study = one folder**: `benchmarks/<study_slug>/`, snake_case slug,
   self-contained (scripts, README, small outputs).
2. **Script naming**:
   - Numbered pipeline scripts in execution order: `01_<step>.R`, `02_<step>.R`
     (cluster wrappers as `01_<step>.sbatch` next to them).
   - Ad-hoc drivers: verb-prefixed snake_case — `run_<thing>.R`,
     `summarize_<thing>.R`, `prof_<thing>.R`, `smoke_<thing>.R`.
   - When a cleaned/numbered script exists, the exact executed version goes in
     `scripts/as-run/` (verbatim copy, path warts included).
3. **README required** in every study folder, with: (a) the question, (b) a
   script table, (c) inputs and where caches come from, (d) outputs and where
   they land, (e) as-run environment caveats (absolute paths, package
   versions, cluster notes).
4. **Outputs**:
   - Commit small figures/tables under the study folder (`outputs/` or
     `figs/`), or as a vignette article when package-facing.
   - Never commit large intermediates (`.rds`/`.hdf5` caches); record their
     path and the script that regenerates them.
   - Mark session-only outputs as such in the README and here.
5. **Load the package, never edit it from experiments**: use
   `pkgload::load_all()` or `library()`; analysis code never goes into `R/`.
6. **Update this index in the same commit** that adds or changes a study
   folder or its outputs.
