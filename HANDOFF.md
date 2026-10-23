# Handoff: GBCD and ImmGenT membership filtering

Date: 2026-10-23. Branch: main.

## Delivered implementation

- `fit_gbcd_programs()` calls optional gbcd 0.2-10 at
  `670e3afb50236303ddf99db33a3394abd4a711bd`. It returns representation,
  signatures, diagnostics and optional native output. It checks memory branches,
  preserves Baseline as a factor and records unknown activity support.
  FLASH indices/PVE/posterior samplers remain unavailable.
- `filter_cells_by_total_membership()` implements the exact ImmGenT iterative
  rule: cutoff 10, 12 passes, maximum scaling, paired effects compensation.
  It returns retained/removed IDs and pass history without modifying the fit.
- Public API: 23 exports, 34 retained internal implementations, 57 source records.
- Workshop model guidance and GoM DE distinctions are linked with attribution
  in `vignettes/advanced-use.Rmd`. No unlicensed workshop source was copied.

## Evidence

Archive: `validation/gbcd-integration-20261023/`.
- Main package check: `check-final.log`, Status OK; full installed suite
  1,575 passes, zero failures/skips, 11 known projection warnings.
- Final metadata/control guards: `gbcd-qualified-tests.rds`, 63 passes,
  zero failures/warnings/skips. Membership tests: 39 passes.
- Final tarball check: `check-artifact.log`, Status OK; tests not repeated.
  It checks examples and rebuilds all four vignettes.
- `artifact-receipt.json` proves root/source-tarball equality and preservation
  of all original provenance rows. Scientific archives stay outside builds.
- Dense native parity: permanent tests; sparse native parity: permanent tests
  plus `lowrank-parity.rds`. Iteration-limited parity does not prove convergence.
- Pre/post and source/installed smoke are identical. Site: 23 references,
  four tutorials, attribution and local links pass.

## Next action

Commit/push complete: `70a30d6b5c10cbabf44ad8d7355dc5c14c16f3c1`.
Live origin/main verified at the same SHA. Inspect hosted CI for that SHA.
Do not re-ask approval for these features.

## Separate pending work

GoM DE itself is not implemented. fastTopics accepts count-topic fits or
row-sum-one nonnegative memberships, with counts; arbitrary signed EBMF
factors do not qualify. Existing distinctive ranking is descriptive.
Plot gallery has a committed plan, not an implemented page.
No scientific calibration study, scheduler job, new release or tag was run.
Preserve `reports/core-review-20261022/` and all existing archives.
