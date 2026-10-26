# Handoff — docs review, planted-truth benchmark and robustness article

**Date:** 2026-10-26 · **Branch:** main · **Status:** 617f6d1 (docs reorder + roxygen markdown) pushed, CI pending; planted-truth/robustness batch committed locally, not pushed

## Goal
Keep flashier.utils scientifically honest and shippable, working down
`docs/plans/2026-10-25-ease-triage.md` in rank order. Ranks 1 and 2 are done.

## Next action
Push the planted-truth/robustness commit, then `gh run list --limit 4`: both it and 617f6d1 need green R CMD check and pkgdown.

## Previous next action
Run `gh run list --limit 4`. The API cleanup push must pass R CMD check (5 jobs) and
pkgdown. Earlier runs 40357365764 (init_fn fix) and 40360265143 (test_gom_de) must also
be green; full-engines installs fastTopics from the pinned Remote. If green,
ask whether to start triage rank 3, the saved-fit comparison API (A1–A2 in
`docs/plans/2026-10-25-complete-package.md`).

## State
- Review reports are committed (user asked, 2026-10-26): core-review-20261022/,
  ebmf-cnmf-adversarial-review-20261025.md, package-review-20261026.md,
  starcat-status-20261025.md and api-review-20261026/. Probe harness and check log stay
  local in validation/api-review-20261026/.
- Tests: `scripts/run-r.sh -e 'devtools::load_all(); testthat::test_dir("tests/testthat")'`
  gives 185 blocks, 0 failures, 0 skips, 11 known projection maxiter warnings, ~4 min.
- Last green CI: bc920ec (R CMD check, 5 jobs, plus pkgdown). Ubuntu dependency setup takes ~40 min.
- Runtime: R4_51 via scripts/run-r.sh. roxygen2 8.1.0 installed, package records 8.1.1; document() still works.
- Local-only (git-excluded): development/, docs/plans/, validation/ (receipts).

## Shipped this thread (one line each)
- 617f6d1 docs site reordered value-first; roxygen markdown enabled (DESCRIPTION `Roxygen:`), all Rd regenerated; stale tracked docs/plans file removed.
- Planted-truth benchmark: archive `validation/ebmf-cnmf-planted-20261026/` (PREREG, JOURNAL, simulate/fit-cnmf/fit-ebmf/evaluate, receipt, SHA256). Derived CSVs in vignettes/articles/data/ebmf-cnmf-planted-*.csv; new section in ebmf-cnmf-comparison.Rmd. Deviations logged: auc_posthoc column added after first evaluation (PREREG rule unchanged); collinear cor 0.89-0.91 vs 0.85 target. All cNMF runs hit the 300-iteration cap (pilot config).
- New website article vignettes/articles/robustness-checks.Rmd; allowlisted in check-site.py + test fixture; navbar entry.
- Review fixes 13c8c16..e52fec5: var_type inheritance, fixed_reference prior policy, stability subsampling, docs, styler.
- 0097ab8 fixed_reference tests for intercept fits and factor subsets; b7f153a metadata extra rows dropped and recorded.
- 2884978 ledgers reconciled (GBCD tracker 150 done, 13 partial); f59a2d0 custom greedy init_fn fits are callback-required.
- Profiling: validation/profiling-20261026/DECISION.md, Rcpp not justified (factor_intervals time is the native sampler).
- ab8ce82 `test_gom_de()` exported fastTopics adapter (pin 5cf9a448); 2cabb52 Rd drift from styler regenerated.
- API cleanup: 12 single-value/reserved args removed from 7 exports, 4 trailing reorders,
  S/var_type/one-sample guards, factor_intervals side default; 322-probe edge matrix.

## Locked decisions
- `prior_policy` defaults to "query"; "fixed_reference" is opt-in (user, 2026-10-26).
- `.match_ids` defaults to `extra = "drop"` with exclusions recorded; sample-keyed holdout folds pass `extra = "fail"` (user, 2026-10-26).
- `test_gom_de` is exported, takes an explicit membership matrix only, never normalizes EBMF scores (user, 2026-10-26).
- basis_id keeps hashing engine versions; only the comment changed (2026-10-26).
- Push once per batch; the last pushed commit must not carry `[skip ci]` when R code changed — it skips CI for the whole push, and R CMD check has no workflow_dispatch.
- Core API migration shipped; do not re-ask (user, 2026-10-23).
- Single-value and reserved arguments are removed from exports; breaking signature changes
  approved by the user (2026-10-26). Do not re-add alignment/bootstrap/refit_inputs/control
  (factor_replicability), resample/alignment (factor_stability), freeze_noise/workers
  (factor_holdout), factor_uncertainty (test_factors), uncertainty (backproject_contrast)
  or engine (plots). Selections stay character IDs only.
- Preserve scientific contracts, original authors, inst/NOTICE, provenance, CREDITS (AGENTS.md).
- GBCD: optional exact pin 670e3afb…; Baseline stays a factor; activity support unknown; copy no algorithms.
- Keep LM/ash statistics separate from EB-SNMF effects; invent no FLASH uncertainty/PVE.
- ImmGenT membership filter defaults: threshold 10, 12 passes; paired rescaling preserves reconstruction.
- No staged comparison API promoted to R/ until contracts are qualified (2026-10-25).

## Findings to keep
- Fixed priors alone do not remove query-set dependence: only fixed priors plus sample-wise
  var_type do (same cell, two query sets: 5e-7 vs 4e-4 to 9e-4).
- fastTopics' matrix path checks row sums one-sided and only warns; test_gom_de checks both sides.

## Dead ends — do not redo
- Pages cleanup in tracked docs/ fails; build under RUNNER_TEMP/pkgdown-site.
- Baseline/GEP labels do not establish GBCD support lineage.
- Descriptive coefficient margins are not conditional GoM DE lfsr.
- Prototype signed reorientation changes cone coverage; fix sign/pole contracts first.
- CE track: reports/CE_G2_LOSS_MEMO.md says no rerun.

## Read first
- `reports/api-review-20261026/REPORT.md:1` — argument cleanup findings and migration.
- `docs/plans/2026-10-25-complete-package.md:1` — A1–A2 contract for rank 3 (staged candidate under validation/implementation-all-20261025/candidate/).
- `docs/plans/2026-10-25-ease-triage.md:1` — remaining ranks 3–5 and follow-ups.
- `R/gom-de.R:49` — newest adapter; pattern for optional engines with `R/gbcd.R`.

## Open items
- GB-005b display order lacks canonical MDS sign and declared group order — fix on request.
- GO_MEMO_SCI004.md and the audited plan define SCI-004a/b differently — reconcile before any release GO.
- GBCD tracker has 13 partial rows (small missing tests) — flashier-gbcd/TRACKER.md.
- Possible follow-up: fix residual variance to the reference (reference tau as S).
- StarCAT (SC-*), SCI-*, VAL-* backlogs unchanged; need data and Slurm.
