# Handoff — package review items implemented, pushed and CI-verified

**Date:** 2026-10-26 · **Branch:** main, in sync with origin · **Status:** suite green; CI green at bc920ec; triage rank 1 done

## Goal
Keep flashier.utils scientifically honest and shippable. This thread implemented every
item of `reports/package-review-20261026.md`.

## Next action
Check CI for the init_fn fix push (`gh run list --limit 2`). Then ask the user whether
to start triage rank 2, the GoM DE adapter (C1 in `docs/plans/2026-10-25-complete-package.md`).

## State
- Commits: fix (projection var_type inheritance, matrix S check, tolerant fixed check,
  metadata extra rows dropped, intervals guard, runtime comment); feat fixed_reference
  prior policy; feat stability subsampling; docs (Rd, LIMITATIONS L1/L5); NEWS; styler.
- Pushed 13c8c16..HEAD to origin/main. CI run 40333350717 (R CMD check, 5 jobs) and
  40333350725 (pkgdown) passed on 2026-10-26; Ubuntu dependency setup took ~40 min.
- A pushed head commit tagged `[skip ci]` skips CI for the whole push; R CMD check has
  no workflow_dispatch, so an untagged commit is needed to trigger it.
- Follow-up (grilled and approved 2026-10-26): fixed_reference tests for intercept fits
  and factor subsets (no bug found); `.match_ids` defaults to `extra = "drop"` with
  exclusions recorded everywhere; sample-keyed holdout folds pass `extra = "fail"`.
- Tests: `scripts/run-r.sh -e 'devtools::load_all(); testthat::test_dir("tests/testthat")'`
  gives 175 blocks, 0 failures, 11 known projection maxiter warnings, ~4 min.
- Untracked, not mine to commit: reports/core-review-20261022/,
  reports/ebmf-cnmf-adversarial-review-20261025.md, reports/starcat-status-20261025.md,
  reports/package-review-20261026.md.
- Runtime: R4_51 via scripts/run-r.sh. Installed roxygen2 is 8.1.0, package records 8.1.1;
  `devtools::document()` still writes correct Rd.

- Triage rank 1 (2026-10-26): ledgers reconciled in 2884978 (GBCD tracker 150 done,
  13 partial, 1 open; memo and STARCAT notes). Local-only edits: development/ ledger,
  docs/plans progress lines. Profiling receipt validation/profiling-20261026/DECISION.md:
  all utility workloads under 1 s at 20k x 1000 x K20; factor_intervals 54 s is the
  native sampler, so no Rcpp.
- GBCD tracker row 33 fixed: caller greedy init_fn fits record refit_policy
  "callback_required" (R/fit.R), so stability refuses automatic replay.
- Known gaps found: GB-005b lacks canonical MDS sign and declared group order;
  GO_MEMO_SCI004.md and the audited plan define SCI-004a/b differently.

## Locked decisions
- `prior_policy` defaults to "query"; "fixed_reference" is opt-in (2026-10-26 plan).
- basis_id keeps hashing engine versions; only the comment was fixed.
- `.metadata_activity()` drops extra metadata rows; no new public argument.
- Core API migration approved and shipped; do not re-ask (2026-10-23).
- Preserve scientific contracts, original authors, inst/NOTICE, provenance (AGENTS.md).
- GBCD: optional exact pin 670e3afb…; Baseline stays a factor; activity support unknown; copy no algorithms.
- Keep LM/ash statistics separate from EB-SNMF effects; invent no FLASH uncertainty/PVE.
- ImmGenT membership filter defaults: threshold 10, 12 passes; paired rescaling preserves reconstruction.
- No staged comparison API promoted to R/ until contracts are qualified (2026-10-25).

## Findings to keep
- Fixed priors alone do not remove query-set dependence. Probe: same cell, two query
  sets, var_type 0: query 4e-4, fixed 9e-4; var_type 1: query 6e-4, fixed 5e-7.
  Only fixed priors plus sample-wise variance are query-independent.
- Mismatched ebnm family under fixed_reference gives a typed error before flashier.

## Dead ends — do not redo
- Pages cleanup in tracked docs/ fails; build under RUNNER_TEMP/pkgdown-site.
- Baseline/GEP labels do not establish GBCD support lineage.
- Descriptive coefficient margins are not conditional GoM DE lfsr.
- Prototype signed reorientation changes cone coverage; fix sign/pole contracts first.
- CE track: reports/CE_G2_LOSS_MEMO.md says no rerun.

## Read first
- `R/holdout.R:11` — `.reference_activity_g()` and the per-factor fixed-g loop.
- `R/projection.R` `project_factor_activity()` — prior_policy wiring and docs.
- `docs/plans/2026-10-25-ease-triage.md:1` — ranked remaining work (GoM DE, saved-fit comparison, B1–B3).
- `reports/starcat-status-20261025.md:1` — StarCAT SC-001..SC-007 still unbuilt.

## Open items
- Possible follow-up: fix residual variance to the reference (pass reference tau as S).
- GoM DE, StarCAT (SC-*), SCI-*, VAL-* backlogs unchanged; Rcpp profiling pending.
