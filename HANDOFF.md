# Handoff — flashier-starcat core (M0–M6) on branch starcat/m0

**Written:** 2026-11-18 · **Branch:** starcat/m0 (from main @ b03ba4c) · **Status:** M0–M6 built, tested and committed locally, not pushed. Clean-clone R CMD check: tests and examples OK, no new NOTE. main's rOpenSci pkgcheck step is parked (D020). CI on b03ba4c: R CMD check, pkgdown and test-coverage all green.

## Check first
Run `git log --oneline 8495e92..HEAD` and `git status --short` on `starcat/m0`. If either shows work this file does not mention, say so before acting.

## Goal
Implement the flashier-starcat design bundle (`flashier-starcat/`) core path. Consensus references across datasets come from `build_factor_reference()`, and native EBMF projection onto them from `project_reference()`. The user chose scope M0–M6 (D020) under `/tusoskill`.

## Next action
Ask the user to choose, then act on the answer:
1. push `starcat/m0` so CI (including pkgdown for the new `_pkgdown.yml` section) runs on it;
2. continue toward the tusoskill floor (ledger at 8/50; `validation/starcat-20261118/`) with R035 starCAT parity, the R100/R103/R105–R107 receipts and R053 plots, each with its own seed and grid asks;
3. or stop this round, which the user must say in their own words before `record-blocker --kind user`.
Also offer a bench-v1 rescore: the recovery cut (congruence ≥ 0.9), 30% effect sparsity and `max_factors = 6` were defaults the user did not approve.

## State
- Commits: bf7b5b7 (M0 snapshot and fixtures), 8495e92 (M1–M6), plus a follow-up commit with edge-case tests, air formatting and this file.
- New exports: `build_factor_reference`, `project_reference`, `reference_programs`, `reference_members`, `reference_diagnostics`; S3 print and summary methods. Files: `R/reference-build.R`, `R/reference-class.R`, `R/reference-projection.R`, `tests/testthat/test-reference.R`, `tests/testthat/test-reference-baseline.R`, `tests/testthat/helper-reference.R`.
- Tracker: 38 tasks DONE in `flashier-starcat/03_EXECUTION_TRACKER.{csv,md}`; records in `flashier-starcat/10_COMPLETION_RECORDS.md`.
- Tests: full suite 0 failures after 8495e92; the reference file has about 95 expectations.
- bench-v1 (R101): `validation/starcat-20261118/benchmark/bench-v1/f1_by_noise.png` and `summary.csv`. F1 is 0.85–0.97 at `min_similarity` 0.5–0.7. There were 0 duplicate-dataset clusters and 0 private-program leaks in 240 builds, with a median build of 0.5 s. The low-noise dip comes from source-fit rotation, not the reference layer. At SD 0.25, 0.70 of 3 truths were recovered in all 3 fits versus 2.55 by the reference, and the reference never fell below what the sources allow (0/80) — `validation/starcat-20261118/benchmark/bench-v1/source_vs_reference.png`.
- Assumed, not run: CI on starcat/m0; the starCAT parity script (R035). The starcat2025 parity test uses datasets with identical feature sets, so renormalization over shared features is not exercised.

## Locked decisions
- D019 (user, 2026-11-18): the spec's generic defaults are confirmed — L1/L2 profiles, congruence, sign flip only signed/signed, intersection features, required `min_similarity`, reciprocal FALSE, density ≥ 2/3, `min_datasets` 2. The starcat2025 constants are also confirmed.
- D021/D022 (user, 2026-11-18): starcat2025 uses starCAT's cross-cluster directed density > 0.666 and the unrenormalized mean of source-SD spectra.
- D023/D024 (agent, recorded): ties break by similarity, then radix member ID. The per-dataset best filter runs only in starcat2025 mode. `reference_id` hashes a permutation-invariant source digest.
- bench-v1 grid (user): seeds 1:20; noise 0.25/0.5/1/2; `min_similarity` 0.5/0.7/0.9.
- The earlier main-branch decisions still hold: see `git show b03ba4c:HANDOFF.md` (inference fixes D1–D8, `prior_policy`, removed arguments, provenance rules).

## Dead ends — do not redo
- `rcmdcheck` from the working tree copies `validation/` (about 10 minutes); check a clean clone in the scratchpad instead.
- Main-branch dead ends: see `git show b03ba4c:HANDOFF.md`.

## Read first
- `flashier-starcat/07_SOURCE_AUDIT_AND_DECISIONS.md` D018–D024.
- `flashier-starcat/10_COMPLETION_RECORDS.md`.
- `flashier-starcat/05_TEST_BENCHMARK_PLAN.md` for the open receipts.

## Open items
- rOpenSci path on main (parked): fresh-clone pkgcheck; user actions — JMLR citation check, `CODECOV_TOKEN`, upstream-issue decision, opening the submission issue.
- Open StarCAT tasks: R035, R053, R100, R101 (rescore), R103–R109, M7–M9, M11, M12 (vignettes, final release gate).
