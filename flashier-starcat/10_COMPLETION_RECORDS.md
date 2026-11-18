# Completion records (08_AGENT_HANDOFF §15)

Branch `starcat/m0`, built 2026-11-18 from `main` at b03ba4c. Scope M0–M6 (D020).

## M0 — R000–R003

- Files: `tests/testthat/helper-reference.R`, `tests/testthat/test-reference-baseline.R`, `tests/testthat/fixtures/m0-legacy-baseline.rds`.
- Public API impact: none.
- Scientific contract: `factor_replicability()`, `project_factor_activity()` (query and `fixed_reference`) and `project_factor_features()` outputs must equal the b03ba4c snapshot within 1e-6.
- Receipt: full suite at b03ba4c, 0 failures; snapshot generator in the local archive `validation/starcat-20261118/benchmark/m0-baseline-generate.R`.

## M1–M3 — R010–R016, R020–R023, R030–R034

- Files: `R/reference-build.R` (internal only).
- Scientific contract: comparison profiles never touch source views. Similarity is all-pairs, with a reason for every ineligible pair. Same-dataset pairs are diagnostic only. Sign equivalence holds only for signed/signed pairs. Directed candidates, accepted edges and merge decisions are stored.
- Tests (`test-reference.R`): source immutability, and order, permutation, rescaling and serialization invariance. Also identical, orthogonal, sign-reversed, Pearson-versus-congruence and low-overlap pairs; equality with `factor_replicability()` congruence; LSAP mocked to fail (R023); and the graph cases (pair, triangle, chain at and below 2/3, starCAT cross-density chain, duplicate dataset, ties, disconnected groups).
- Open: R035 (parity script against starCAT Python) is not done.

## M4–M5 — R040–R046, R050–R052

- Files: `R/reference-class.R`.
- Public API: `build_factor_reference()`, `reference_programs()`, `reference_members()`, `reference_diagnostics()`, plus print and summary methods.
- Scientific contract: one member per dataset, constraint-incompatible clusters excluded from programs, singletons kept as diagnostics, and an explicit reference feature basis. `reference_id` hashes scientific content only (D024).
- Provenance: five rows in `inst/provenance.csv`; CREDITS.md entry; NEWS.md.
- Open: R053 (plots) not done.

## M6 — R060–R065

- Files: `R/reference-projection.R`.
- Public API: `project_reference()`.
- Scientific contract: refuses references that are not projection-ready, and requires `feature_scale` to match the projection manifest. Needs a declared activity prior with matching support, keeps the programs fixed (checked inside `.native_projection()`), and returns L1 activity only for nonnegative priors.
- R065: no refactor was needed. `.native_projection()` is reused unchanged, and the M0 snapshot test proves the legacy outputs did not change.

## M10 — R102 (split and merged factors)

- Receipt: a deterministic unit test in `test-reference.R`, not an external benchmark. D2 splits P1 into two halves and D3 merges P1 and P2. No cluster holds two D2 factors. The merged D3 factor (congruence 0.707 with P1) stays a singleton at `min_similarity = 0.8`. P1's program counts 2 datasets, not 3.

## Not started

R035, R053, M7–M9, M10 receipts other than bench-v1 (R101), M11, M12.
