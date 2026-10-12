# Revised Development Milestones

The plan is split into a **core reference/projection track**, a **conditional residual-program track**, and an **experimental prior-transport track**.

Core release does not wait for experimental prior transport.

---

## M0 — Baseline, source snapshot and regression guard

**Goal:** freeze current behavior before introducing reference abstractions.

### Tasks

- [ ] Record repository HEAD, `DESCRIPTION`, `NAMESPACE`, existing exported functions and dependency pins.
- [ ] Freeze regression tests for `factor_replicability()`.
- [ ] Freeze regression tests for `project_factor_activity()` and `project_factor_features()`.
- [ ] Freeze `.native_projection()` behavior through public callers.
- [ ] Add deterministic signed, semi-nonnegative and non-negative fixtures.
- [ ] Record current package check receipt in the external validation archive.
- [ ] Add planned provenance rows as `planned` only in development notes, not package registry yet.

### Acceptance

- Existing tests pass unchanged.
- Added regression tests capture current projection/alignment behavior.
- No public behavior changes.

### Depends on

None.

---

## M1 — Comparison-profile contract

**Goal:** derive comparable program profiles without mutating source inferential coordinates.

### Tasks

- [ ] Implement source collection/validation.
- [ ] Implement member eligibility table with explicit exclusions.
- [ ] Require zero offsets for schema-v1 reference building.
- [ ] Implement generic non-negative L1 profiles.
- [ ] Implement generic signed/semi-nonnegative L2 profiles.
- [ ] Implement constraint-aware sign-orientation rules.
- [ ] Implement pairwise feature-set policy.
- [ ] Implement final reference-feature policy.
- [ ] Implement `starcat2025` source-statistics validation and profile transform.
- [ ] Add invariance tests for source factor rescaling and factor permutation.

### Acceptance

- Source `factor_representation` objects serialize identically before and after profile construction.
- Equivalent rescalings of the same source factor produce identical comparison profiles.
- Semi-nonnegative factors are never sign-flipped.
- Invalid STARCAT-style inputs fail with actionable errors.

### Depends on

M0.

---

## M2 — All-pairs cross-dataset similarity

**Goal:** compute the raw material needed for a many-dataset consensus graph.

### Tasks

- [ ] Extract reusable feature-matching and similarity primitives from `alignment.R` without changing `factor_replicability()` behavior.
- [ ] Implement internal all-pairs similarity table.
- [ ] Support congruence, Pearson and cosine.
- [ ] Implement sign-equivalence rules from source constraints.
- [ ] Record pairwise feature coverage and exclusions.
- [ ] Exclude same-dataset pairs from graph eligibility.
- [ ] Confirm no LSAP assignment is used in consensus graph construction.

### Acceptance

- Pairwise congruence agrees with existing `.align_programs()` before assignment.
- Similarity output is invariant to dataset/factor input ordering.
- Every ineligible pair has a reason.

### Depends on

M1.

---

## M3 — Program graph and deterministic clustering

**Goal:** implement auditable reciprocal-neighbor graph construction and density-constrained clustering.

### Tasks

- [ ] Build directed candidate-neighbor table.
- [ ] Build accepted reciprocal undirected edge table.
- [ ] Implement similarity threshold and optional rank restriction.
- [ ] Implement STARCAT-2025 defaults.
- [ ] Implement deterministic descending-edge cluster merging.
- [ ] Enforce one member per dataset per cluster.
- [ ] Implement explicit accepted-edge density calculation.
- [ ] Record rejected merge decisions and reasons.
- [ ] Add hand-verifiable graph fixtures.
- [ ] Add STARCAT graph-parity benchmark script outside package builds.

### Acceptance

- Cluster output is invariant to input order.
- Duplicate-dataset membership is impossible.
- Graph and merge decisions are fully reconstructable from stored tables/settings.
- STARCAT-mode toy/parity fixtures match the predeclared algorithm contract.

### Depends on

M2.

---

## M4 — Consensus estimation and `factor_reference`

**Goal:** create a stable, serializable consensus reference.

### Tasks

- [ ] Recompute member profiles on the final reference feature universe.
- [ ] Align member signs only when constraints permit.
- [ ] Implement arithmetic-mean consensus.
- [ ] Renormalize consensus profiles after averaging.
- [ ] Reject constraint-incompatible clusters from projection-ready programs while preserving diagnostics.
- [ ] Compute cluster support/similarity/heterogeneity diagnostics.
- [ ] Implement `factor_reference` constructor and validator.
- [ ] Implement deterministic `reference_id`.
- [ ] Implement `print.factor_reference()` and `summary.factor_reference()`.
- [ ] Implement `reference_programs()`, `reference_members()` and `reference_diagnostics()`.
- [ ] Add serialization round-trip tests.

### Acceptance

- Consensus means match hand calculations.
- All final programs use one explicit common feature basis.
- `reference_id` is unchanged by input ordering or serialization.
- Invalid schema objects fail closed.

### Depends on

M3.

---

## M5 — Top-level reference builder and diagnostics

**Goal:** expose one coherent user workflow.

### Tasks

- [ ] Implement `build_factor_reference()` orchestration.
- [ ] Implement generic defaults.
- [ ] Implement STARCAT-2025 mode validation/defaults.
- [ ] Implement projection-ready versus discovery-only references; schema-v1 native projection requires one compatible activity-support contract across selected programs.
- [ ] Store compact source provenance rather than full fits by default.
- [ ] Add optional diagnostic plots only where needed.
- [ ] Add end-to-end synthetic example.
- [ ] Add source-credit/provenance hooks.

### Acceptance

- One call builds a valid reference from multiple source fits/views.
- Build is deterministic without RNG.
- Discovery-only reference is valid but refuses projection.
- Projection-ready reference has an explicit query transform contract.

### Depends on

M4.

---

## M6 — Native EBMF reference projection

**Goal:** project `factor_reference` programs into a novel dataset using existing native `flashier` machinery.

### Tasks

- [ ] Implement `reference_projection` return class.
- [ ] Implement `project_reference()` using `.native_projection()` or a shared refactored equivalent.
- [ ] Add exact/intersection feature matching.
- [ ] Add one-to-one gene mapping.
- [ ] Preserve dense-memory budget checks.
- [ ] Require explicit or reference-declared projection prior family compatible with the reference activity-support contract.
- [ ] Keep query-side prior/noise recalibration explicit in metadata.
- [ ] Add optional L1 normalized activity for non-negative activity support.
- [ ] Refactor shared legacy projection internals only after regression coverage exists.
- [ ] Preserve legacy public APIs.

### Acceptance

- Fixed reference programs remain unchanged within current projection tolerance.
- Raw activity is always retained.
- Feature coverage/missing/excluded sets are explicit.
- Legacy projection regression tests pass.
- Query preprocessing incompatibility fails before fitting.

### Depends on

M5.

---

## M7 — Residual-program native feasibility gate

**Goal:** determine whether a clean public residual-program extension can be built using exported `flashier` APIs.

### Tasks

- [ ] Probe fixed feature programs through `flash_greedy()`.
- [ ] Probe fixed feature programs through subsequent `flash_backfit()`.
- [ ] Verify reference activities can update while reference feature effects remain fixed.
- [ ] Determine safe null-check behavior for mixed fixed/reference and newly added factors.
- [ ] Verify factor IDs/order can be mapped deterministically after greedy additions/removals.
- [ ] Record a go/no-go decision.

### Acceptance

All mandatory probes in `06_FLASHIER_FEASIBILITY_PROBES.md` pass, or the feature is explicitly deferred.

### Depends on

M6.

---

## M8 — Residual-program discovery (conditional)

**Goal:** expose query-residual factors after fitting the fixed reference.

### Tasks

- [ ] Implement `discover_residual_programs()` only if M7 is GO.
- [ ] Add residual factors with native greedy fitting.
- [ ] Backfit without altering reference feature programs.
- [ ] Report known-only and full reconstruction.
- [ ] Report incremental SSE/fit-R2 improvement.
- [ ] Align residual programs back to the reference.
- [ ] Add descriptive residual evidence statuses.
- [ ] Add optional stability hooks.

### Acceptance

- Reference feature programs remain fixed.
- Simulated omitted programs can be recovered as residual structure.
- Complete-reference simulations do not systematically produce large stable residual programs.
- Output never labels a residual factor as biologically novel without external evidence.

### Depends on

M7 GO.

---

## M9 — Stability/support integration

**Goal:** combine within-source stability evidence with cross-dataset reproducibility without conflating the two.

### Tasks

- [ ] Accept existing `factor_stability()` summaries for member factors.
- [ ] Add member stability fields to reference diagnostics.
- [ ] Add consensus-level descriptive stability summaries.
- [ ] Add source-dataset leave-one-out reference robustness analysis in benchmark scripts.
- [ ] Keep stability weighting off by default.

### Acceptance

- Stability data are optional.
- Missing stability does not block reference construction.
- Cross-dataset support and within-dataset recovery remain separately reported.

### Depends on

M5.

---

## M10 — Scientific validation and benchmarks

**Goal:** establish scientific behavior and regression receipts before release.

### Tasks

- [ ] Exact synthetic invariance tests.
- [ ] Shared versus dataset-specific program simulations.
- [ ] Split/merged factor stress test.
- [ ] Feature-overlap stress test.
- [ ] Source-scale/preprocessing mismatch stress test.
- [ ] STARCAT graph parity test on identical similarity/profile inputs.
- [ ] Query missing/extra reference-program simulation.
- [ ] Small-query projection benchmark.
- [ ] Signed/semi-nonnegative benchmark.
- [ ] Leave-one-dataset-out real-data benchmark.
- [ ] Existing pancreas multi-technology benchmark where appropriate.
- [ ] Runtime/memory benchmark.
- [ ] Residual-program benchmarks if M8 implemented.

### Acceptance

See `05_TEST_BENCHMARK_PLAN.md` release gates.

### Depends on

M5 and M6; residual-specific items depend on M8.

---

## M11 — Experimental prior transport (deferred research track)

**Goal:** investigate, not assume, whether reusable activity priors improve transport.

### Tasks

- [ ] Probe extraction of fitted priors with `flash_fit_get_g()`.
- [ ] Probe replay with `flash_ebnm(g_init=..., fix_g=TRUE)`.
- [ ] Determine whether per-program priors can be expressed cleanly for multi-factor projection.
- [ ] Define how source activities must be rescaled to the final consensus basis.
- [ ] Define dataset-balanced prior estimation if pooling is attempted.
- [ ] Benchmark fixed versus query-re-estimated priors under composition shift.
- [ ] Only then design a public API, if justified.

### Acceptance

No public fixed-reference-prior API is added without a positive benchmark and an explicit transportability rationale.

### Depends on

M6. Does not block core release.

---

## M12 — Documentation, provenance and release

**Goal:** make the core workflow reproducible and agent/user accessible.

### Tasks

- [ ] Add consensus-reference vignette.
- [ ] Add projection vignette.
- [ ] Add residual-program vignette only if M8 exists.
- [ ] Update `inst/provenance.csv`.
- [ ] Update `CREDITS.md` and `vignettes/credits.Rmd` consistently.
- [ ] Update `_pkgdown.yml`, `NAMESPACE`, manuals and `NEWS.md`.
- [ ] Ensure package checks do not download benchmark datasets.
- [ ] Store benchmark receipts outside package builds per repository instructions.
- [ ] Run full package check on supported CI targets.

### Acceptance

- All exported functions documented.
- Credits/provenance are internally consistent.
- No required remote downloads during package check.
- Core benchmark report is frozen and linked from development/release notes.
- Full checks pass.

### Depends on

M10 for core release. M8 documentation is conditional.

---

# Release tracks

## Core reference/projection release

```text
M0 -> M1 -> M2 -> M3 -> M4 -> M5 -> M6 -> M10 -> M12
```

## Residual-program extension

```text
M6 -> M7 -> [GO] -> M8 -> residual M10 -> M12
```

## Experimental prior transport

```text
M6 -> M11
```

This track is explicitly non-blocking.
