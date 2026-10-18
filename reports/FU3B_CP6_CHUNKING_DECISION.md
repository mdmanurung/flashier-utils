# FU-3B / CP-6 — Chunked projection feasibility decision

**Date:** 2026-10-19
**Scope:** `flashier.utils::project_factor_activity()` /
`project_factor_features()` (FU-3B) and
`cellprograms::project_programs_to_cells()` (CP-6).
**Question:** can query cells be processed in row chunks with results
identical to a single joint projection?

## Separability audit

Row-chunk equivalence requires every fitted quantity for row *i* to depend
only on row *i*'s data. Auditing the projection path:

| Component | Row-separable? | Evidence |
|---|---|---|
| Preprocessing replay (cellprograms: per-gene OLS residualization with fit-time coefficients, then fit-time center/scale) | Yes | `R/cell-projection.R`: replay uses only fit-time `residualize_coefficients`, `center`, `scale`; no cross-row statistics |
| Fixed effects / offsets | Yes | `B` and offsets are frozen from the reference fit (`R/projection.R`) |
| Activity prior hyperparameters | **No** | `.projection_priors()` resolves prior *functions* (e.g. `ebnm_point_laplace`); the native backfit re-estimates the prior `g` by empirical Bayes pooling over **all query rows in the fit**. A chunk gets its own `g` → different shrinkage → different activities |
| Residual variance, `var_type = 0` (constant) | **No** | single scalar pooled over all elements |
| Residual variance, `var_type = 1` (by row) | Yes | per-row estimation |
| Residual variance, `var_type = 2` (by column) | **No** | pooled over rows |
| Convergence / iteration count | Approximately | stopping is global to each fit; chunk fits can terminate in different local states |

## Decision

**NO-GO for exact chunking under the current API.** Exact chunked projection
would require frozen *fitted* priors (fixed `g`, not re-estimated) together
with frozen or per-row variance. flashier does not expose fixed-`g`
projection through this path; that is an upstream feature, out of scope for
this milestone.

**Approximate chunking is not implemented.** No registered workflow has
demonstrated a memory need: oversize queries are refused pre-flight by
`.projection_budget()` against `dense_budget_bytes` (default 256 MB) with an
explicit error naming the remedy. If a future need is demonstrated, a
chunked mode must (i) be opt-in and labelled APPROXIMATE, (ii) record
chunking in `analysis_metadata`, (iii) ship equivalence-bound tests that
report `max |A_chunked − A_joint|` rather than asserting zero.

## Memory guarantee for v1 (CP-6)

Bounded memory is delivered as **refusal, not silent chunking**: peak dense
allocations are estimated and checked before any allocation; queries that
exceed the declared budget fail fast. Per-cell-type loops in
`project_programs_to_cells()` bound peak usage to one cell type at a time.

## Recorded limitation (feeds REL-1 register)

Projected activities are empirical-Bayes estimates pooled over the query
set: **the same cell projected in two different query sets can receive
different activities** (different estimated prior, and — at
`var_type = 0/2` — different variance). This holds with or without chunking
and is now stated in the function documentation.
