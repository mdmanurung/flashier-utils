# Architecture decisions

## ADR-001 — Public upstream APIs only

Numerical work stays in `flashier`, `ebnm`, and `gbcd`.

## ADR-002 — Semantic activity/effects orientation

Never expose raw native row/column semantics as the user contract.

## ADR-003 — Projection is a new native fit

A projection result must carry the native indices, PVE, multipliers and scaling
of the **projection fit**, not those of the source fit.

Source basis information is retained separately.

## ADR-004 — Projection fits are not automatically refittable

Projection is a specialized fixed-side procedure. Stability requires an explicit
callback until that procedure is serialized as a first-class refit recipe.

## ADR-005 — Edited fits require explicit refit callbacks

Post-hoc backfit, removal, fixing and unfixing change the reference procedure in
ways not captured by the original standard configuration.

Pure reordering is the exception.

## ADR-006 — Do not serialize arbitrary EBNM closures

Unrecognized callbacks remain explicit-callback workflows.

## ADR-007 — Semi-NMF is a thin wrapper

The numerical algorithm already exists in FLASH.

## ADR-008 — FLASH intercept is an offset

The native fixed all-ones intercept is removed from semantic factor axes and
becomes `factor_representation$offset`.

## ADR-009 — Offset is derived from raw native means

Never derive an offset from a normalized LDF intercept column.

## ADR-010 — Native indices are authoritative

Once any native factor is hidden semantically, every posterior/PVE/sampler
operation must use `manifest$native_indices`.

## ADR-011 — Intercept stays native factor 1 in v1

Adaptive factors may reorder after it; the intercept itself does not move.

## ADR-012 — Feature projection must estimate new offsets

If a source model includes a feature baseline, projecting to new features
without estimating their baselines is an incomplete model.

## ADR-013 — Initialized fitting requires both sides

No hidden least-squares completion of a missing initialization side.

## ADR-014 — Initialized values must respect known support

A nonnegative prior and negative fixed/initial values is rejected early.

## ADR-015 — Initialized factors are protected from nullcheck

Only newly discovered adaptive factors are automatically nullchecked.

## ADR-016 — Nullcheck new initialized-fit factors from high to low

This preserves stable user factor IDs when native factor removal shifts
subsequent indices.

## ADR-017 — GBCD uses `kmax`, not `max_factors`

Upstream Kmax is not a strict final rank ceiling.

## ADR-018 — GBCD Baseline is a factor, not an offset

Its membership varies by cell and its covariance prior is signed.

## ADR-019 — GBCD global activity support is conservative

Because Baseline and GEPs have different activity support, the global
two-state/three-state representation declaration uses `"unknown"` and
component-specific support is stored separately.

## ADR-020 — GBCD effect and uncertainty stages remain distinct

EB-SNMF feature effects and conditional LM+ash z/lfsr are preserved with their
different statistical origins.

## ADR-021 — GBCD gets a memory preflight

The wrapper should prevent accidental dense cell × cell covariance allocation
beyond an explicit budget when the pinned upstream branch predicts that path.

## ADR-022 — Multistart owns initialization

Caller-supplied greedy init callbacks conflict with multistart and are rejected.

## ADR-023 — Multistart preserves sign constraints explicitly

A seeded custom wrapper around `flash_greedy_init_default()` must pass the
correct sign constraints because overriding `init_fn` bypasses FLASH's automatic
constraint injection.

## ADR-024 — softImpute/IRLBA multistart starts are signed-only in v1

Do not claim constraint-safe nonnegative initialization from helpers that do not
accept sign constraints.

## ADR-025 — ELBO is an optimization criterion

Highest ELBO is not labeled "biologically best."

## ADR-026 — Group prevalence extends existing APIs

No new parallel group-structure API.

## ADR-027 — Split-signed initialization remains experimental

Native semi-NMF is the baseline.

## ADR-028 — Ordinary old workflows retain basis identity

New metadata is not injected into unchanged zero-offset ordinary
representations.

Projection basis IDs may change only where the metadata itself is corrected.
