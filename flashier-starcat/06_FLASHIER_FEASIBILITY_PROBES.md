# `flashier` Feasibility Probes

These probes answer implementation questions before public APIs are committed.

They are executable development experiments, not package features.

---

# A. Core projection capability — already expected to pass

Current `flashier-utils` already demonstrates:

```text
flash_init
-> flash_factors_init
-> flash_factors_fix(feature side)
-> flash_backfit
```

with verification that the fixed program matrix remains unchanged.

The new reference projection layer should reuse this path.

No upstream change is expected for core `project_reference()`.

---

# B. Residual-program discovery probe

## Question B1

Can `flash_greedy()` add new factors to a `flash` object containing reference factors whose feature side is fixed?

### Procedure

1. Simulate `Y = A_ref B_ref^T + A_new B_new^T + E`.
2. Initialize `B_ref` and estimated `A_ref`.
3. Fix the feature side of all reference factors.
4. Save exact reference feature posterior means.
5. Call `flash_greedy(Kmax = ...)` for additional factors.
6. Verify reference feature posterior means unchanged.
7. Verify new factors are added and explain omitted signal.

### Pass

Reference feature programs unchanged within `1e-12` or the existing package's fixed-factor tolerance.

---

## Question B2

Does subsequent `flash_backfit()` respect the existing fixed feature-side flags while updating reference activities and new factors?

### Procedure

Continue B1 and run backfit.

Verify:

- reference feature programs fixed;
- reference activities may update;
- new factors may update on both sides;
- ELBO/reconstruction does not degrade unexpectedly.

### Pass

All fixed reference feature programs remain unchanged.

---

## Question B3

How does `flash_nullcheck()` behave with fixed reference factors plus newly added factors?

### Procedure

Test whether native nullcheck can:

- remove a fixed reference factor;
- renumber factors in a way that loses identity;
- remove only newly added weak factors safely.

### Decision

If fixed reference identity cannot be guaranteed, do not call global nullcheck in `discover_residual_programs()`. Use explicit diagnostics for residual factors instead.

---

## Question B4

Are factor indices and fixed flags stable enough to map reference versus residual factors after greedy fitting/backfit?

### Pass

The implementation can deterministically reconstruct:

```text
reference IDs
residual IDs
native indices
fixed dimensions
```

without inspecting undocumented internal slots beyond currently supported public/native object fields already used by the package.

---

# C. Fitted-prior transport probe — experimental only

Public `flashier` exposes:

```r
flash_fit_get_g(f, n)
```

and `flash_ebnm()` accepts `g_init` with `fix_g=TRUE`.

This makes a fixed-prior experiment technically plausible but does **not** establish a scientifically valid consensus prior.

## Question C1

Can a fitted activity prior be extracted, serialized and replayed for the same prior family?

Test:

```r
g <- flashier::flash_fit_get_g(flashier::flash_fit(fit), activity_dimension)
fn <- flashier::flash_ebnm(prior_family = ..., g_init = g[[k]], fix_g = TRUE)
```

with a minimal one-factor projection.

## Question C2

Can different fixed priors be assigned to different reference factors cleanly using public APIs?

If this requires factor-by-factor initialization, test the exact construction and factor-order stability.

## Question C3

What is the effect of changing factor scale before replaying the prior?

A fitted prior is not meaningful without the exact factor coordinate scale. Demonstrate scale dependence explicitly.

## Question C4

How should multiple source priors be combined for one consensus factor?

No implementation is allowed until a principled answer is benchmarked.

Possible research options:

- no pooling; use query re-estimation;
- fit a new prior to source activities after mapping each member to the final consensus basis;
- dataset-balanced mixture prior;
- hierarchical prior across sources.

These are research choices, not utility glue.

---

# D. Residual-noise freeze

Current package documentation explicitly treats fully frozen trained noise/hyperparameters as unsupported in holdout evaluation.

Do not add a `freeze_noise=TRUE` public option to the new reference layer unless a separate native probe and validation plan is completed.

---

# E. Go/no-go record

Create a development receipt containing:

```text
flashier version
flashier commit if known
ebnm version
probe code hash
fixture hash
B1 PASS/FAIL
B2 PASS/FAIL
B3 behavior
B4 PASS/FAIL
C1-C4 experimental notes
final decision
```

`discover_residual_programs()` may be implemented only when B1, B2 and B4 pass and B3 has a safe documented policy.
