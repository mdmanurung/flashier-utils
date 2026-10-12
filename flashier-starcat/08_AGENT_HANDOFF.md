# Coding-Agent Handoff

## Mission

Implement the revised consensus-reference and projection layer in `mdmanurung/flashier-utils` while preserving the package's existing scientific contracts, provenance rules and public APIs.

The coding agent must treat `01_IMPLEMENTATION_SPEC.md`, `04_API_AND_DATA_CONTRACTS.md` and `02_MILESTONES.md` as normative.

---

# 1. Before editing code

Read these repository files in full:

```text
AGENTS.md
DESCRIPTION
NAMESPACE
CREDITS.md
inst/provenance.csv
R/representation.R
R/alignment.R
R/projection.R
R/holdout.R
R/stability.R
R/fit.R
R/provenance.R
tests/testthat/test-alignment.R
tests/testthat/test-projection.R
tests/testthat/test-representation.R
tests/testthat/test-holdout.R
```

Then:

1. record current HEAD;
2. run tests using the repository-prescribed R runtime;
3. run the package check used by the repository;
4. complete M0 before refactoring shared internals.

Do not assume the audit SHA in this bundle is still current.

---

# 2. Core architectural rule

Keep three concepts separate:

```text
source factor_representation   = inferential coordinates
comparison profile             = matching-only normalized feature vector
factor_reference               = final consensus program basis
```

Never overwrite source activity/effect coordinates merely to obtain a comparable profile.

---

# 3. Core implementation sequence

```text
M0 baseline
  -> M1 comparison profiles
  -> M2 all-pairs similarity
  -> M3 graph/clustering
  -> M4 consensus + factor_reference
  -> M5 top-level builder
  -> M6 native projection
  -> M10 validation
  -> M12 release/docs
```

Residual-program work:

```text
M6 -> M7 feasibility -> GO only -> M8
```

Prior transport:

```text
M11 experimental only
```

Do not allow M11 to delay or complicate core projection.

---

# 4. Public API budget

Do not export low-level helpers unless a user-facing need is demonstrated.

Core planned exports:

```r
build_factor_reference()
project_reference()
reference_programs()
reference_members()
reference_diagnostics()
print.factor_reference()
summary.factor_reference()
```

Conditional:

```r
discover_residual_programs()
```

Keep graph/similarity/consensus helpers internal initially.

---

# 5. Scientific fail-closed rules

Stop with an actionable error rather than guessing when:

- feature scale is unknown/incompatible;
- source offsets are nonzero for schema v1;
- preprocessing compatibility is undeclared;
- sign equivalence is requested without both sides signed;
- a feature map is many-to-one;
- a consensus mixes incompatible constraint types;
- reference programs are rank deficient on query overlap;
- query prior family cannot be resolved;
- the reference is not projection-ready;
- a required STARCAT source statistic is missing.

Isolated unusable source factors may be excluded only when their exclusion is recorded and enough valid material remains.

---

# 6. STARCAT implementation rule

Use STARCAT as a conceptual and validation reference.

Do not mechanically port Python source into R.

The specialized `starcat2025` mode must follow the documented contract in this bundle and receive explicit conceptual attribution.

Any discrepancy between:

- the paper;
- the audited STARCAT source;
- the independent R implementation;

must be recorded in `07_SOURCE_AUDIT_AND_DECISIONS.md` or a new decision entry before merge.

---

# 7. Existing alignment code

Important distinction:

- existing `factor_replicability()` uses congruence + one-to-one LSAP assignment;
- reference construction requires all-pairs similarities.

Reuse only shared primitives such as:

- source/view validation;
- feature matching;
- constraint-aware sign handling;
- congruence computation.

Do not use the assignment output as graph input.

Add a regression/adversarial test that would fail if LSAP accidentally pruned graph edges.

---

# 8. Existing projection code

Prefer refactoring toward one shared native fixed-program helper rather than duplicating `.native_projection()`.

Regression requirements:

- `project_factor_activity()` unchanged;
- `project_factor_features()` unchanged;
- `factor_holdout()` unchanged;
- reference projection verifies fixed effects exactly/tolerantly using the same policy as current code;
- warning and iteration diagnostics remain available.

Do not use `flash_update_data()` as an unvalidated shortcut.

---

# 9. Reference object size

Default reference objects must remain lightweight.

Do not serialize:

- full source matrices;
- full source activities;
- complete native fits;
- large stability replicate fits;

unless the caller explicitly requests a debug object.

Store hashes and compact source metadata instead.

---

# 10. Determinism

No RNG should be needed for reference construction.

All of these must be invariant to input ordering:

```text
pairwise similarity rows
accepted edges
cluster membership
program order
reference_id
```

Every tie rule must be explicit.

---

# 11. Residual-program stop condition

Do not begin M8 unless M7 records GO.

If native greedy/backfit cannot preserve reference feature programs cleanly using supported interfaces:

- mark M8 `BLOCKED` or `DEFERRED`;
- create an upstream design note if appropriate;
- ship core consensus/reference projection without residual discovery.

Partial core completion is preferable to fragile native-object manipulation.

---

# 12. Prior transport stop condition

Do not infer a consensus fitted prior from source factors.

M11 must answer both:

1. **technical:** can a prior be extracted/replayed cleanly?;
2. **scientific:** what distribution is justified for a consensus program after source-to-reference rescaling and dataset composition differences?

Without both answers, no public fixed-prior option is added.

---

# 13. Files expected to change

Likely new files:

```text
R/reference-build.R
R/reference-class.R
R/reference-projection.R
R/reference-residual.R       # conditional
```

Likely modified files:

```text
R/alignment.R                # shared pre-assignment primitives only
R/projection.R               # shared projection helper/refactor
R/provenance.R               # if registry plumbing requires
NAMESPACE
DESCRIPTION                   # only if dependencies actually change
inst/provenance.csv
CREDITS.md
vignettes/credits.Rmd
_pkgdown.yml
NEWS.md
```

Do not add an NNLS dependency for the core implementation.

---

# 14. Pull-request sequence

Recommended PRs:

### PR 1 — M0
Regression baseline and fixtures only.

### PR 2 — M1 + M2
Comparison profiles and all-pairs similarity.

### PR 3 — M3
Graph, clustering, audit tables and parity fixture.

### PR 4 — M4 + M5
Consensus estimation, `factor_reference`, builder and accessors.

### PR 5 — M6
Reference projection plus safe shared-internal refactor.

### PR 6 — M9/M10 core
Stability metadata and validation receipts.

### PR 7 — M12 core
Documentation/provenance/release.

### Separate experimental PR
M7/M8 residual programs after GO.

### Separate research branch
M11 prior transport.

---

# 15. Per-task completion record

For every tracker item marked DONE, record:

```text
Task ID:
Files changed:
Public API impact:
Scientific contract:
Tests added/changed:
Benchmark/receipt:
Backward compatibility:
Provenance/credit changes:
Known limitations:
Package check result:
```

Do not mark a task DONE when only code exists but its stated acceptance artifact does not.

---

# 16. Completion check

Before declaring the core project complete:

- run the tracker validator;
- confirm no core task is BLOCKED/TODO;
- run full package checks;
- verify CREDITS and provenance registry consistency;
- verify benchmark receipts are outside package build artifacts;
- verify final reference/projection vignettes build without remote data downloads.
