# Revised Test and Benchmark Plan

The validation strategy separates **exact software invariants** from **scientific performance benchmarks**.

Scientific benchmark receipts remain outside package builds, following repository instructions.

---

# 1. Exact unit invariants

These are release-blocking.

## Source immutability

After reference building:

```r
identical(serialize(source_before, NULL), serialize(source_after, NULL))
```

for source representations/fits used in test fixtures.

## Factor-scale invariance

For any positive scalar `c`:

```text
A_k -> c A_k
B_k -> B_k / c
```

must not change:

- comparison profile;
- pairwise similarity;
- reference cluster assignment;
- final consensus program, except numerical tolerance.

## Factor permutation invariance

Permuting factors within a source dataset must not change final memberships or `reference_id`.

## Dataset-order invariance

Permuting source dataset input order must not change graph edges, clusters, programs or `reference_id`.

## Serialization invariance

`saveRDS()` -> `readRDS()` preserves `reference_id` and all scientific content.

---

# 2. Sign/constraint unit tests

Cases:

1. signed activity + signed effects: sign-reversed factors may align;
2. nonnegative activity + signed effects: sign reversal is not allowed;
3. nonnegative/nonnegative: negative effects fail profile validation;
4. unknown support: no automatic sign flip;
5. mixed incompatible constraints in one proposed cluster: excluded from final projection-ready programs and recorded.

---

# 3. Feature-policy tests

Test independently:

- exact common feature basis;
- pairwise intersection;
- explicit reference feature set;
- source discovery/HVG sets;
- STARCAT pairwise union-of-source-sets restricted to shared features;
- final source-set union restricted to globally shared features;
- duplicate/mapped feature rejection;
- zero or nonfinite source SD in `starcat2025` mode.

---

# 4. Similarity tests

Hand-constructed vectors:

- identical -> expected 1;
- orthogonal -> expected 0 for cosine/congruence;
- sign reversed -> expected -1 raw and +1 only when sign equivalence allowed;
- linear affine-shift examples distinguish Pearson from congruence;
- low overlap -> ineligible with explicit reason.

Confirm graph similarity path never invokes LSAP assignment.

---

# 5. Graph/clustering adversarial tests

Construct tiny graphs where the expected result is exact.

Required cases:

- simple two-dataset pair;
- triangle with all edges;
- chain A-B-C with missing A-C edge;
- cluster merge just below density threshold;
- cluster merge exactly at density threshold;
- duplicate-dataset collision;
- equal-similarity tie;
- isolated singleton;
- multiple disconnected consensus groups.

Freeze `density >= 2/3` as the schema-v1 STARCAT-mode rule and test the exact two-thirds case explicitly.

---

# 6. STARCAT graph parity benchmark

Purpose: validate the independently written graph implementation against STARCAT's published/catalog code path without copying its source into the R package.

Run outside package checks.

Use identical prepared program profiles / similarity inputs and compare:

- directed candidate neighbors;
- reciprocal accepted edges;
- terminal memberships;
- singleton assignments.

The benchmark must record:

```text
STARCAT commit
flashier-utils commit
profile inputs hash
threshold
rank_max
density rule
membership diff
```

Any intentional deviation is documented in `07_SOURCE_AUDIT_AND_DECISIONS.md`.

---

# 7. Exact synthetic consensus recovery

Generate noiseless or near-noiseless factors with separated shared programs and dataset-specific programs.

At the exact/noiseless setting require:

- no false shared merges;
- all intended shared programs recovered;
- dataset-specific programs remain singleton/unreferenced by default;
- consensus profiles equal the expected normalized mean within tolerance.

Then stress with increasing noise and record degradation rather than imposing an arbitrary universal threshold.

---

# 8. Split/merged-factor stress test

Simulate:

- one true program recovered once in dataset A;
- split into two factors in dataset B;
- merged with another program in dataset C.

Questions:

- does one-member-per-dataset prevent dataset B dominating support?
- are ambiguous factors visible in graph/merge diagnostics?
- does density filtering prevent chain-like overmerging?

No automatic split/merge repair is required in the first release.

---

# 9. Preprocessing/scale mismatch stress test

Construct source profiles that have the same latent biology but incompatible feature-wise scaling.

Test that:

- a missing comparison contract fails;
- declared incompatible feature scales fail;
- `starcat2025` source-SD standardization restores the intended profile comparison in a controlled example;
- generic mode does not silently apply STARCAT-specific scaling.

---

# 10. Partial feature-overlap benchmark

Remove reference features at:

```text
0%, 10%, 25%, 50%, 75%
```

Measure:

- activity correlation to full-feature projection;
- RMSE / fit R2;
- rank failures;
- program conditioning;
- coverage diagnostics.

Below `min_coverage`, projection must fail explicitly.

---

# 11. Missing-reference-program projection simulation

Generate query data with one true program omitted from the reference.

Core projection expectations:

- known-program activity remains recoverable where identifiable;
- reconstruction degrades relative to complete reference.

If residual discovery is implemented:

- combined reconstruction improves materially;
- at least one residual factor matches the omitted truth more strongly than any included reference program;
- reference feature programs remain unchanged.

---

# 12. Extra-reference-program simulation

Include reference programs absent from query.

Compare activity assigned to:

- truly present programs;
- absent reference programs.

The benchmark evaluates shrinkage behavior under the chosen projection prior family but does not impose a universal zero threshold.

---

# 13. Small-query benchmark

Query sizes:

```text
1, 5, 10, 25, 50, 100, 500
```

Measure:

- activity recovery when simulation truth is available;
- reconstruction;
- convergence behavior;
- runtime;
- sensitivity to projection prior family.

This benchmark motivates, but does not presuppose, experimental fixed-prior transport.

---

# 14. Signed and semi-nonnegative benchmark

Simulate separately:

- signed/signed factors;
- nonnegative activity + signed feature effects;
- fully nonnegative factors.

Verify:

- correct sign-equivalence rules;
- correct profile normalization;
- consensus cluster recovery;
- projection support;
- L1 activity normalization only for nonnegative activity.

---

# 15. Leave-one-dataset-out real-data validation

For each source dataset:

1. construct the reference using all other datasets;
2. project the held-out dataset;
3. compare projected reference activities with independently discovered held-out programs;
4. compare profile matches and reconstruction;
5. rebuild every data-derived reference quantity without using the held-out dataset.

This is a strict anti-leakage requirement.

Primary real-data question:

> Do programs with stronger independent cross-dataset support transfer more reproducibly to a held-out dataset?

Use descriptive/association analysis with appropriate biological units; do not treat cells as independent replicates for inferential claims.

---

# 16. Pancreas multi-technology benchmark

Use the existing Jamboree-derived pancreas resources already supported by the repository where feasible.

Potential design:

- independently fit factors in multiple technologies;
- build a technology-spanning reference;
- hold out one technology;
- test program transfer and cell-type-associated activity patterns.

This benchmark is preferred for development because the package already has data preparation/provenance machinery for these resources.

It remains outside package checks if the data are downloaded.

---

# 17. STARCAT projection comparison

Do not add a production NNLS dependency solely for parity.

External benchmark compares:

```text
STARCAT fixed-program projection
flashier-utils native EBMF projection
```

on a nonnegative reference/query setting.

Report:

- raw usage/activity association;
- row-normalized usage association where appropriate;
- top-program agreement;
- reconstruction;
- query-size behavior;
- absent-program behavior;
- partial-gene-overlap behavior.

Interpret differences as method differences, not implementation errors, unless the same mathematical objective is intentionally being tested.

---

# 18. Residual-program stability benchmark

Only if M8 is implemented.

Repeatedly resample independent biological units when available and refit residual factors.

Report:

- residual recovery rate;
- median matched similarity;
- reference-like versus distinct residual classification stability.

A one-dataset residual factor is never called a consensus biological program solely because it has high PVE or reconstruction improvement.

---

# 19. Runtime and memory benchmark

Scale:

```text
datasets:          2-20
factors/dataset:   5-100
features:          500-20,000
query samples:     10-100,000
```

Measure separately:

- comparison profile construction;
- all-pairs similarity;
- graph/clustering;
- consensus estimation;
- projection;
- residual discovery if enabled.

Respect the existing dense-memory budget contract.

---

# 20. Core release gates

Core release fails if any of the following is true:

- source representation changes during reference construction;
- factor rescaling changes cluster membership;
- input order changes cluster membership or `reference_id`;
- a cluster contains multiple factors from the same dataset;
- consensus programs mix incompatible factor supports;
- final reference feature universe is implicit;
- projection changes fixed reference programs;
- query preprocessing/feature scale can be omitted;
- raw activity is silently replaced by normalized activity;
- legacy alignment/projection regression tests fail;
- serialization changes scientific content;
- provenance cannot trace a consensus program to source factors.

Scientific release additionally requires completed receipts for:

- exact synthetic recovery;
- graph parity;
- feature-overlap stress;
- missing/extra program projection;
- signed/semi-NMF simulation;
- leave-one-dataset-out real-data validation.
