# API and Data Contracts

This file is normative for object schemas and high-level function behavior.

---

# 1. Naming and identity

## Source member ID

```text
<dataset_id>::<factor_id>
```

Both components must be non-empty and `dataset_id` must be unique across source datasets.

## Reference program ID

Default deterministic IDs:

```text
C001
C002
...
```

Program IDs are assigned after final clustering using a deterministic ordering:

1. decreasing `n_datasets`;
2. decreasing `median_similarity`;
3. lexical smallest member ID.

The implementation may instead use stable hash-backed IDs if that improves incremental reference updates, but the choice must be versioned and tested.

---

# 2. Source table

Required source-level table:

```text
dataset_id
basis_id
fit_id
feature_scale
preprocessing_hash
activity_support
effect_support
n_factors
n_features
```

Rules:

- one row per source dataset;
- unique `dataset_id`;
- unique `basis_id` is preferred but not required if the same basis is intentionally reused;
- no source sample-level activity matrix is embedded by default.

---

# 3. Member eligibility table

```text
member_id
dataset_id
factor_id
basis_id
activity_support
effect_support
n_features
zero_activity
zero_effect
eligible
reason
```

`reason` is empty only when `eligible=TRUE`.

Allowed support values follow current package semantics:

```text
signed
nonnegative
unknown
```

---

# 4. Comparison manifest

Minimum schema:

```r
list(
  mode = "generic" | "starcat2025",
  feature_scale = <non-empty string>,
  preprocessing_compatibility = <non-empty string>,
  profile_normalization = "auto" | "l1" | "l2",
  feature_standardization = "none" | "source_sd",
  algorithm_version = "1.0.0"
)
```

For `starcat2025`:

```text
profile_normalization must be compatible with nonnegative fixed-sum spectra
feature_standardization = source_sd
source_feature_sets required
source_feature_sd required
metric = pearson
```

The manifest must be serializable without callbacks/functions.

---

# 5. Comparison-profile transform record

For every member factor store:

```text
member_id
normalization
normalization_constant
orientation
anchor_feature
sign_multiplier
feature_standardization
```

`sign_multiplier` may be `-1` only if both activity and feature supports are signed.

---

# 6. Pairwise similarity table

Required columns:

```text
member_i
dataset_i
factor_i
member_j
dataset_j
factor_j
metric
raw_similarity
similarity
orientation
shared_features
coverage_i
coverage_j
eligible
reason
```

Rules:

- each unordered pair occurs once;
- `member_i < member_j` under lexical ordering for storage;
- `raw_similarity` retains sign;
- `similarity` is the clustering score;
- `orientation` is `1`, `-1`, or `NA`;
- same-dataset rows may exist but `eligible=FALSE` for reference graph construction;
- nonfinite similarity cannot be eligible.

---

# 7. Graph contracts

## Directed candidates

```text
from
to
similarity
rank_from
rank_scope
above_threshold
eligible
reason
```

## Accepted undirected edges

```text
member_a
member_b
similarity
reciprocal
rank_a
rank_b
```

Rules:

- no self edges;
- no duplicate undirected edges;
- all endpoints occur in the member table;
- accepted edges are always cross-dataset;
- accepted edge sorting is deterministic.

---

# 8. Merge audit table

Every considered accepted edge should yield one audit row:

```text
step
member_a
member_b
cluster_a_before
cluster_b_before
similarity
proposed_size
proposed_density
dataset_collision
decision
reason
```

Allowed decisions:

```text
merge
already_same_cluster
reject_dataset_collision
reject_density
```

This table is central to reproducibility and debugging.

---

# 9. Cluster table

```text
cluster_id
member_id
dataset_id
factor_id
role
```

`role` is:

```text
member
singleton
excluded_constraint
```

Rules:

- a member belongs to exactly one terminal cluster;
- no ordinary consensus cluster contains two members from the same dataset;
- excluded/incompatible clusters remain auditable.

---

# 10. Program constraints

`program_constraints` is a data frame:

```text
program_id
activity_support
effect_support
```

A reference program is projection-ready only if its consensus members imply one compatible support pair.

In schema v1, a whole `factor_reference` is projection-ready only when all selected reference programs share a single compatible `activity_support`, because the native projection path initializes the fixed programs jointly with one activity-side prior family.

---

# 11. Reference programs table

```text
program_id
cluster_id
n_members
n_datasets
n_source_datasets
support_fraction
mean_similarity
median_similarity
min_similarity
cluster_density
constraint_compatible
singleton
member_stability_median
member_stability_min
```

`support_fraction` denominator is `n_source_datasets`, not an implicit per-feature denominator.

If a future algorithm uses dataset-specific eligibility denominators, it must add a separate named field rather than silently changing this definition.

---

# 12. Reference members table

```text
program_id
cluster_id
member_id
dataset_id
factor_id
basis_id
orientation_to_consensus
similarity_to_consensus
activity_support
effect_support
stability_recovery_rate
stability_median_similarity
```

Stability fields may be `NA`.

---

# 13. Feature diagnostics table

```text
program_id
feature
consensus_effect
member_mean
member_sd
member_mad
min_member_effect
max_member_effect
n_members_observed
```

In schema v1 all final reference features are common to all contributing members, so `n_members_observed` should normally equal `n_members`.

---

# 14. `factor_reference`

Required structure:

```r
list(
  programs,
  program_constraints,
  programs_table,
  members_table,
  feature_table,
  clusters,
  similarity,
  graph = list(
    directed_candidates = ...,
    accepted_edges = ...,
    merge_audit = ...
  ),
  excluded_members,
  feature_sets,
  comparison_manifest,
  projection_manifest,
  projection_ready,
  source_table,
  build_settings,
  reference_id,
  schema_version = "1.0.0",
  provenance
)
```

Class:

```r
c("factor_reference", "list")
```

## Matrix orientation

```text
rows    = reference features
columns = reference programs
```

## Offset

Schema v1 reference offset is implicitly zero and is not estimated by consensus building.

---

# 15. `reference_id`

Deterministic digest over canonical ordered scientific fields.

Before hashing:

- sort source table by `dataset_id`;
- sort members by `member_id`;
- sort programs by deterministic program order;
- preserve program matrix row/column names in canonical order;
- omit timestamps and warnings.

Serialization round-trip must preserve `reference_id` exactly.

---

# 16. Projection manifest

Minimum projection-ready schema:

```r
list(
  feature_scale = <non-empty string>,
  scope = <non-empty string>,
  transformation = <non-empty string>,
  compatibility = <non-empty string>,
  projection_prior_family = <optional serializable identifier>,
  version = "1.0.0"
)
```

Do not store an R function as the only persistent description of the projection prior.

---

# 17. `project_reference()` contract

Required arguments:

```r
reference
query
sample_side
feature_scale
preprocessing
```

Projection must fail before native fitting when:

- `reference$projection_ready` is false;
- feature scale is missing or incompatible;
- preprocessing contract is missing;
- feature map is not one-to-one;
- coverage is below threshold;
- selected reference programs are rank deficient on shared features;
- required projection prior family cannot be resolved;
- dense allocation exceeds budget.

---

# 18. `reference_projection`

Required fields:

```r
list(
  activity,
  activity_normalized,
  zero_sum_activity,
  representation,
  reconstruction,
  coverage,
  missing_reference_features,
  excluded_query_features,
  warnings,
  native_iteration_diagnostics,
  reference_id,
  analysis_metadata,
  provenance
)
```

Class:

```r
c("reference_projection", "list")
```

## Coverage table

```text
n_reference_features
n_query_features
n_shared_features
reference_fraction
query_fraction
```

## Reconstruction metrics

```text
rmse
mae
fit_r2
correlation
n_scored_entries
```

These are in-sample conditional fit diagnostics.

---

# 19. L1 normalized activity

When requested and allowed:

\[
A^*_{ik} = A_{ik}/\sum_k A_{ik}
\]

Rules:

- all projected activity support must be non-negative;
- negative posterior means outside numerical tolerance are an error;
- zero-sum rows remain zero;
- `zero_sum_activity` marks those rows;
- raw activity remains unchanged.

---

# 20. Residual-program return contract

If implemented after M7 GO:

```r
list(
  reference_projection,
  residual_representation,
  combined_representation,
  reconstruction_comparison,
  residual_evidence,
  warnings,
  analysis_metadata,
  provenance
)
```

## Reconstruction comparison

```text
known_sse
full_sse
delta_sse_residual
known_fit_r2
full_fit_r2
delta_fit_r2_residual
```

## Residual evidence

```text
factor
best_reference_program
similarity
shared_features
stability
factor_metric
status
```

Allowed descriptive statuses:

```text
reference_like_residual
distinct_residual
unstable_residual
ambiguous
```

---

# 21. Error policy

Fail rather than guess for:

- incompatible source/query units;
- generic reference building without an explicit `min_similarity`;
- nonzero consensus-source offsets in schema v1;
- many-to-one feature mappings;
- unsupported sign equivalence;
- mixed incompatible factor constraints within a consensus;
- rank-deficient fixed reference programs;
- unsupported schema versions;
- missing projection contract.

Record/exclude rather than abort the whole build for isolated zero factors or isolated pairwise insufficient overlap, provided enough valid material remains to build the requested reference.

---

# 22. Backward-compatibility contract

No first-release changes to the documented meaning of:

```text
factor_replicability
factor_stability
project_factor_activity
project_factor_features
factor_holdout
factor_representation
```

Shared internal refactoring is allowed only with frozen regression tests.
