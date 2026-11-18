# Projection onto new data and consensus references across datasets

Read this when the user wants to score fitted programs in a new dataset, add a
new assay to fitted samples, or build one program catalog from several
independently fitted datasets (the STARCAT idea, done with EBMF).

## Score fitted programs on new samples

Feature effects stay fixed; query activities are estimated by native flashier.

```r
projected <- project_factor_activity(fit, query, sample_side = "rows",
  feature_scale = "log1p_normalized_counts",          # must equal the reference
  preprocessing = list(scope = "within_dataset",
                       description = "log1p library-size normalized counts"),
  feature_match = "intersection", min_coverage = 0.8, maxiter = 200)
projected$representation$activity     # query samples x factors
projected$reconstruction              # in-sample fit_r2, rmse; not prediction
```

- The query must be preprocessed like the reference; the call only checks the
  declared labels.
- `feature_match = "exact"` (default) requires every reference feature;
  `"intersection"` uses shared features down to `min_coverage`. Map gene names
  with `gene_mapping = c(query_name = "reference_name")`.
- Needs at least two query samples. Rank-deficient shared support is refused.
- `prior_policy = "query"` (default) re-estimates activity priors on the query,
  so a sample's activity depends on the rest of the query set.
  `prior_policy = "fixed_reference"` reuses the reference fit's priors (needs
  the native fit, not a view). Only `"fixed_reference"` plus sample-wise
  `var_type` (1 for row samples, 2 for column samples) makes each sample
  independent of the others.

New features measured on the same samples (for example protein after RNA):

```r
project_factor_features(view, new_assay, sample_side = "rows",
  feature_scale = "clr_adt", ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace))
```

## Consensus references across datasets (experimental)

`build_factor_reference()` compares the factors of several independently fitted
datasets pair by pair, links similar factors in a graph, merges greedily with at
most one factor per dataset, and averages each group into a consensus program on
a common feature basis. Conceptual method: STARCAT (Kotliar et al. 2025, Nature
Methods 22:1964); independent R code. **Experimental:** arguments may change.

```r
fits <- list(cohortA = fit_a, cohortB = fit_b, cohortC = fit_c)   # names = dataset IDs
ref <- build_factor_reference(fits,
  comparison_manifest = list(mode = "generic",
                             feature_scale = "log1p_normalized_counts",
                             preprocessing_compatibility = "user_declared"),
  min_similarity = 0.8,                     # required in generic mode; no universal value
  projection_manifest = list(feature_scale = "log1p_normalized_counts",
                             scope = "external_reference", transformation = "log1p",
                             compatibility = "exact_declared_contract",
                             activity_prior = "ebnm_point_exponential"))
ref                                         # programs, singletons, projection readiness
reference_programs(ref)                     # features x programs
reference_members(ref)                      # each source factor: cluster, program, eligibility
reference_diagnostics(ref, "programs")      # support, mean/min similarity, member fit
reference_diagnostics(ref, "merges")        # every merge decision and its density
```

| Mode | Similarity | Edges | Merge rule | Consensus |
|---|---|---|---|---|
| `"generic"` | congruence (or `"pearson"`) on shared features | all pairs above `min_similarity`; `rank_max`, `reciprocal` optional | merged-set density ≥ `min_cluster_density` (2/3) | mean L1 (nonneg effects) or L2 profile, renormalized |
| `"starcat2025"` | Pearson on source-SD-scaled spectra, union of source feature sets | reciprocal top `rank_max` (default n datasets), > 0.5, best per dataset | starCAT cross-cluster density > 0.666 | mean of SD-scaled spectra, no renormalization |

- Sources need a declared `feature_scale` matching the manifest and zero offsets
  (no intercept fits).
- `starcat2025` needs nonnegative effects plus `source_feature_sets` and
  `source_feature_sd`: lists named by dataset with each dataset's overdispersed
  genes and named per-gene SD.
- Sign flips are allowed only between factors with signed activity and effects.
  Semi-NMF and NMF factors keep their sign.
- Clusters spanning fewer than `min_datasets` (default 2) datasets stay as
  singletons in the diagnostics. Clusters that mix supports are recorded as
  `constraint_incompatible`.
- `reference_id` hashes the scientific content and ignores dataset and factor
  order. A modified object fails validation.

Project a query onto the consensus programs:

```r
out <- project_reference(ref, query, sample_side = "rows",
  feature_scale = "log1p_normalized_counts",   # must equal projection_manifest$feature_scale
  preprocessing = list(scope = "identity"),
  feature_match = "intersection", normalize_activity = "l1")
out$activity; out$activity_normalized; out$reconstruction
```

- A reference built without `projection_manifest` is for inspection only and
  is refused.
- The activity prior is the `ebnm_fn` argument or `projection_manifest$activity_prior`.
  Its support must match the programs' activity support. No source prior is
  inherited.
- `normalize_activity = "l1"` is allowed only with a nonnegative prior; raw
  activity is always returned too.
- Cross-dataset recurrence is descriptive. It does not prove a biological
  program; check held-out transfer before claiming one.
