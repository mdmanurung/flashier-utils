# Fitting variants and preprocessing

Read this when the user needs batch correction before fitting, nonnegative or
semi-nonnegative programs, GBCD, grade-of-membership DE on counts, membership
filtering, or a Seurat embedding of factors.

## Batch correction before fitting

```r
corrected <- correct_feature_batch(X, "rows", metadata, batch_col = "batch",
                                   assignments = "celltype", protect = ~ condition)
```

This is Harmony's per-cluster ridge regression (Korsunsky et al. 2019), applied
to the feature matrix as in cNMF `Preprocess`. Here it is independently written
in R. `assignments` is a metadata column of hard clusters or a soft
samples × clusters matrix. `harmony_assignments(pcs, metadata, "batch", seed = 1)`
gets soft clusters from the harmony package. `protect` keeps covariates of
interest from being removed. Values are not clipped, so negative values can
appear; use a signed prior on that side.

## Prior families and supports

`fit_ebmf(X, sample_side, ebnm_fn = ...)`. A single function applies to both
sides. A list of two follows native row/column order.

| Programs | Row prior (samples as rows) | Column prior |
|---|---|---|
| Signed (PCA-like) | `ebnm_normal` or `ebnm_point_laplace` | `ebnm_point_laplace` |
| Semi-NMF: nonnegative activity, signed effects | `ebnm_point_exponential` | `ebnm_point_laplace` |
| NMF | `ebnm_point_exponential` | `ebnm_point_exponential` |

- `S` supplies known noise SDs. `var_type` sets which residual variances are
  estimated: 0 = one, 1 = per row, 2 = per column.
- `max_factors` is the greedy upper bound.
- `intercept = "feature"` adds a per-feature intercept, which shows up as the
  view's `offset`.
- Compare seeds and initializations; agreement between them is not calibration.

## GBCD gene programs

```r
programs <- fit_gbcd_programs(X, "rows", kmax = 5, feature_scale = "natural_log",
  preprocessing = list(scope = "within_dataset",
                       description = "log1p library-size normalized counts", log_base = "e"))
view <- programs$representation          # factors are GEP1, GEP2, ...
programs$signatures; programs$diagnostics
```

- Needs `gbcd` (stephenslab) at the audited commit 670e3af.
- `kmax` controls initialization only; it does not cap the final rank.
- `z_score` and `lfsr` come from separate gene-wise regression plus ashr,
  conditional on the memberships.
- No native PVE or posterior sampling is available, and projection priors
  cannot be inferred.

## Grade-of-membership DE on counts

```r
test_gom_de(counts, memberships, sample_side = "rows", seed = 1)
```

This calls `fastTopics::de_analysis()`. Counts must be nonnegative integers,
and memberships must be nonnegative with rows summing to one. The default
`lfc_stat = "le"` reports least-extreme log2 fold changes.
`factor_features(select = "distinctive")` is a descriptive margin, not GoM DE.

## Filter cells by total membership

```r
filtered <- filter_cells_by_total_membership(nonnegative_fit, max_val = 10, numiter = 12)
filtered$removed_samples; filtered$history
```

This follows the ImmGenT-GP rule: repeatedly drop samples whose total
max-scaled membership exceeds `max_val`. Activities must be nonnegative. The
rule selects broad membership; it does not detect poor cell quality.

## Seurat and embeddings

- `plot_factor_embedding(view, reduction = umap_matrix, factors = c("F1", "F2"))` colours an embedding by
  factor activity.
- With a Seurat object, the object's cell IDs must equal the activity sample IDs.
