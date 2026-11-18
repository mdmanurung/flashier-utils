# Inference on factors: summaries, tests, contrasts, uncertainty, stability

Read this when the user asks whether a factor differs between groups, how a
group difference looks in feature units, how certain the estimates are, or
whether factors are reproducible.

The examples assume `view <- standardize_factors(fit)` and a `metadata` data
frame with one row per sample: a `sample` column matching the sample IDs, the
biological unit (here `donor`) and the design columns.

## Descriptive group summaries

```r
summarize_factor_activity(view, metadata, by = "group", unit_col = "donor")
contrast <- compare_factor_groups(view, metadata, group = "group",
                                  contrast = c("A", "B"), unit_col = "donor")  # A minus B
plot_factor_activity(view, metadata = metadata, group = "group", style = "distribution")
```

Summaries weight observations by default; cells are not balanced by donor.
They describe; they do not test.

## Association tests

```r
association <- test_factor_association(view, metadata, formula = ~ group + age,
  engine = "lm", unit_col = "donor", independence = "independent")
association$table   # one row per factor x term: estimate, std.error, statistic, df, p.value, adjusted p
```

| Design | `engine` | `independence` | Formula |
|---|---|---|---|
| One sample per donor | `"lm"` (or `"limma"` with tens of factors) | `"independent"` | `~ group` |
| Several samples per donor | `"dream"` (needs variancePartition) | `"repeated"` | `~ group + (1 | donor)` |

- The random-effect grouping must map one to one onto `unit_col`;
  `(1 | donor:visit)` is refused.
- dream reports its own Satterthwaite df and p; `moderate = TRUE` adds
  empirical Bayes moderation (pooled df, anti-conservative with few factors).
- limma always moderates and needs factor residual variances on a similar scale;
  use `"lm"` for raw-scale views.
- Multiple testing: `adjust = "BH"` over `adjust_scope = "all_tests"` by default.
- Pass `contrasts =` a numeric matrix (rows named by design terms) for custom contrasts.
- Non-estimable rows get `p.value = NA`, not 0.

## Back-project a coefficient into feature units

```r
reconstructed <- backproject_contrast(view, association$table,
  input = "coefficients", factor_basis = view$manifest, term = "groupB")
reconstructed$effects   # per-feature B minus A
```

- `factor_basis` must be the manifest of the view the test ran on.
- Point coefficients carry no uncertainty. For intervals pass
  `coefficient_covariance` (from `control = list(joint_coefficient_covariance = TRUE)`
  in `lm`) or compatible draws; factor posterior SDs cannot substitute.
- `validate_backprojection()` compares against an observed feature contrast with
  matching direction, covariates and units. Agreement is descriptive.

## Posterior uncertainty

```r
factor_uncertainty(fit, side = "effects", sample_side = "rows",
                   representation = view, include_lfsr = TRUE)
factor_intervals(fit, side = "effects", sample_side = "rows",
                 representation = view, level = 0.9, nsamp = 200, seed = 1)
```

Pass the **native fit** plus `representation = view`: the view supplies the
coordinates, the fit supplies the posterior. These quantities are conditional
on the estimated priors and rank; they are not standard errors of a regression
and say nothing about the biological identity of a factor. `factor_intervals()`
warns when fewer than 5 draws fall in a tail.

## Stability under resampling of biological units

```r
stab <- factor_stability(X, reference_fit = fit, metadata = metadata,
  unit_col = "donor", seed = 4, B = 100,
  scheme = "subsample", fraction = 0.8, min_congruence = 0.9)
stab$summary
```

- `scheme = "subsample"` avoids duplicated units, which bias empirical Bayes
  estimates under `"bootstrap"` (the default; it warns and documents this).
- `fit_fn = NULL` replays the reference fit's recorded `fit_ebmf()` recipe;
  custom or edited fits need a callback
  `function(X, metadata, sample_side, seed, settings, observation_aux)`.
- `workers > 1` uses forked workers on Unix; `checkpoint_dir` makes long runs resumable.

## Held-out prediction by biological unit

```r
folds <- data.frame(sample = metadata$sample,   # all samples of a donor share a fold
  fold = ifelse(metadata$donor %in% test_donors, "test", "train"))
callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
  fit_ebmf(X, sample_side, max_factors = 10, seed = seed, verbose = 0,
           feature_scale = "log1p_normalized_counts")
}
factor_holdout(X, metadata = metadata, unit_col = "donor", folds = folds,
  fit_fn = callback, seed = 4, evaluation = "masked",
  activity_features = genes_a, scoring_features = genes_b,
  control = list(sample_side = "rows", preprocessing_manifest = list(scope = "identity"),
                 evaluated_folds = "test"))$metrics
```

- `evaluation = "projection"` scores the same entries used to estimate test
  activities (optimistic); `"masked"` scores disjoint withheld features or entries.
- Metrics include `prediction_r2` and `skill` (1 − SSE / SSE of the training-mean
  baseline; 0 means no better than the baseline).
- Preprocessing fitted on all data makes the result conditional; use
  `preprocess_fn = list(fit = ..., transform = ...)` to fit it on training folds.

## Matching factors between fits

```r
factor_replicability(view_a, view_b, min_congruence = 0.9)$matches   # one-to-one, Tucker congruence
compare_prior_fits(view_nmf, view_semi)                              # coverage, residual profiles
```

Both need matching `feature_scale`. `factor_replicability()` is pairwise
one-to-one assignment; for many datasets use `build_factor_reference()`
(see projection-and-references.md).
