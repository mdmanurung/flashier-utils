---
name: flashier-utils
description: Interpret, test and transfer empirical Bayes matrix factorizations fitted with the R package flashier, using the R package flashier.utils. Use when the user fits flashier/EBMF/EBNMF factors (or GBCD gene programs) to a named matrix such as single-cell or bulk expression, and wants a named activity/effects view, top features per factor, activity plots by metadata, factor-metadata association tests, a group contrast back-projected into feature units, posterior uncertainty, resampling stability, held-out prediction, projection of fitted programs onto new samples, or a consensus program reference across several datasets (STARCAT-like). Do not use for cNMF, NMF in Python, MOFA, or Python packages; do not use when no flashier fit or EBMF is involved.
license: MIT
---

# flashier.utils

flashier.utils turns a flashier fit into a **view**: a samples × factors
`activity` matrix and a features × factors `effects` matrix that share factor IDs
and reconstruct the fit. Every downstream function takes a native fit or a view,
records its coordinate basis, and refuses to mix bases.

Website: <https://mdmanurung.github.io/flashier-utils/> (full docs for agents:
<https://mdmanurung.github.io/flashier-utils/llms.txt>).

## Setup

```r
# R >= 4.3; flashier >= 1.0.60 (CRAN flashier 1.0.7 is too old)
remotes::install_github("mdmanurung/flashier-utils")
library(flashier.utils)
```

Optional engines are in Suggests and are checked at call time: psych and clue
(`factor_replicability()`, `build_factor_reference()` comparisons), limma and
variancePartition (`test_factor_association()` engines), gbcd, fastTopics, harmony.

## Core workflow

```r
fit <- fit_ebmf(X, sample_side = "rows",                 # rows of X are samples
  ebnm_fn = list(ebnm::ebnm_point_exponential,          # row prior (= samples here)
                 ebnm::ebnm_point_laplace),             # column prior (= features here)
  max_factors = 20, backfit = TRUE, seed = 1, verbose = 0,
  feature_scale = "log1p_normalized_counts")            # label only; never applied
view <- standardize_factors(fit)                        # default scaling = "ldf"
plot_factor_pve(fit)                                    # look first
factor_features(view, factor = "F1", select = "distinctive", n = 20)
plot_factor_activity(view, metadata = meta, group = "celltype")
```

A bare `flashier::flash()` fit works too: pass `sample_side` (and
`feature_scale`) to `standardize_factors()`. `fit_ebmf()` records them for you.

### Declarations the package will not guess

| Declaration | Where | Why it is required |
|---|---|---|
| `sample_side` (`"rows"`/`"columns"`) | `fit_ebmf()`, bare fits | Orientation is never inferred |
| Unique row and column names | the input matrix | Results are joined by ID, never by position |
| `feature_scale` (free text) | `fit_ebmf()` / `standardize_factors()` | Comparisons and projection need matching units; `"unknown"` is refused there |
| `unit_col` + `independence` | association, stability, holdout | The biological unit (e.g. donor), not the cell, is the replicate |
| `seed` | `factor_stability()`, `factor_holdout()`, `test_gom_de()` | Required integer; no hidden RNG |
| `preprocessing = list(scope = ...)` | `project_factor_activity()`, `project_reference()` | You declare the query was processed like the reference |

### Coordinates

`standardize_factors(fit, scaling = ...)`: `"ldf"` (default, unit-norm factors,
scale on effects), `"raw"` (flashier's own), `"rms"` (unit RMS activity per
factor, comparable across factors). All reconstruct the same matrix. Keep the
activity and effects of one view together; never combine coefficients from
different views. Reconstruct with
`sweep(tcrossprod(view$activity, view$effects), 2, view$offset, "+")`.

### Choosing priors

`ebnm_fn = list(row_prior, column_prior)` is in **native row/column order, not
sample/feature order**; with `sample_side = "columns"` reverse it. Nonnegative
activity plus signed effects (semi-NMF) suits expression programs: with samples
as rows, `list(ebnm::ebnm_point_exponential, ebnm::ebnm_point_laplace)`. Signed both sides:
`ebnm_normal` / `ebnm_point_laplace`. The declared support decides what is
allowed later (sign flips, L1 activity normalization, anchors).

## Task map

| Goal | Function | Reference |
|---|---|---|
| Top features per factor | `factor_features()` (`select = "largest"`, `"distinctive"`, `"both"`) | this file |
| Plots | `plot_factor_pve`, `plot_factor_weights`, `plot_factor_features`, `plot_factor_activity`, `plot_factor_scatter`, `plot_factor_data`, `plot_factor_embedding` | this file |
| Matrices | `factor_activity()`, `factor_effects()` | this file |
| Group summaries and contrasts | `summarize_factor_activity()`, `compare_factor_groups()` | [inference.md](references/inference.md) |
| Factor ~ metadata tests | `test_factor_association()` (`"lm"`, `"limma"`, `"dream"`) | [inference.md](references/inference.md) |
| Contrast in feature units | `backproject_contrast()`, `validate_backprojection()` | [inference.md](references/inference.md) |
| Posterior SD, lfsr, intervals | `factor_uncertainty()`, `factor_intervals()` | [inference.md](references/inference.md) |
| Stability, held-out prediction, matching fits | `factor_stability()`, `factor_holdout()`, `factor_replicability()`, `compare_prior_fits()` | [inference.md](references/inference.md) |
| Score fitted programs on new samples | `project_factor_activity()`, `project_factor_features()` | [projection-and-references.md](references/projection-and-references.md) |
| Consensus programs across datasets (STARCAT-like) | `build_factor_reference()`, `project_reference()`, `reference_*()` | [projection-and-references.md](references/projection-and-references.md) |
| Batch correction, GBCD, GoM DE, membership filter | `correct_feature_batch()`, `fit_gbcd_programs()`, `test_gom_de()`, `filter_cells_by_total_membership()` | [fitting.md](references/fitting.md) |
| Who wrote what | `factor_provenance()` | this file |

## Traps

- **Sign is arbitrary for signed factors.** Flipping a factor on both sides
  leaves the fit unchanged, so a negative effect is not repression. Only
  semi-NMF or NMF supports fix the sign.
- **Contrast direction differs by function.** `compare_factor_groups(contrast = c("A", "B"))`
  is A − B. An `lm` term `groupB` from `test_factor_association()` is B − A.
- **Cells are not replicates.** With several samples per donor, use
  `independence = "repeated"` with `engine = "dream"` and the donor as a random
  effect, e.g. `formula = ~ group + (1 | donor)`, `unit_col = "donor"`.
- **limma moderation needs comparable scales.** Use `engine = "lm"` for raw-scale
  views; limma warns when factor residual variances differ more than tenfold.
- **Posterior SDs are conditional** on the fitted priors and rank (mean-field,
  usually too narrow). They are not regression standard errors; keep the native
  fit, because a view alone has no posterior.
- **Back-projection needs the matching basis**: pass `factor_basis = view$manifest`
  from the same view the coefficients were fitted on.
- **`feature_scale` is never applied.** Normalize and log-transform before fitting.
- **Projection activities depend on the query set** (priors and pooled noise are
  re-estimated). For per-sample independence use
  `prior_policy = "fixed_reference"` with the native fit and sample-wise
  `var_type` (1 when samples are rows, 2 when columns).
- **Display normalization is not inference.** `display_normalization` in plots
  only rescales colours.
- Errors are deliberate refusals (unknown units, mixed bases, missing unit
  column). Fix the declaration; do not work around the check.

## Credits

Native FLASH: Jason Willwerscheid, Peter Carbonetto, Wei Wang, Matthew Stephens.
Use `factor_provenance("<function>")` for audited sources of any function and
cite the upstream methods it reports.
