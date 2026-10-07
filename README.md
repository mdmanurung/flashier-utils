# flashier utils

Interpret a [flashier](https://willwerscheid.github.io/flashier/) matrix
factorization. See where factors are active,
identify their strongest features and relate them to sample metadata.

## Installation

Requires R 4.3 or later.

```r
install.packages("remotes")
remotes::install_github("mdmanurung/flashier-utils")
```

## Quick start

The GTEx example bundled with flashier contains association z-scores across
tissues. Here, tissues are the samples and variant-gene pairs are the features.

```r
library(flashier.utils)
data("gtex", package = "flashier")

fit <- fit_ebmf(gtex, sample_side = "columns", max_factors = 3,
                backfit = TRUE, seed = 1, verbose = 0,
                feature_scale = "association_z_score")
view <- standardize_factors(fit)

plot_factor_activity(view)
features <- factor_features(view, factor = "F1", select = "largest", n = 5)
features[, c("feature", "estimate")]
```

`view$activity` contains one row per sample. `view$effects` contains one row
per feature. Matching columns describe the same factor.

For your own data, supply a numeric matrix with unique row and column names.
Choose `sample_side` explicitly and preprocess the data before fitting.
`feature_scale` records the input units. The package does not transform them.

## Tutorials

Choose the task you want to accomplish. Start with **Fit and extract factors**
if you are new to the package.

| Task | Vignettes |
| --- | --- |
| Fit and extract factors | [Fit and interpret factors](https://mdmanurung.github.io/flashier-utils/articles/first-fit.html); [Factor scales and measurement units](https://mdmanurung.github.io/flashier-utils/articles/representation-units.html); [Reorder, remove, fix and refit factors](https://mdmanurung.github.io/flashier-utils/articles/edit-fits.html) |
| Identify features and pathways | [Find defining features](https://mdmanurung.github.io/flashier-utils/articles/largest-distinctive.html); [Annotate factors with pathways](https://mdmanurung.github.io/flashier-utils/articles/enrichment.html) |
| Compare groups and conditions | [Compare groups and test associations](https://mdmanurung.github.io/flashier-utils/articles/contrasts.html); [Translate contrasts into feature effects](https://mdmanurung.github.io/flashier-utils/articles/backprojection.html); [Compare effects across conditions](https://mdmanurung.github.io/flashier-utils/articles/multi-cohort.html) |
| Assess uncertainty, stability and prediction | [Understand uncertainty](https://mdmanurung.github.io/flashier-utils/articles/uncertainty.html); [Check factor stability](https://mdmanurung.github.io/flashier-utils/articles/stability-replicability.html); [Evaluate held-out prediction](https://mdmanurung.github.io/flashier-utils/articles/holdout.html) |
| Visualize single-cell factors | [Display factors in Seurat](https://mdmanurung.github.io/flashier-utils/articles/single-cell.html) |
| Find sources and cite methods | [Original authors and sources](https://mdmanurung.github.io/flashier-utils/articles/credits.html); [Find sources and cite methods](https://mdmanurung.github.io/flashier-utils/articles/provenance.html) |

Browse all groups in the [tutorial index](https://mdmanurung.github.io/flashier-utils/articles/index.html).
See the [function reference](https://mdmanurung.github.io/flashier-utils/reference/index.html)
for arguments and return values.

Group summaries are descriptive. Association tests condition on the fitted
factors and require you to identify the independent biological units, such as
donors. Factor posterior uncertainty does not replace a regression standard error.

### Where each tutorial's data come from

Every tutorial has a "Data and objects" section or an opening paragraph that says
which objects it uses and where they come from. In short:

| Data | Tutorials | How to get it |
| --- | --- | --- |
| GTEx association z-scores shipped with flashier | [Fit and interpret](https://mdmanurung.github.io/flashier-utils/articles/first-fit.html); [Scales](https://mdmanurung.github.io/flashier-utils/articles/representation-units.html); [Edit fits](https://mdmanurung.github.io/flashier-utils/articles/edit-fits.html); [Backprojection](https://mdmanurung.github.io/flashier-utils/articles/backprojection.html); [Uncertainty](https://mdmanurung.github.io/flashier-utils/articles/uncertainty.html); [Stability](https://mdmanurung.github.io/flashier-utils/articles/stability-replicability.html) | `data("gtex", package = "flashier")`; nothing to download |
| Published human-pancreas programs (Jamboree CEL-Seq2) | [Pancreas programs](https://mdmanurung.github.io/flashier-utils/articles/pancreas-programs.html); [Defining features](https://mdmanurung.github.io/flashier-utils/articles/largest-distinctive.html); [Annotations](https://mdmanurung.github.io/flashier-utils/articles/program-annotation.html); [Diagnostics](https://mdmanurung.github.io/flashier-utils/articles/program-diagnostics.html); [Transfer](https://mdmanurung.github.io/flashier-utils/articles/dataset-transfer.html); [Pathways](https://mdmanurung.github.io/flashier-utils/articles/enrichment.html) | Run `.github/scripts/prepare-pancreas.R` once; it writes `pancreas-worked-example.rds` |
| Pancreas cells from four technologies (Jamboree) | [Contrasts](https://mdmanurung.github.io/flashier-utils/articles/contrasts.html); [Held-out prediction](https://mdmanurung.github.io/flashier-utils/articles/holdout.html); [Conditions](https://mdmanurung.github.io/flashier-utils/articles/multi-cohort.html); [Seurat](https://mdmanurung.github.io/flashier-utils/articles/single-cell.html); [Pathways](https://mdmanurung.github.io/flashier-utils/articles/enrichment.html); [Transfer](https://mdmanurung.github.io/flashier-utils/articles/dataset-transfer.html) | The same script also writes `pancreas-cells.rds` |
| Simulated | [Protein interpretation](https://mdmanurung.github.io/flashier-utils/articles/protein-interpretation.html) | Created inside the article |

The pancreas data are not shipped with the package. From a clone of this repository,
run the following once (it downloads about 86 MB), then set the same variable
before running the tutorial code:

```sh
FLASHIER_UTILS_PANCREAS_DIR=~/pancreas-jamboree-example Rscript .github/scripts/prepare-pancreas.R
```

## Original authors and citations

flashier is by Jason Willwerscheid, Peter Carbonetto, Wei Wang and Matthew
Stephens.

ebnm is by Jason Willwerscheid, Matthew Stephens and Peter Carbonetto.
The adapted singlecelljamboreeR feature helpers are by Peter Carbonetto and
Matthew Stephens. Plotting and reconstruction workflows credit David Zemmour
and the Zemmour Lab.

Cite the methods used in your analysis with `citation("flashier")`,
`citation("ebnm")` and the relevant optional packages.
[Original authors and sources](CREDITS.md) lists the full credits and licenses.
Mikhael Manurung maintains flashier.utils.
