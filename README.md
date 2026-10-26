# flashier utils

Interpret a [flashier](https://willwerscheid.github.io/flashier/) matrix
factorization. After one call you get a *view*: a named samples × factors
activity matrix and a named features × factors effects matrix that share factor
columns and reconstruct the fit. From a view you can rank the features that
drive each factor, plot activity across samples, summarize and test factors
against sample metadata, back-project a group contrast into feature units, and
check posterior uncertainty, stability and held-out prediction.

## Installation

Requires R 4.3 or later. The install pulls the pinned flashier commit.

```r
install.packages("remotes")
remotes::install_github("mdmanurung/flashier-utils")
```

Optional packages unlock some functions: `psych` and `clue` for
`factor_replicability()` and `factor_stability()`; `limma` or `variancePartition`
for the matching `test_factors()` engines; `gbcd` for `fit_gbcd_programs()`;
`fastTopics` for `test_gom_de()`. Each page in the function reference names what
it needs.

## Quick start

Use the GTEx association z-scores bundled with flashier. Tissues occupy columns;
SNP–gene pairs occupy rows.

```r
library(flashier.utils)
data("gtex", package = "flashier")
set.seed(1)
fit <- flashier::flash(gtex, greedy_Kmax = 3, backfit = TRUE, verbose = 0)
view <- standardize_factors(fit, sample_side = "columns",
                            feature_scale = "association_z_score")
plot_factor_activity(view)
factor_features(view, factor = "F1", select = "largest", n = 5)
```

Your matrix needs unique row and column names. Tell `standardize_factors()`
which side holds samples and what the values measure; the package records both
and never preprocesses or rescales your data. Keep the native `fit` as well as
the `view`: uncertainty functions need the fit, everything else takes the view.

## Workflow

| Step | Functions |
| --- | --- |
| Fit, or standardize an existing fit | `fit_ebmf()`, `standardize_factors()`, `factor_activity()`, `factor_effects()`, `factor_pve()` |
| Interpret factors | `factor_features()`, `plot_factor_features()`, `plot_factor_activity()` |
| Relate factors to metadata | `summarize_factor_activity()`, `compare_factor_groups()`, `test_factors()`, `backproject_contrast()` |
| Check uncertainty and reproducibility | `factor_uncertainty()`, `factor_intervals()`, `factor_stability()`, `factor_replicability()`, `factor_holdout()` |
| Transfer programs to new data | `project_factor_activity()`, `project_factor_features()` |

## Four focused tutorials

| Task | Tutorial |
| --- | --- |
| Beginner: native fit, extraction, features, one plot | [Get started](https://mdmanurung.github.io/flashier-utils/articles/first-fit.html) |
| Coordinates, uncertainty, group contrasts, nonnegative programs | [Advanced use](https://mdmanurung.github.io/flashier-utils/articles/advanced-use.html) |
| Different priors: prediction, program agreement and known truth | [Compare priors](https://mdmanurung.github.io/flashier-utils/articles/compare-priors.html) |
| Real-data interpretation and reference transfer | [Pancreas programs](https://mdmanurung.github.io/flashier-utils/articles/pancreas-programs.html) |

The first uses bundled GTEx data; the middle two use explicit simulations.
The pancreas example uses attributed published Jamboree data, prepared once
outside package builds (about 86 MB downloaded):

```sh
FLASHIER_UTILS_PANCREAS_DIR=~/pancreas-example scripts/run-r.sh .github/scripts/prepare-pancreas.R
```

Set the same variable before running the case study. No data are downloaded by
package builds. The [plot gallery](https://mdmanurung.github.io/flashier-utils/articles/plot-gallery.html)
shows twelve plotting recipes from one simulated fit, and the
[function reference](https://mdmanurung.github.io/flashier-utils/reference/index.html)
lists complete arguments.

## Limits

Unnamed fits are rejected. Native intercepts without role metadata remain
ordinary factors. Group summaries are descriptive; association tests require
explicit independent units. Conditional posterior uncertainty does not replace
regression standard errors or biological replication. Compatibility is tested
against the pinned flashier version, not every historical or future object
format.

## Migrating from 1.x

Version 2.0 exposes 24 public functions and makes the rest internal. Nine core
functions changed positional argument order, so use named arguments, such as
`factor_features(view, factor = "F1", select = "largest", n = 5)`. Defaults and
calculations are unchanged; [NEWS](https://mdmanurung.github.io/flashier-utils/news/index.html)
lists the removed exports and affected signatures.

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
