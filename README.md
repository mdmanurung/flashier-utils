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

## Migrating to the 2.0 development API

This development version exposes 23 public functions. The remaining 34 functions
are internal; their implementations, authorship and provenance records remain.
Calls to removed exports must migrate to the core or the upstream package, or
stay on the previous version. Internal functions are not supported APIs.

Nine core functions also change positional argument order. Use named scientific
arguments, such as `factor_features(view, factor = "F1", select = "largest", n = 5)`.
Defaults and scientific calculations are unchanged. See [NEWS](https://mdmanurung.github.io/flashier-utils/news/index.html) for the
removed functions and affected signatures, and the function references for current
argument order.


## Quick start

Use the GTEx association z-scores bundled with flashier. Tissues occupy columns;
SNP–gene pairs occupy rows. A fit created directly with flashier works through
one explicit standardization step:

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

Activities have samples in rows; effects have features in rows, with matching
factor columns. Supply unique row/column IDs and declare sample orientation and
measurement units. The package does not preprocess your measurements.
Keep the native fit for uncertainty; most downstream utilities consume `view`.
Unnamed fits are rejected. Native intercepts without role metadata remain
ordinary factors. Compatibility is tested against the pinned flashier version,
not every historical or future object format.

## Four focused tutorials

| Task | Tutorial |
| --- | --- |
| Beginner: native fit, extraction, features, one plot | [Get started](https://mdmanurung.github.io/flashier-utils/articles/first-fit.html) |
| Coordinates, conditional uncertainty and group contrasts | [Advanced use](https://mdmanurung.github.io/flashier-utils/articles/advanced-use.html) |
| Different priors: prediction, program agreement and known truth | [Compare priors](https://mdmanurung.github.io/flashier-utils/articles/compare-priors.html) |
| Real-data interpretation and reference transfer | [Pancreas programs](https://mdmanurung.github.io/flashier-utils/articles/pancreas-programs.html) |

The first uses bundled GTEx data; the middle two use explicit simulations.
The pancreas example uses attributed published Jamboree data, prepared once
outside package builds (about 86 MB downloaded):

```sh
FLASHIER_UTILS_PANCREAS_DIR=~/pancreas-example scripts/run-r.sh .github/scripts/prepare-pancreas.R
```

Set the same variable before running the case study. No data are downloaded by
package builds. See the [function reference](https://mdmanurung.github.io/flashier-utils/reference/index.html)
for complete arguments. Group summaries are descriptive; association tests require
explicit independent units. Conditional posterior uncertainty does not replace
regression standard errors or biological replication.

[Plot gallery](https://mdmanurung.github.io/flashier-utils/articles/plot-gallery.html)
shows twelve plotting recipes from one simulated fit.

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
