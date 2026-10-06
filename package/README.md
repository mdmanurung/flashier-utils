# flashier utils

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("mdmanurung/flashier-utils", subdir = "package")
```

Interpret signed empirical Bayes matrix factorizations using matched sample
activities and feature programs. No preprocessing or biological annotation is
performed. Feature units are the units of the supplied matrix.

```r
library(flashier.utils)
set.seed(3)
X <- tcrossprod(matrix(rnorm(24), 12, 2), matrix(rnorm(16), 8, 2)) +
  matrix(rnorm(96, sd = 0.2), 12, 8)
dimnames(X) <- list(paste0("s", 1:12), paste0("g", 1:8))
fit <- fit_ebmf(X, sample_side = "rows", seed = 1, max_factors = 2,
                feature_scale = "arcsinh_intensity", verbose = 0)
view <- standardize_factors(fit)
factor_features(view, factor = "F1", select = "both", n = 3)
plot_factor_activity(view)
beta <- setNames(rep(1, ncol(view$activity)), colnames(view$activity))
backproject_contrast(view, beta, input = "coefficients",
                     factor_basis = view$manifest)
provenance("standardize_factors")
```

Activity and effects satisfy `tcrossprod(A, B) == fitted(fit)` with samples in
rows. Native LDF D is placed once. A coefficients' activity manifest is required
for detached backprojection inputs; altered scale/sign/order requires matching
coefficients and a matching manifest. Point coefficients produce point feature
effects without p-values. Posterior factor SD describes conditional variational
uncertainty, not a regression standard error.

Signed heatmaps are the default. Structure proportions require nonnegative
scores and explicit display normalization. Distinctiveness depends on the
recorded sign and scale convention; it is a descriptive competitor margin.

Code origins are available through `provenance()` and `inst/NOTICE`. The two
MIT singlecelljamboreeR feature helpers are locally adapted, so enrichment/Stan/
Seurat are unnecessary for this core workflow. Maintainer: Mikhael Manurung
(mikhael.manurung@gmail.com). GitHub project name: flashier-utils; R package
identifier: flashier.utils (R package names cannot contain spaces or hyphens).

Development checks from the repository root:

```r
devtools::test("package")
rcmdcheck::rcmdcheck("package", args = "--no-manual")
```

## Later modules

The 34 exports cover explicit lm/limma/dream association, native fgsea and
decoupleR enrichment, donor bootstrap/checkpoints, named program assignment,
new-sample projection and masked scoring, native mash/correlation, and exact-ID
Seurat/uwot displays. Optional dependencies load only when requested.

All downstream inference conditions on estimated factor scores unless stated
otherwise. Repeated-donor benchmark coverage was below nominal; calibration is
not established. brms/uncertain-response and SuSiE modes are deferred. CorShrink
supports its audited default only; native dosym=TRUE fails and is refused.

Source notices are installed with the package. Executable vignettes describe
assumptions, units, donor replication and display-only transforms. Cross-platform
validation results are available in GitHub Actions. The scientific development evidence
archive is retained separately from this source distribution.
