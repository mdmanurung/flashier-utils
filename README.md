# flashier utils

An R toolkit for interpreting empirical Bayes matrix factorizations using public
flashier interfaces. R package identifier: `flashier.utils`.

Core methods are by Jason Willwerscheid, Peter Carbonetto, Wei Wang and Matthew
Stephens (flashier), and Jason Willwerscheid, Matthew Stephens and Peter Carbonetto
(ebnm), with their package contributors. Adapted singlecelljamboreeR feature
helpers are by Peter Carbonetto and Matthew Stephens. Full method, workflow and
dependency author lists are in [Original authors and sources](CREDITS.md).

```r
install.packages("remotes")
remotes::install_github("mdmanurung/flashier-utils")
```

Documentation: [reference and tutorials](https://mdmanurung.github.io/flashier-utils/).
Maintainer: Mikhael Manurung <mikhael.manurung@gmail.com>.

Source version: **1.0.0**, with 34 exports. Signed decoupleR
ULM/MLM networks annotate estimated feature programs; fgsea accepts named sets.
Network scores and association tests condition on the fitted factor basis.
SuSiE, brms/uncertain-response, NNLM/cycling, CorShrink dosym/bootstrap and
frozen-noise projection remain unsupported.

[Release 1.0.0](https://github.com/mdmanurung/flashier-utils/releases/tag/v1.0.0)
is tagged at `794ccfd`. Its exact source passed 380 assertions, a clean R4_51
R CMD check (34 examples, 12 vignettes), and isolated core installation.
[All five CI jobs](https://github.com/mdmanurung/flashier-utils/actions/runs/37432844465)
passed: Linux release/devel, macOS, Windows and full engines, including optional
dependency loading. The published source tarball checksum and installation were
verified. R >= 4.3 is declared; the oldest supported R has not been tested.

Software verification and scientific calibration are separate. Expanded
known-activity calibration completed 10,000/10,000 runs without warnings or
failures, with all 200 native audits passing at tolerance 1e-8. F1 was the primary
target: 500 seeds (2101–2600) per scenario, 4/20 factors, null/planted effect of
one, homogeneous/heterogeneous residuals for lm/limma, and 20 donors with two
visits for dream (Satterthwaite df, serial workers).

| Engine | Successful / attempted | 95% interval coverage | Null rejection | Power |
|---|---:|---:|---:|---:|
| lm | 4,000 / 4,000 | 95.4–95.6% | 4.4–4.6% | 50.8–86.4% |
| limma | 4,000 / 4,000 | 93.6–95.6% | 4.4–6.4% | 53.6–87.4% |
| dream | 2,000 / 2,000 | 95.2% | 4.8% | 85.2–86.6% |

Rates are ranges across separate 500-seed scenarios, without pooling correlated
factors. Dream coverage was 476/500 in each scenario (Monte Carlo SE 0.96
percentage points; exact 95% binomial CI 92.94–96.90%). Limma's lowest coverage
was 468/500 under heterogeneous noise (MCSE 1.09 points; CI 91.08–95.58%).
The earlier dream finding, 21/25 seeds (**84%**, MCSE 7.33 points), is preserved
and reproduced by the first 25 seeds. All expanded scenario coverage CIs include
95%; these simulations assess these known-activity settings. They do not
establish learned-factor calibration or account for uncertainty in an EBMF fit.

## Getting started

Interpret signed empirical Bayes matrix factorizations using matched sample
activities and feature programs. No preprocessing is performed. Biological
annotation uses explicitly supplied feature sets or networks. Feature units are
the units of the supplied matrix.

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
devtools::test()
rcmdcheck::rcmdcheck(".", args = "--no-manual")
```


## Signed network annotation

Annotate EBMF feature weights with signed decoupleR networks (ULM by default):

```r
# network is a data frame with source, target (feature ID), and mor (signed weight).
# For gene-based fits, e.g. network <- decoupleR::get_collectri(organism = "human")
annotation <- enrich_factors(view, network, engine = "decoupleR",
                             control = list(method = "ulm", minsize = 5L))
annotation$table  # source, condition (factor ID), score, native p_value, basis_id
annotation$analysis_metadata$pathway_sizes  # input and overlapping target counts
```

Use `method = "mlm"` for multivariate network scoring. Targets must match the
fit's feature IDs; `feature_map` provides an explicit one-to-one mapping.
Scores describe estimated factor programs and depend on factor orientation.
The complete mapped feature universe is retained; network hashes, excluded
targets, mapping losses and annotation metadata are recorded. Duplicate edges
and nonfinite weights are rejected. Native filtering and errors are retained;
centering across factors defaults to FALSE. Native p-values do not propagate
EBMF uncertainty or establish biological replication. Supply custom marker networks for non-gene features such as CyTOF
markers; gene regulatory networks require gene-level features.

## Original authors and sources

This toolkit wraps and builds on methods by their original authors. Every function
manual names its sources; `provenance()` provides audited authors, source pins
and reuse classifications. See [full credits](CREDITS.md), the
[credits article](https://mdmanurung.github.io/flashier-utils/articles/credits.html),
and the installed `NOTICE` and license texts.
