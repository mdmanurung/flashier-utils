# flashier utils

An R toolkit for interpreting empirical Bayes matrix factorizations using public
flashier interfaces. R package identifier: `flashier.utils`.

```r
install.packages("remotes")
remotes::install_github("mdmanurung/flashier-utils", subdir = "package")
```

See [package/README.md](package/README.md) for examples and assumptions.
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
