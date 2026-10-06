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

Software verification and scientific calibration are separate. The earlier
Linux R 4.5.1 baseline passed 343 assertions, 34 examples and 12 vignettes with
zero R CMD check errors, warnings or notes. The published source snapshot
`c35f6b9` passed [Linux release/devel, macOS, Windows and full-engine CI](https://github.com/mdmanurung/flashier-utils/actions/runs/37423198808),
including optional-dependency loading. The signed-network candidate passes 380
assertions under R4_51 and installs with only core dependencies. These receipts
describe their respective snapshots;
each release requires successful checks of its exact source commit. R >= 4.3
is the declared minimum; the oldest supported R has not been tested.

Known-activity repeated-donor dream intervals covered 21/25 seeds (**84%**,
Monte Carlo SE 7.33%) under both null and planted effects. That small simulation
does not establish nominal coverage or learned-factor calibration. Expanded
500-seed calibration is scheduled after release; it will preserve the earlier
results and report failures and Monte Carlo uncertainty. Release notes and
tagged artifacts must identify the exact checked source commit.
