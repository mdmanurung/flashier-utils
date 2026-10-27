# rOpenSci-style review of flashier.utils, 2026-10-27

AI-assisted review passes run with the rOpenSci review tools on commit `22c8877`
(`main`). Prompt, skill text, method, model and what was not run are recorded in
each report's transparency header, following the shared template of
`ropensci-review-tools/ropensci-skills` (commit `9dccd2e1`).

| Report | Pass | Verdict |
|---|---|---|
| `Review_AI_tools.qmd` | pkgcheck, goodpractice, autotest, srr, package-standards gates | amber: 5 gates fail, none in R code |
| `Review_AI_test_audit.qmd` | run-test-audit | green/amber: 88.3 percent coverage, all exports tested, depth gaps in 2 files |
| `Review_AI_complexity.qmd` | run-complexity | amber: 5 functions over cc 95, no duplication |
| `Review_AI_dependencies.qmd` | run-dependency-review | green: nothing prunable; flashier pin is the issue |
| `Review_AI_performance.qmd` | run-performance-review | green: cost is native; guards in place |

`evidence/` holds the derived tables each report cites (complexity per function,
coverage by file, dependency tree, goodpractice results and lint counts, pkgcheck
checks, pkgstats percentiles, local R CMD check output, autotest results,
sessionInfo). Large tool objects stay in the local scratch area and are not
committed.

No package files were modified by the review. The five gates in `Review_AI_tools.qmd`
were fixed in the following commit (examples, package doc, CONTRIBUTING, oldrel CI,
quoted 'flashier'); the reports describe the state at `22c8877`. Fixes are listed per report under
"Suggested fix priority" and are applied in separate commits, if at all.
