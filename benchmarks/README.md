# Benchmarks and experiment scripts

Ad-hoc experiment, baseline, and smoke-test scripts that exercise this package.
These are kept separate from the package source (`R/`, `tests/`, `man/`) so the
package itself stays clean; nothing here is run by `R CMD check`.

| Script | What it does |
|---|---|
| `fu_baseline.R` | Loads the package via `pkgload::load_all()` and runs the full testthat suite file by file, tabulating passed/failed/error/skipped per file. Used as the pre/post baseline when swapping the backend or dependency versions (e.g. flashier 1.0.60 vs 1.0.7). |
| `sync3_smoke.R` | SYNC-3 release-gate smoke test: refits one Stephenson cell type (`CD4_naive`, K30 masked pseudobulk) through the migrated flashier.utils backend and checks parity with the committed masked fit (n factors, ELBO, max abs diff of `L_pm`/`F_pm`/scores, fallback flag). |

## Notes

- Scripts are provenance records of what was actually run, including
  environment-specific absolute paths (`/workspace/.Rlib`, sibling repo
  checkouts, cached `.rds` inputs). Update paths before rerunning elsewhere.
- `sync3_smoke.R` also loads the cellprograms package from a sibling checkout;
  the companion benchmark harness used in that sync effort lives in
  `mdmanurung/cellprograms` under `benchmarks/cmp_bench/`.
