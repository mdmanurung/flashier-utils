# CREDITS and NOTICE review, 2026-11-03

Every package the code calls is declared and credited, and every pinned or audited
upstream commit is named in CREDITS.md or inst/NOTICE. Rerun with:

```sh
scripts/run-r.sh reports/credits-review-20261103/check-credits.R reports/credits-review-20261103
```

## What the script checks

| Table | Rows | Check |
|---|---|---|
| `packages.tsv` | every Imports, Suggests, `pkg::` call and `.check_engine()` name | declared in DESCRIPTION, called in `R/`, credited in CREDITS.md |
| `provenance.tsv` | 63 rows of `inst/provenance.csv` | function defined in `R/`; source commit named in CREDITS.md or NOTICE |
| `pins.tsv` | 4 `Remotes:` pins | commit named in NOTICE and CREDITS.md |

## Findings and fixes

The first run found four gaps. Each is fixed in the same commit as this receipt.

- `scales` became an Import on 2026-11-03 and had no credit entry. Added
  `### scales (1.4.0)` from its installed DESCRIPTION and citation.
- The CorShrink pin `a9f6ba0` was named in neither file. Added to the CorShrink
  entry in CREDITS.md and as a call-only paragraph in NOTICE (GPL (>= 2), no
  source copied).
- `irlba` and `softImpute` are user-selectable `fit_multistart_ebmf()`
  initializers, reached through `flashier::flash_greedy_init_*()`, with no credit.
  Added entries that name the route.
- `patchwork` was in Suggests and credited, but no code, test or vignette uses
  it. Removed from Suggests, NOTICE and CREDITS.md.

After the fixes all three tables are clean. `irlba` and `softImpute` show
`called_in_R = FALSE` because `R/` calls them through flashier; that is expected.
`knitr`, `rmarkdown` and `testthat` are build tools and carry no credit entry.
The bundled licence file `inst/licenses/singlecelljamboreeR.txt` covers the only
copied code (modified `rank_effects` and `compute_le_effects`), as NOTICE states.

## Not checked

- Author lists inside existing CREDITS.md entries against current upstream
  DESCRIPTION files; entries record the documentation snapshot versions.
- Whether the audited workflow repositories (ZemmourLib, ImmGenT-GP, STARCAT)
  changed their licences after their recorded commits.
