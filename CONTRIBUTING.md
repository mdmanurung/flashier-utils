# Contributing to flashier.utils

Thanks for considering a contribution. Bug reports, documentation fixes and
scientific questions are all welcome.

## Report a bug or ask a question

Open an issue at <https://github.com/mdmanurung/flashier-utils/issues>. Include
the package version (`packageVersion("flashier.utils")`), the `flashier` and
`ebnm` versions, a small reproducible example with named matrices, and the
output of `provenance(result = <your result>)` if the question concerns a
function's origin or runtime versions.

## Propose a change

1. Open an issue first for anything beyond a typo, so the scope and the
   scientific contract can be agreed before code is written.
2. Fork, branch from `main`, and keep changes focused.
3. Keep the contracts: unique row and column IDs, explicit `sample_side`,
   factor selections as character IDs, no silent normalisation, uncertainty
   labelled with what it conditions on. Do not add single-value or reserved
   arguments to exports.
4. Credit original authors. Any method adapted from another package needs its
   authors in the roxygen `@author` field, a row in `inst/provenance.csv`, an
   entry in `CREDITS.md`, and the upstream license notice in `inst/NOTICE`
   when source is reused. Do not invent individual authors where attribution
   is unresolved.
5. Add or update tests in `tests/testthat/` and run the suite:

   ```r
   devtools::document()
   devtools::test()
   devtools::check()
   ```

   Keep every exported function with `@param`, `@return`, `@examples` and
   `@author`; the site checker enforces attribution and link integrity.
6. Add a NEWS.md bullet for user-visible changes.
7. Open a pull request against `main`. CI runs R CMD check on Linux (release,
   devel, oldrel), macOS and Windows plus a job with every optional engine.

## Code style

Base R first; `ggplot2` for plots; dependencies only in `Suggests` and gated with
`.check_engine()` unless the whole package needs them. Lines up to 120
characters; run `styler::style_pkg()` on files you touch.

## Conduct

Be respectful. Reviews and discussions follow the
[rOpenSci code of conduct](https://ropensci.org/code-of-conduct/).
