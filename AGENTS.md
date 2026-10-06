# Repository instructions

Use the R4_51 runtime through scripts/run-r.sh for local R work.
The R package lives at the repository root. Keep the existing public API and
scientific contracts unless the user requests a change.

Always credit original function and package authors in source documentation,
manuals, tutorials and website content. Preserve inst/NOTICE, all upstream
license notices, and inst/provenance.csv. Cite upstream methods and distinguish
wrappers, modified source and conceptual workflow adaptations. Do not invent
individual authors where source attribution is unresolved. Keep CREDITS.md and
vignettes/credits.Rmd consistent when attribution changes.

Scientific benchmarks, frozen outputs and validation receipts stay in the local
archive and outside package builds. Do not rewrite published release tags.
