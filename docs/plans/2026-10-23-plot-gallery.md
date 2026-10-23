# Plot gallery plan

Date: 2026-10-23
Project: /exports/para-lipg-hpc/mdmanurung/flashier-utils
Baseline: main @ dcb5c471f1624035f6e11e095927f10fd544f7db
Status: not started; planning only

## Objective

Add a visual catalogue of supported plotting calls from flashier and flashier.utils.
Use the R Graph Gallery density page's library sections, preview cards and linked examples.
Generate original figures and recipes; do not copy its assets or page source.

Assumption: add `vignettes/plot-gallery.Rmd` as a fifth vignette, alongside the four
existing tutorials. The user's optional choice between a fifth vignette and a
website-only article was requested; no answer was received before this draft.
This plan assumes a fifth vignette because the user requested a vignette page.
If the user chooses website-only, use `vignettes/articles/plot-gallery.Rmd`, exclude
that directory from package builds, and keep four installed vignettes. Do not remove
a tutorial or silently relax all article limits.

## Page structure

1. Title: Plot gallery. One sentence inviting the reader to choose a preview.
2. Compact links to Native flashier, flashier.utils, Shared setup and Credits.
3. Native flashier: six preview cards, ordered from fit overview to sample structure.
4. flashier.utils: six preview cards, ordered from activities to feature effects.
5. Shared setup and twelve anchored R recipes. Each card links to its recipe.
6. Credits and links to current function references and the four tutorials.

Each card contains a real preview, a short title, one use-case sentence and a
package label. Make the image/title a keyboard-accessible link to its recipe.
Recipes show the minimal plotting call and its value/normalization meaning.
Link to function help for full argument lists. Do not add one article per card.

Use a scoped CSS grid: three columns on desktop, two on tablets and one on phones.
Keep figures at a consistent landscape ratio. Use existing theme colors, readable
labels and descriptive image alt text. Scoped CSS must also work in the installed
html_vignette, which does not share pkgdown's Bootstrap 5 layout.
No JavaScript filters, modal viewers, new framework or plotting API is needed.

## Verified coverage

The installed flashier 1.0.60 exports these six dedicated plotting functions.
The plot.flash method dispatches to them with plot_type; explain that alias once,
rather than adding duplicate cards or using deprecated include_scree/include_pm.

| Card | Supported call | What the preview shows |
|---|---|---|
| Factor variance | flashier::flash_plot_scree() | Native factor PVE; sorted labels preserve factor identity. |
| Native heatmap | flashier::flash_plot_heatmap() | Loadings across named samples. |
| Native bars | flashier::flash_plot_bar() | Individual posterior means for selected components. |
| Native histogram | flashier::flash_plot_histogram() | Loading distributions, with explicit grouping. |
| Native scatter | flashier::flash_plot_scatter() | Loadings against a supplied, aligned numeric covariate. |
| Native structure | flashier::flash_plot_structure() | Nonnegative, max-absolute-scaled loadings stacked by component. |
| Activity heatmap | plot_factor_activity(style = "heatmap") | Activity in declared representation units. |
| Activity points | plot_factor_activity(style = "point") | Grouped sample activity, faceted by factor. |
| Activity boxplots | plot_factor_activity(style = "distribution") | Boxplots of activity by group; not density or violin plots. |
| Activity structure | plot_factor_activity(style = "structure") | Explicit row-sum-normalized nonnegative activity. |
| Feature heatmap | plot_factor_features(style = "heatmap") | Selected signed feature effects. |
| Feature bars | plot_factor_features(style = "bar") | Selected feature effects with their sign. |

Keep the current 21-function public API. Do not use triple-colon calls or re-export
internal plotters. Optional internal diagnostics/network/gate/embedding plots are
outside this public gallery.

## Shared example and scientific contracts

Use one small, seeded synthetic matrix with named sample rows and feature columns,
at least two recovered factors, nonnegative loadings and signed feature effects.
Fit directly with flashier::flash and explicit appropriate priors; declare
sample_side = "rows", feature_scale and prior_support in standardize_factors.
Reuse the existing vignette's simulation conventions; do not add a fitting helper.
Provide sample/group metadata and a numeric plotting covariate with exact ID alignment.
Use only existing dependencies and local computation; no downloads or scheduler jobs.

Fit once and reuse the same fit/view across all cards. Keep raw activity/effect
coordinates unchanged. Use named arguments under the 2.0 API. Label simulated data.
Native posterior-mean displays apply max-absolute scaling; utilities use resolved
units unless an explicit display normalization is requested. Native stacked values
and utilities row-sum proportions are different quantities and must be labeled.
Structure examples require nonnegative activity and must disclose any zero-row
omission. Native PVE is not an additive PCA partition for nonorthogonal factors.
Keep factor identities visible when PVE ordering changes.

Store plotting calls once in named knitr chunks. Use ref.label with eval = FALSE
for the corresponding displayed recipe instead of duplicating or rerunning calls.
Do not use eval(parse(...)). Render thumbnail images from the actual plot outputs.
Keep file/chunk names explicit and check their paths in both rendering formats.
A clearly marked shared setup supplies the fit/view/metadata for each short recipe.

## Review findings and resolutions

| Finding | Resolution |
|---|---|
| Four tutorials are currently enforced as a four-article total. | Preserve them; allow only the named gallery in addition, subject to the placement choice above. |
| Only two plotting functions remain public in utilities. | Cover their six supported styles; omit retained internal APIs. |
| Native and utilities structure plots have different scales. | Label normalization on each card/recipe; do not present them as interchangeable. |
| Utilities distribution style is a boxplot. | Label the card as boxplots; no invented density option. |
| Bootstrap-only cards would fail in installed vignette HTML. | Use scoped CSS grid that works in both outputs. |
| Repeated fits or duplicated snippets could drift. | One shared fit and named chunks; recipes reuse chunk source. |
| Native mean-side names differ from sample/feature semantics. | Explicit pm_which and sample_side; align IDs before grouping/covariate use. |
| Original plotting authors must remain visible. | Cite package methods and existing source credits; link CREDITS and use package citations. |
| docs/ contains ignored generated site output. | Save the plan here as requested by the planning workflow; render implementation previews to validation/, never clean docs/ while this plan is active. |

## Progress summary

- [ ] G1: Prepare and check one shared plotting example — not started.
- [ ] G2: Build the twelve-card vignette and recipes — not started.
- [ ] G3: Connect navigation and update article checks — not started.
- [ ] G4: Render, inspect and qualify both outputs — not started.

## Implementation steps

### [ ] G1: Prepare and check the shared example

- Status: not started. Dependencies: confirm the placement assumption before creating the article.
- Outcome/actions: put the minimal seeded fit/view/metadata setup inside the chosen Rmd; verify actual native arguments against installed help and print all twelve ggplot objects.
- Affected areas: new gallery Rmd only; no package implementation, exports or dependencies.
- Validation: R4_51 through scripts/run-r.sh; assert rank >= 2, names/metadata alignment and valid nonnegative activity for structure plots. Check generated plot data, not empty placeholders.
- Evidence: pending; retain runtime/version, warnings and rendered example outputs under validation/plot-gallery-20261023/.

### [ ] G2: Build the gallery page

- Status: not started. Dependencies: G1.
- Outcome/actions: add the two card grids, anchored recipes, reusable chunks, shared setup, scale labels and credits. Use R Markdown/knitr and minimal scoped CSS.
- Affected areas: the chosen Rmd; existing article links only where needed.
- Validation: twelve unique previews and twelve working recipe targets; visible package/function names; working function references; alt text; no signed structure input or hidden internal calls. Confirm plotting leaves fit/view coordinates unchanged.
- Evidence: pending; HTML and figure assets in the validation preview directory.

### [ ] G3: Connect navigation and article validation

- Status: not started. Dependencies: G2 and placement decision.
- Outcome/actions: add Plot gallery to _pkgdown.yml as a distinct gallery entry, preserve the four tutorial entries, and add a concise README link.
- Affected areas: _pkgdown.yml, README.md, scripts/check-site.py, scripts/test-check-site.py; .Rbuildignore only for the website-only alternative.
- Validation: count at most four narrative tutorials and allow the specifically named gallery; never change the guard to an unrestricted higher number. Add small fixtures for four tutorials plus gallery passing, an extra narrative tutorial failing, and a broken gallery image/recipe link failing. Inventory nested article sources if the website-only alternative is selected.
- Evidence: pending; checker fixture output and new site receipt. The scripts are local helpers currently excluded from Git; decide tracked CI integration only if delivery is requested.

### [ ] G4: Render and qualify

- Status: not started. Dependencies: G3.
- Outcome/actions: render standalone vignette HTML and a pkgdown preview using scripts/run-r.sh; keep outputs under validation/plot-gallery-20261023/. Use existing supported offline asset/metadata handling only if required.
- Validation: inspect desktop/tablet/mobile layouts; validate all images, local links, recipe anchors and attribution. Build/check the package with vignettes and cached pancreas input. Confirm five installed vignette outputs for the fifth-vignette choice, or four installed outputs plus the website gallery for the alternative. Retain the 21 exports and all existing notices/registry records. Do not rerun unrelated scientific studies.
- Evidence: pending; preview HTML, relevant screenshots when available, check logs, source/asset hashes and final acceptance receipt. Report manual-only visual checks if browser screenshots are unavailable.

## Credits and sources

- Layout reference: https://r-graph-gallery.com/density-plot.html — page organization only.
- Native API: installed flashier 1.0.60 NAMESPACE and help for flash_plot_* and plot.flash; inspection output /tmp/flashier-gallery-api.txt is temporary, so re-query installed help when implementing.
- Utilities API: R/plots.R and man/plot_factor_activity.Rd / man/plot_factor_features.Rd at the baseline commit.
- Authorship: retain FLASH/flashier and plotting-method credits through citation("flashier"), provenance(), existing manuals, CREDITS.md, inst/NOTICE and inst/provenance.csv. Do not copy upstream source or invent per-function authors.

## Scope limits and separate findings

This is a documentation plan. No new plotter, export, scientific filtering rule,
public release, scheduler run, commit or push is part of planning.
The preceding core review's approved implementation is complete. Non-lm offsets
remain explicitly unsupported; full convergence/calibration and wider upstream-version
qualification remain separate evidence work. Cross-platform hosted checks already passed.

The linked ImmGenT cell filter is separate: it repeatedly removes samples above
max-scaled total membership 10, rescales and runs twelve passes, then compensates
feature effects on the retained cells. Our internal membership score is one-pass;
the iterative removal/rescaling workflow is not implemented. Do not add it to this gallery.
Source: https://raw.githubusercontent.com/AgueroZZ/immgenT-GP-analysis/main/analysis/Methods_FlashierFit.Rmd

## Decisions and evidence log

- 2026-10-23 — Checked the reference page, installed native plotting exports/arguments/help, current utilities styles, navigation and article guard. Reviewed and patched this plan; implementation has not started.
- 2026-10-23 — Placement preference requested asynchronously; fifth-vignette assumption remains explicit pending any reply.
