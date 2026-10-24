# Plot gallery plan

Date: 2026-10-23
Project: /exports/para-lipg-hpc/mdmanurung/flashier-utils
Baseline: main @ 09bca1c; package implementation @ 70a30d6
Last updated: 2026-10-24
Status: complete and delivered; source 952ebed pushed, hosted checks passed, gallery live
Estimated implementation time: 4–6 hours, excluding hosted CI queue time.

## Objective

Add a visual catalogue of supported plotting calls from flashier and flashier.utils.
Use the R Graph Gallery density page's library sections, preview cards and linked examples.
Generate original figures and recipes; do not copy its assets or page source.

Decision for this plan: use `vignettes/articles/plot-gallery.Rmd` as a
website-only article. Keep the four installed tutorials. The earlier fifth-vignette
assumption had no confirmed user choice; the four-vignette preference still binds.
Exclude `vignettes/articles` explicitly from package builds. Check the installed
pkgdown's support for this placement during G1 before adding the full page.

## Why this milestone comes next

The existing plotting APIs cover all twelve cards. The work needs no new export,
model, dependency, scientific GO, external dataset or scheduler job.
A local rendered gallery is a bounded deliverable with visible acceptance checks.
Publishing the site remains a later delivery action.

| Candidate | Current evidence | Disposition |
|---|---|---|
| Close existing CI | Run 39992327476 completed successfully; all five jobs passed. | Prerequisite cleared. |
| Plot gallery | Plan exists; article does not. Four tutorials and public plotters exist. | Next milestone. |
| GB-005 summaries and ordering | Present in R/metadata.R and R/plots.R; tests in test-metadata.R and test-plots.R. | Already implemented; do not repeat. Tests were not rerun for planning. |
| STARCAT/reference builder and GoM DE | No current public implementation; reference/scientific contracts need separate review. GoM DE needs counts and valid memberships. | Larger separate work. |
| CE G2 and release | CE_G2_LOSS_MEMO.md recommends no rerun. GO_MEMO_SCI004.md still names release prerequisites. | Neither follows from green CI; defer. |

The historical SCI-004 memo does not establish scientific recovery/calibration GO
for the broader reference plan. Its title and git subject do not prove those gates.
Keep claim-specific reference qualification separate from this documentation milestone.

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
labels and descriptive image alt text. Keep CSS scoped to the gallery; use ordinary links and no required JavaScript.
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

Keep the current 23-function public API and 57 provenance records. Do not use triple-colon calls or re-export
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
Keep file/chunk names explicit and check their paths in the standalone article preview and pkgdown site.
A clearly marked shared setup supplies the fit/view/metadata for each short recipe.

## Review findings and resolutions

| Finding | Resolution |
|---|---|
| Four tutorials are currently enforced as a four-article total. | Addressed: keep four installed tutorials; allow only the named website gallery in addition. |
| Only two plotting functions remain public in utilities. | Cover their six supported styles; omit retained internal APIs. |
| Native and utilities structure plots have different scales. | Label normalization on each card/recipe; do not present them as interchangeable. |
| Utilities distribution style is a boxplot. | Label the card as boxplots; no invented density option. |
| A fifth installed vignette conflicts with the retained tutorial limit. | Addressed: website-only article with an explicit build exclusion and tarball audit. |
| Repeated fits or duplicated snippets could drift. | One shared fit and named chunks; recipes reuse chunk source. |
| Native mean-side names differ from sample/feature semantics. | Explicit pm_which and sample_side; align IDs before grouping/covariate use. |
| Original plotting authors must remain visible. | Cite package methods and existing source credits; link CREDITS and use package citations. |
| docs/ contains tracked plans and generated site output. | Addressed: render under validation/ locally and RUNNER_TEMP/pkgdown-site in CI; never clean docs/. |
| Old plan says 21 exports and unimplemented filtering. | Addressed: 23 exports; filtering shipped at 70a30d6. Link the existing workflow without adding filter behavior. |
| Site checker ignores nested article sources and URL fragments. | Addressed: G3 adds the named article and verifies recipe anchors as well as image files. |
| Local checker scripts are ignored and absent from CI. | Addressed: move the checker and fixtures to .github/scripts; use the tracked copies locally and in CI. |
| Local installed utilities may differ from the checkout. | Addressed: G1 checks package identity; G4 installs the final built candidate before rendering and checking. |

## Progress summary

- [x] G0: Reconcile delivery and select the milestone — complete.
- [x] G1: Prepare and check one shared plotting example — complete.
- [x] G2: Build the twelve-card website article and recipes — complete.
- [x] G3: Connect navigation and update article checks — complete.
- [x] G4: Render, inspect and qualify both outputs — complete.

## Implementation steps

### [x] G0: Reconcile delivery and scope

- Status: complete. Dependencies: HANDOFF.md and current source inventory.
- Outcome/actions: verified CI run 39992327476; reconciled GB-005 against source/tests; selected the website-only gallery.
- Affected areas: this plan, HANDOFF.md and existing delivery records.
- Validation: live GitHub run reports success for full-engines, Linux release/devel, macOS release and Windows release. Source has 23 exports and four tutorial Rmds.
- Evidence: https://github.com/mdmanurung/flashier-utils/actions/runs/39992327476 ; R-devel completed 2026-10-23T17:45:32Z. Delivery prerequisite verified; implementation is now complete below.

### [x] G1: Prepare and check the shared example

- Status: complete. Dependencies: G0. Estimate: 45–75 minutes.
- Outcome/actions: create the website-only article and its build exclusion. Put one small seeded native fit/view/metadata setup in the Rmd. Inspect installed native help; verify all twelve calls before designing cards.
- Affected areas: vignettes/articles/plot-gallery.Rmd and .Rbuildignore. Read R/plots.R, man/plot_factor_*.Rd, existing vignette setups and native plotting help before editing.
- Validation: use R4_51 through scripts/run-r.sh. Verify installed utilities match the checkout. Assert rank >= 2, exact sample/metadata/covariate IDs, nonnegative activity and signed effects. Build each plot and inspect nonempty finite layer data. Compare saved fit/view coordinates before and after plotting. Prove pkgdown discovers the nested article with a small render first.
- Evidence: validation/plot-gallery-20261023/render-example.log, example.rds and api-inspection.txt. Twelve plots have nonempty finite layers; serialized fit/view bytes remain identical. Installed 2.0.0.9000 source code matches the qualified baseline. pkgdown discovers the nested article and emits articles/plot-gallery.html.
- Failure handling: if rank/support requirements fail, fix the example setup. Do not alter inference or mask errors to force a preview.

### [x] G2: Build the gallery page

- Status: complete. Dependencies: G1. Estimate: 90–120 minutes.
- Outcome/actions: add two six-card grids, twelve anchored recipes, shared setup, scale labels and credits. Store plotting calls once in named chunks; recipes reuse source without refitting.
- Affected areas: vignettes/articles/plot-gallery.Rmd; scoped CSS inside the article only.
- Validation: exactly twelve distinct real previews and twelve unique recipe targets; visible package/function names; descriptive alt text; keyboard links; readable axes. Confirm boxplots are named correctly, native PVE is qualified, structure normalizations differ explicitly and any excluded zero rows are disclosed. Use public calls only.
- Evidence: validation/plot-gallery-20261023/example-final.rds and standalone-final/plot-gallery.html; twelve finite, nonempty plots, matching recipes and exact fit/view invariance. Browser and screenshot checks passed in both formats at 375/768/1280px. A long inline expression, duplicate image links and crowded labels were fixed before final qualification.

### [x] G3: Connect navigation and enforce checks

- Status: complete. Dependencies: G2. Estimate: 45–60 minutes.
- Outcome/actions: add a gallery navbar link and article entry; retain the four tutorial entries. Add one README link. Extend the existing checker without adding a dependency. Inventory nested sources, allow only plot-gallery as the extra article, and verify local URL fragments including same-page recipe links.
- Affected areas: _pkgdown.yml, README.md, scripts/check-site.py and scripts/test-check-site.py for local work. For durable delivery, move the two currently ignored helpers to .github/scripts/check-site.py and .github/scripts/test-check-site.py; update script references. Add a step in .github/workflows/pkgdown.yaml before artifact upload. Keep receipts under validation/ locally or RUNNER_TEMP in CI.
- Validation: fixtures accept four tutorials plus this gallery; reject a fifth tutorial, an unapproved nested article, missing image, missing recipe anchor, duplicate recipe ID and stale article link. Preserve attribution and exported-reference checks. Run the checker on the actual site. CI must invoke the same tracked scripts tested locally.
- Evidence: validation/plot-gallery-20261023/site-fixtures-final.log (16 passing tests) and site-final-receipt.json (23 exports, five articles, credits and links). The helpers now live under .github/scripts, outside the local scripts exclusion; the workflow invokes them before Pages artifact upload. No active script callers retain the old helper paths.

### [x] G4: Render and qualify the final candidate

- Status: complete. Dependencies: G3. Estimate: 60–105 minutes.
- Outcome/actions: build a fresh package tarball with the four vignettes, install that candidate, and render standalone gallery HTML plus a pkgdown preview outside docs/. Use cached pancreas inputs; document any offline asset handling.
- Affected areas: validation/plot-gallery-20261023/ only, plus fixes to the scoped gallery/navigation/checker files if checks fail.
- Validation: R CMD check with vignettes and tests through scripts/run-r.sh. Audit the tarball: four installed tutorials, no website-only gallery, no study outputs/receipts/plans. Confirm the 23 exports, provenance rows and license notices remain intact. Check all twelve assets, anchors, local function links and credits. Inspect layouts at 375, 768 and 1280 pixels; retain screenshots if browser tools exist. Record visual checks as unverified if unavailable.
- Evidence: validation/plot-gallery-20261023/receipt.json and REPORT.md. R CMD check status OK; 1,577 assertions, zero failures/skips, 11 known projection warnings. Sixteen site fixtures pass. Four tutorials install; the gallery and archives are excluded. Both final pages pass browser checks at 375/768/1280px with 1/2/3 columns, twelve images/targets, visible keyboard focus and working recipe navigation. Screenshots were inspected. Final source hashes and twelve distinct image hashes are frozen in the receipt.
- Failure handling: preserve failed attempts. Rerun affected checks after fixes; qualify the final source bytes. Revert only milestone-owned edits if rollback is needed; preserve pre-existing HANDOFF.md edits and reports/core-review-20261022/.

## Milestone acceptance and delivery boundary

Completion means a checked local site has twelve working previews/recipes while
the package still installs four tutorials. Public exports, scientific coordinates
and attribution remain unchanged. No calibration or release claim follows.
New plot behavior, filter changes, GBCD/GoM DE additions and study reruns are outside scope.

Implementation and local qualification are complete. The user requested delivery;
CI/gallery commits 5a2a581 and 952ebed are pushed and remote-verified.
Package run 40079575208 passed all five jobs. Website run 40079575198 passed
build and deployment. The live gallery and all twelve PNGs passed HTTP checks.
When delivery is separately requested, inspect the final diff, stage only milestone
files, and verify hosted package/site CI for the delivered source before claiming
publication. Green CI at 70a30d6 does not qualify a future gallery revision.

## Credits and sources

- Layout reference: https://r-graph-gallery.com/density-plot.html — page organization only.
- Native API: installed flashier 1.0.60 NAMESPACE and help for flash_plot_* and plot.flash; inspection output /tmp/flashier-gallery-api.txt is temporary, so re-query installed help when implementing.
- Utilities API: R/plots.R and man/plot_factor_activity.Rd / man/plot_factor_features.Rd at the baseline commit.
- Authorship: retain FLASH/flashier and plotting-method credits through citation("flashier"), provenance(), existing manuals, CREDITS.md, inst/NOTICE and inst/provenance.csv. Do not copy upstream source or invent per-function authors.

## Scope limits and separate findings

This documentation milestone is complete locally. It adds no plotter, export or
scientific filtering rule. No public release, scheduler run, commit or push ran.
The preceding core review's approved implementation is complete. Non-lm offsets
remain explicitly unsupported; full convergence/calibration and wider upstream-version
qualification remain separate evidence work. Earlier cross-platform hosted checks passed for 70a30d6. Current gallery source 952ebed also passed all five hosted package jobs and Pages build/deployment; the live gallery is verified.

The ImmGenT iterative filter is already public as
`filter_cells_by_total_membership()` (threshold 10, 12 passes). Its checked
workflow is in advanced-use. Link that workflow if useful; add no thirteenth card
or filtering logic. Retained reconstruction is preserved by paired rescaling.

## Decisions and evidence log

- 2026-10-23 — Checked the reference page, installed native plotting exports/arguments/help, current utilities styles, navigation and article guard. Reviewed and patched this plan; implementation has not started.
- 2026-10-23 — Earlier draft requested a placement preference and assumed a fifth vignette. Superseded by the reviewed website-only decision.

- 2026-10-23 — Reviewed against main @ 09bca1c, current plotters, navigation, checker fixtures, CI workflow and delivery receipts. Replaced stale fifth-vignette/21-export/filter assumptions; converted G0–G4 into the current progress tracker. Earlier placement-request entries are historical and superseded by the website-only decision above.

- 2026-10-23 — Runtime inspection: R 4.5.1, flashier 1.0.60, ebnm 1.1.42, knitr 1.52, rmarkdown 2.32 and pkgdown 2.2.1 are available. All six native flash_plot_* exports and signatures were inspected. No plots, package tests or scientific studies ran during planning.

- 2026-10-23 — Implementation authorized through all plan tasks. The default validation/library package is stale (1.0.0.9000); install this checkout into a dedicated milestone library before testing plots.

- 2026-10-23 — G1 passed. pkgdown identifies the nested article as articles/plot-gallery in navigation configuration while emitting articles/plot-gallery.html; the config uses that verified name. Fifteen site-check fixtures pass. Final renders and browser checks remain underway.

- 2026-10-23 — Package check status OK; tests passed 1,577 assertions with zero failures/skips and 11 known projection warnings. Tarball audit found four installed tutorials, no gallery/study payload, and exact package code/notice/provenance equality. Final browser checks found one phone-width overflow from an inline metadata expression; the prose now names its attribute and field separately.

- 2026-10-24 — G4 complete after resuming. Final audit passed source/tarball inventory, all package tests/vignettes, 16 checker fixtures, twelve distinct image bytes, both browser formats and all three viewport sizes. Preserved failed attempts and original root tarball. Source snapshots and receipts stay outside package builds.

- 2026-10-24 — User requested commit and push. Created CI commit 5a2a581 and gallery commit 952ebed, pushed to main and independently verified the remote ref. All eight committed implementation files match the qualified hashes. Hosted package/Pages runs are active; qualification notes will follow after their results.

- 2026-10-24 — Delivery verified for source 952ebed: package run 40079575208 passed full-engines, Linux release/devel, macOS and Windows; website run 40079575198 passed build/deployment. The live gallery has twelve working PNGs and recipe targets. Delivery receipts are in validation/plot-gallery-20261023/delivery.json and live-gallery-receipt.json. The final notes-only commit skips duplicate CI because all eight implementation file hashes remain unchanged.
