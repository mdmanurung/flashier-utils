# Handoff — plot gallery complete; broader backlog retained

**Date:** 2026-10-24 · **Branch:** main; source 952ebed · **Status:** Gallery delivered at source 952ebed; five package jobs and Pages passed; live gallery verified.

## Goal
The plot-gallery milestone is implemented and verified. Preserve the broader
remaining-task inventory below; its scientific/reference work remains separate.

## Next action
Open https://mdmanurung.github.io/flashier-utils/articles/plot-gallery.html.
Gallery implementation and delivery are complete. Read
`validation/plot-gallery-20261023/delivery.json` for commit/run receipts.
Future package/reference/scientific work below requires its own scoped request.

## State
- Unrelated untracked reports/core-review-20261022/ remain untouched and uncommitted.
- API migration shipped at dcb5c47; GBCD and ImmGenT filter shipped at 70a30d6; Pages repaired at 23ab919.
- API: 23 exports, 34 retained internal implementations, 57 provenance records.
- Current candidate: R CMD check OK, including all four vignettes and 1,577 passing assertions; zero failures/skips and 11 known projection warnings. Receipt: validation/plot-gallery-20261023/receipt.json.
- Prior GBCD/membership focused receipts remain 63/39 passes. Current tarball matches all 142 R/man/inst/tests files plus package metadata; NOTICE, provenance and 23 exports are unchanged.
- Live CI verification 2026-10-23: run 39992327476 completed successfully; all five jobs passed, including R-devel (completed 17:45:32Z). Delivery records refreshed.
- Website run 39994153869 passed per delivery receipt; package source unchanged since 70a30d6.
- Gallery: vignettes/articles/plot-gallery.Rmd has twelve original previews/recipes. Both HTML formats pass 375/768/1280px browser checks, image decoding, anchors, keyboard focus/navigation and screenshot review. Four tutorials still install.
- Inventory is not a fresh scientific audit: old unchecked rows can describe completed work.
- In particular, git log names commit 271f583 as a SCI-004a/b GO and CE G2 loss/rerun memo; read receipts before declaring those gates unrun.

## Locked decisions
- Core API migration and signature changes were approved and shipped; do not request approval again.
- Use R4_51 through `scripts/run-r.sh`; preserve scientific contracts, original authors and provenance/notices.
- GBCD uses optional exact pin 670e3afb50236303ddf99db33a3394abd4a711bd; copy no algorithms.
- GBCD Baseline remains a factor; activity supports remain unknown after upstream sorting.
- Keep LM/ash statistics separate from EB-SNMF effects; invent no FLASH uncertainty/PVE.
- ImmGenT defaults: threshold 10, 12 passes; paired rescaling preserves retained reconstruction.
- User authorized gallery implementation, then explicitly requested commit and push. CI/gallery commits 5a2a581 and 952ebed are pushed and remote-verified. All five package jobs and Pages build/deploy passed; the public gallery and all twelve images were verified. No study, SLURM job, package release or tag was requested.
- Reviewed plan selects a website-only gallery at vignettes/articles/plot-gallery.Rmd; preserve four installed tutorials. The prior fifth-vignette assumption had no confirmed user choice.

## Dead ends — do not redo
- Pages cleanup fails in tracked docs/; build under RUNNER_TEMP/pkgdown-site.
- Baseline/GEP labels do not establish GBCD support lineage; retain unknown support.
- Descriptive coefficient margins are not conditional GoM DE lfsr.
- Software parity does not establish convergence or scientific calibration.

## Read first
- `validation/gbcd-integration-20261023/delivery.json:1` and `REPORT.md:1` — shipped scope, CI and evidence limits.
- `docs/plans/2026-10-23-plot-gallery.md:1` — G0–G4 complete with evidence; read validation/plot-gallery-20261023/REPORT.md and receipt.json first.
- `development/major-improvements-progress.md:1` — later progress overrides the original unchecked ledger.
- `reports/GBCD_STARCAT_BACKLOG_MEMO.md:1` — reconcile historical remaining tasks and scientific gate status.
- `docs/plans/2026-10-07-gbcd-starcat-audited-implementation-plan.md:183` — full task contracts/dependencies, not current completion status.
- `docs/plans/2026-10-08-consensus-ebmf-plan.md:89` — separate CE track; reconcile with later G2 evidence.

## Gallery state and delivery boundary
- G0–G4 complete. CI checker/fixtures live under .github/scripts and run before Pages artifact upload.
- Eight milestone source files changed; the plan and this handoff also changed. Preserve unrelated reports/core-review-20261022/.
- Runtime: R4_51 through scripts/run-r.sh. Current validated install is validation/plot-gallery-20261023/final-library; the default validation/library is stale (1.0.0.9000). Use the final library for 2.0 examples.
- Local builds used cached pancreas inputs and integrity-checked cached CDN assets. Browser tests used local Firefox; all logs and attempts are archived.
- Source 952ebed passed package run 40079575208 (five jobs) and website run 40079575198 (build/deploy). Live HTTP checks passed twelve cards, twelve recipe targets and twelve distinct PNGs. The notes-only follow-up skips duplicate CI; implementation hashes remain unchanged.
- GoM DE remains unimplemented; requires counts plus nonnegative row-sum-one memberships or a topic fit.

## Open items — native/reference backlog (reconcile before implementation)
- GB-005a/b already implemented: R/metadata.R and R/plots.R; regression tests in test-metadata.R and test-plots.R. Do not repeat stale unchecked ledger rows; the full package suite was rerun for gallery qualification.
- SC-001a source eligibility/identity; SC-001b generic profiles; SC-001c source-faithful STARCAT profiles.
- SC-002 all-pairs similarities without assignment pruning; SC-003a deterministic candidates/edges; SC-003b generic/STARCAT merges.
- SC-004 consensus/schema/stable IDs; SC-005 reference builder/accessors/diagnostics; SC-006 native reference projection.
- SC-007 optional source-stability import; verify exact matching source basis.

## Open items — scientific and final validation status to reconcile
- SCI-001 complete protocol stress/conditional cases; SCI-002 remaining profile/stress/residual/real fixtures.
- SCI-003a broader native/analytic qualification; SCI-003b candidate profile/graph qualification; SCI-003c reconcile scientific scope against shipped GBCD engine parity.
- SCI-004a/b: reports/GO_MEMO_SCI004.md describes defect/release prerequisites; it does not prove the broader reference recovery/calibration gates. Resolve claim-specific evidence before reference work.
- SCI-005 conditional interval/sign calibration report; SCI-006 frozen cross-technology report; SCI-008 integrated claim bridge/allowed claims.
- Extend independent native/real oracles to uncovered branches; reconcile FND-003 completion after GBCD delivery.
- VAL-001 final native/GBCD invariance after outstanding additions; VAL-002 integrated reference oracles; VAL-003 frozen integrated real replay; VAL-004 final expanded package/site/CI qualification.
- DOC-001 reconcile site-validator enforcement against recent checks; DOC-002 document future APIs/qualified claims and preserve credits; PDF manual remains locally unchecked.
- Reconcile stale ledgers with final receipts; retain failed attempts and frozen scientific outputs.

## Open items — separate conditional/research/performance tracks
- CE: reports/CE_G2_LOSS_MEMO.md records lost artifacts and recommends no rerun. Keep CE deferred; no recovery/calibration GO follows from that memo.
- SCI-007 residual truth/prediction study; RES-001 production gate; RES-002 residual API after GO; RES-003 integrated replay after GO.
- EXP-001 frozen-prior transport research remains deferred/blocked pending defined target and data.
- Rcpp: profiling pending; implement only if measurements justify it.
