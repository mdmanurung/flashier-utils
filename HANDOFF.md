# Handoff — flashier.utils 2.0 API migration

**Date:** 2026-10-23 · **Branch:** main · **Status:** Locally qualified; user authorized commit and push on 2026-10-23.

## Goal
The correctness review, four-tutorial consolidation and approved API migration are complete.
Preserve scientific contracts, attribution, existing work and validation archives.

## Next action
Read `reports/core-review-20261022/MIGRATION.md:1` and the final integrity receipt below.
Read validation/core-api-migration-20261023/delivery/receipt.json for committed and remote-verified state.
If that receipt is absent, finish the already authorized commit/push after checking Git state; do not re-ask permission.
Do not ask again for approval of the export reduction or nine signature reorders.

## State
- Development version 2.0.0.9000: 21 public functions; 34 implementations/manuals retained internally; all nine approved argument orders applied.
- Commit sequence: correctness fixes, tutorial consolidation, approved API migration; review reports and validation archives remain local.
- Source and final installed suites: 1,475 assertions each; 0 failures/errors/skips; 11 known native projection iteration warnings.
- Final tarball: `validation/core-api-migration-20261023/flashier.utils_2.0.0.9000.tar.gz`; R CMD check --no-manual: Status OK, examples and four rebuilt vignettes including cached pancreas.
- Fresh native matrix 514/514; 120 function bodies/defaults unchanged by migration; 111 existing direct caller bindings preserved; forwarded dots exercised by tests.
- Exact pre/post migration and source/installed smoke parity for activities, effects, offsets, factor IDs and coordinate manifests.
- New site passes 21 public references/four articles, attribution/local links; six checker fixtures pass. All 55 source records and original attribution bytes retained.
- No new release, tag rewriting, scheduler study, external reviewer call or alarm; inspect hosted checks separately after delivery.

## Locked decisions
- User approved both the 55-to-21 export reduction and nine proposed signature reorders: "approve both" (2026-10-23); approval is saved in the archive.
- Use R4_51 through scripts/run-r.sh; scientific calculations/defaults remain unchanged by the migration.
- Preserve original authors, CREDITS.md, inst/NOTICE, upstream licenses and inst/provenance.csv; internalization does not remove authorship.
- Scientific outputs/receipts remain outside package builds; preserve prior candidates and failed attempts; never rewrite published tags.
- User explicitly requested commit and push on 2026-10-23; no new studies, release, tag rewriting or destructive action follows.

## Dead ends — do not redo
- Earlier root/1.1 tarballs predate this migration; use the final 2.0 tarball/receipt.
- The first 2.0 candidate passed checks but had a README NEWS.md link rendered as missing NEWS.html; final link/site/tarball are corrected and rechecked.
- The provenance plan is stale: regenerating only from it loses five records and newer pins. The fixed local generator preserves all 55 registry records byte-for-byte.
- CSV test reports omit the nested result column; full RDS is preserved. Build-generated DESCRIPTION fields need DCF comparison, not raw bytes.
- Cached site assets and documented offline metadata overrides are required locally; iteration-limit parity and the small prior study do not establish convergence/superiority.

## Read first
- `reports/core-review-20261022/MIGRATION.md:1` — applied API change and external-caller migration.
- `reports/core-review-20261022/proposed-arguments.csv:1` — pre-migration and approved signature orders.
- `reports/core-review-20261022/REPORT.md:1` — correctness findings and current qualification.
- `validation/core-api-migration-20261023/integrity-receipt.json:1` — final source/tarball/artifact hashes and limits.
- `validation/core-api-migration-20261023/approval.json:1` — explicit approval and proposal hashes.

## Open items
- Verify delivery/receipt.json against live remote; local scripts/review reports/validation archives are outside the committed package scope.
- PDF manual and hosted/cross-platform CI remain unchecked; old positional calls and removed exports require external-caller migration.
- Broader GBCD/STARCAT/consensus scientific gates are separate; consult flashier-gbcd/TRACKER.md and flashier-starcat/08_AGENT_HANDOFF.md only if resumed.
- Historical scheduler counters are not current evidence; do not restart old controllers or frozen studies from previous handoffs.

Prior handoff and pre-migration package/candidates are preserved in the validation archives.
