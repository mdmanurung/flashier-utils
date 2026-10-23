# Handoff — GBCD and ImmGenT filtering

**Date:** 2026-10-23 · **Branch:** main @ 23ab919 · **Status:** Features pushed; verify remaining R-devel CI.

## Goal
GBCD and ImmGenT filtering are implemented and shipped to main.
Finish CI verification before choosing separate GoM DE or gallery work.

## Next action
Run `gh run view 39992327476 --repo mdmanurung/flashier-utils --json status,conclusion,jobs`.
R-devel last remained in dependency installation; record its final result in the delivery receipt.

## State
- Source: `70a30d6`; handoff: `00a5be4`; Pages fix: `23ab919`; all pushed and remote-verified.
- Working tree before this handoff: only unrelated untracked `reports/core-review-20261022/`.
- API: 23 exports, 34 retained internal implementations, 57 provenance records.
- Full installed suite: 1,575 passes, zero failures/skips, 11 known projection warnings.
- Final GBCD checks: 63 passes without warnings/skips; membership checks: 39 passes.
- Full candidate R CMD check: OK; final tarball check without repeated tests: OK.
- Source/tarball equality: 147 files; original provenance rows preserved; archives excluded.
- Pre/post and source/installed smoke identical; local Pages build and links pass.
- CI 39992327476: full-engines, release Linux/macOS/Windows passed; R-devel last still running.
- Website run 39994153869: build and deployment passed on 23ab919.
- Later commits change handoff/workflow only; package source equals qualified 70a30d6.

## Locked decisions
- User approved the core API migration and signature changes; do not ask again.
- User requested full GBCD support and the exact ImmGenT selection rule.
- Use R4_51 through `scripts/run-r.sh`; preserve scientific contracts and original authors.
- GBCD calls optional upstream pin 670e3afb50236303ddf99db33a3394abd4a711bd; copy no algorithms.
- GBCD Baseline remains a factor; all activity supports stay unknown after upstream sorting.
- Keep LM/ash signature statistics separate from EB-SNMF effects; invent no FLASH uncertainty/PVE.
- ImmGenT filter defaults: threshold 10, 12 passes; paired rescaling preserves retained reconstruction.
- User authorized commit/push, including this handoff; preserve unrelated reports and archives.
- No new scientific study, SLURM job, release or tag follows from these feature requests.

## Dead ends — do not redo
- Building Pages into tracked docs/ fails pkgdown cleanup; build in RUNNER_TEMP/pkgdown-site.
- Inferring GBCD support from Baseline/GEP labels loses prior lineage; retain unknown support.
- Treating distinctive coefficient margins as GoM DE lfsr conflates descriptive and conditional inference.
- Bounded engine parity proves software agreement, not convergence or scientific calibration.

## Read first
- `validation/gbcd-integration-20261023/delivery.json:1` — remote SHAs, CI IDs and remaining checks.
- `validation/gbcd-integration-20261023/REPORT.md:1` — validation scope and Pages repair.
- `validation/gbcd-integration-20261023/artifact-receipt.json:1` — final tarball/source hashes.
- `R/gbcd.R:1` — adapter, memory preflight, mapping and conditional signature contracts.
- `R/membership.R:1` — exact iterative selection and paired output.
- `vignettes/advanced-use.Rmd:1` — credited teaching guidance and examples.
- `docs/plans/2026-10-23-plot-gallery.md:1` — committed gallery plan.
- `development/major-improvements-progress.md:1` — broader backlog; GB-004 software qualified.

## Open items
- R-devel CI — verify completion; diagnose a failure rather than assuming success.
- GoM DE — not implemented; requires counts and nonnegative row-sum-one memberships or a topic fit.
- Plot gallery — plan saved; page implementation remains separate work.
- Scientific calibration — not run; software tests and CI do not establish it.
