# Handoff: flashier.utils improvements 1–3

## Reality check 2026-10-16 — session-local evidence trees are gone

The following HANDOFF references are NOT in the repository, and their
session-local copies are lost: `validation/` (all gbcd-starcat implementation
receipts; the consensus-ebmf CE study roots including `ce04-g2-ceiling-v2`
with its 180-fold results, receipts and SHA256SUMS manifests),
`development/major-improvements-progress.md`, `docs/plans/*` (including the
2026-10-07 audited 46-task gbcd-starcat plan and the consensus-ebmf
plan/proposal docs), `scripts/run-r.sh`, and the `benchmarks/consensus-ebmf-*.R`
plus `benchmarks/workflow/consensus-ebmf-*.smk` helpers. Consequences:

- The CE G2 ceiling study has no recoverable results; any verdict requires a
  full reconstruction and rerun (~550 CPU-h originally reserved). No
  scientific verdict from CE G2 exists or is claimable.
- `flashier-gbcd/` and `flashier-starcat/` (tracked) are the authoritative
  specs/trackers; SCI-004a/b evidence must be regenerated from them.
- The five delivered commits (91c9fc8..0e0bdc0) are on this main; their
  FND-001/FND-002 closure is being verified by fresh oracles rather than by
  the lost receipts.

## Committed implementation and delivery (2026-10-14)

The user authorized commit and push. The completed software stages after
1dbf193 are now recorded as five focused commits:

- 91c9fc8: projection coordinates, metadata identity and foundation boundaries.
- 54c80a3: feature intercepts, affine projection and semi-nonnegative fitting.
- f09fc74: initialized and one-side anchored native fitting.
- 76a1fb3: versioned factor support, zero-rank guards and intercept replay.
- 0e0bdc0: seeded native multistart comparison.

Each staged tree was exported without changing the working files and tested
through scripts/run-r.sh. Focused checks passed 180, 216, 109, 305 and 105
assertions respectively, with zero failures/errors/skips. Only the foundation
checks emitted the 11 existing native projection iteration warnings. All 151
package source, manual, test, vignette and attribution files match the sealed
multistart snapshot that passed the full 1,170-assertion suite. Scientific
archives and validation receipts remain excluded from Git and package builds.

Delivery evidence is local at
validation/gbcd-starcat-implementation/commit-push-20261014-v1/; its final receipt
records HEAD, remote verification and the clean-tree check after delivery.
Check that receipt and the live remote before resuming publication work.
The limited package-check and rendering limitations below still apply. This
delivery does not establish scientific GO or release qualification.

The user also asked about Rcpp feasibility. Source inspection identified
.draw_summary() and correlate_factor_markers() as profiling candidates. No
performance profile, speedup measurement or C++ implementation has been made.
Profile representative data before replacing kernels; retain scientific results,
quantile definitions, missing-value behavior and upstream attribution.

## Verified multistart stage (2026-10-14)

User resumed during foundation closeout; continued directly to GB-006. The full
three-improvement goal remains open. At stage closeout, main had inherited and
new uncommitted changes after 1dbf193; no commit, push or release had been performed.

`fit_multistart_ebmf()` is implemented around public native FLASH initializers and
fits. Unique initializer seeds and explicit native sign constraints are recorded.
softImpute/IRLBA require signed priors; IRLBA requires dense complete input. Candidate
failures, native ELBO, rank and numerical/initializer digests are retained. Maximum
finite ELBO wins, exact ties choose the smallest start ID, and select=none leaves
no selected best fit. keep_fits=FALSE keeps the selected fit and diagnostics.
Initializer traces include native validation probes. All-fail conditions expose
candidate reasons. All candidate fits require explicit stability replay callbacks;
legacy serialized multistart metadata also refuses automatic ordinary replay.
Initializer arguments are forced eagerly so saved callbacks cannot retain the
caller data/candidate list. Large-versus-small payload checks and native-fit RDS
round trips pass. The earlier functional pass and failed size-threshold attempt
are preserved; developer source metadata is not mistaken for retained input.

Final source-loaded R4_51 suite: **1,170 passed**, zero failures/errors/skips,
11 unchanged native projection iteration warnings. Native parity covers both
orientations, constrained starts and installed spectral helpers. Initializer seeds
and constraints are observed at actual native calls; failures, RNG restoration,
finite selection/ties and numerical reproducibility have independent fixtures.
Current engines: pinned flashier 1.0.60 and ebnm 1.1-42; optional softImpute 1.4-3
and irlba 2.3.7. Original authors are credited in source, manuals, provenance,
CREDITS and the credits/first-fit articles. NOTICE is unchanged.

Evidence: `validation/gbcd-starcat-implementation/multistart-20261014-v1/`.
Legacy real replay `legacy-v10` matches `legacy-v2` at 1e-8 for all published
coordinates, rankings, contrasts and query projection. Final staged package build
and content/metadata audit passed. Final R CMD check completed with zero errors
and two warnings for unbuilt vignette outputs/absent inst/doc. Installation,
code/docs and all examples passed; tests/manual/vignettes were skipped there,
and the complete source suite is separately recorded. delivery-receipt.json and
SHA256SUMS bind final source/evidence; verify manifests before continuation.
The support/refit archive below is sealed and checksum-verified; do not modify it.
No site render, scientific GO, release or CE study action follows.

Next: GB-005 activity prevalence and bounded within-group MDS ordering, then excluded
fixed-B/profile candidates before SCI-004a/b-gated core reference integration.
Optional GBCD still needs exact-pin installed-engine qualification and nonnative
mapping contracts. Follow development/major-improvements-progress.md; preserve the
full original 1–3 scope and every frozen study/failed attempt.

---

## Verified support/refit stage (2026-10-14)

Continue the authorized three major improvements; the full goal remains open.
Existing uncommitted changes on main after 1dbf193 were preserved, and the next
foundation stage has been implemented. No commit, push or release was performed.

Completed in this continuation: versioned per-factor native support/role records
for intercept, initialized and projected fits, with exact ID/native-position
validation, conservative global summaries and hashed serialization. Ordinary legacy
manifests and basis hashes remain unchanged. Native edits preserve records by ID;
projection retains fixed-side support and records the estimated-side prior.
Support/role overrides in a supplied view are refused for a versioned native fit.
Nonnative-engine mapping remains with the optional GBCD adapter; do not invent
FLASH indices for such outputs.

FND-005/GB-002 zero-rank and replay boundaries are complete: supported extraction
and summaries stay typed/empty; selection-based feature ranking, association,
backprojection, network/display/reduction and stability fail early without biological
factors. Default intercept-only native backfitting remains available. Automatic
intercept bootstrap and pure reorder match direct native FLASH recipes in both
orientations. An incomplete stored intercept recipe requires an explicit callback.

Final source-loaded R4_51 suite: **1,065 passed**, zero failures/errors/skips,
11 unchanged native projection iteration-limit warnings. Current source, logs,
RDS/TSV results and preserved attempts are under
`validation/gbcd-starcat-implementation/support-refit-20261014-v1/`.
Real replay `legacy-v8` matches frozen `legacy-v2` at 1e-8 for reference/signed
coordinates, rankings, contrasts and query projection, with all 14,528 genes and
48/24 cells. These real data have no independent donors or affine examples.
Final staged source build and exact R/man/inst/tests content audit passed, retaining
NOTICE and excluding scientific/planning payload. Final R CMD check completed
with zero errors and two warnings for unbuilt vignette outputs/absent inst/doc.
Installation, code/docs and all examples passed; tests/manual/vignettes were skipped
there. The complete source suite is recorded separately. No site or release
qualification. delivery-receipt.json and SHA256SUMS bind final sources and evidence.

Next: follow `development/major-improvements-progress.md` for GB-006 multistart,
then GB-005 and excluded fixed-B/profile candidates. SCI-004a/b GO is still required
before core reference integration; optional GBCD needs exact-pin installed-engine
qualification. Do not narrow the full 1–3 goal to completed wrappers/foundations.

Live scheduler inspection found no queued/running CE study jobs. Its G2 completion
marker is absent; no final scientific verdict/accounting audit was run. Frozen CE
roots and source were not modified, restarted or resubmitted. Older G2 counts below
are historical. The evaluator binds frozen hashes despite reading the live checkout.

---

## Previous verified package progress (2026-10-13)

User authorization: "work on the above major improvements from 1-3 please",
followed by "resume". Scope remains correctness repairs, expanded native fitting
(semi-nonnegative, intercept, initialized/anchored, multistart and optional GBCD),
and reusable cross-dataset program references/projection. Scientific gates in the
audited plan still apply. This authorization supersedes older statements below
that package implementation was not authorized. Keep the full goal open.

Current source changes remain uncommitted on main after 1dbf193. The latest
source-loaded R4_51 full suite passed **942 checks**, zero failures/errors/skips,
and the same 11 pre-existing native projection iteration-limit warnings. Evidence:
`validation/gbcd-starcat-implementation/native-fitting-20261013-v1/tests.log`,
`tests.rds`, `tests.tsv` and receipt.txt. Documentation and the initialized-fit
manual example ran successfully. Source and failed attempts are archived with SHA256SUMS. A package-only staged
source build passed exact R/man/inst/tests content parity and retained NOTICE, with
no scientific/archive/planning payload. Root R CMD build was interrupted during
its whole-checkout archive copy before exclusions; its log remains. R CMD check
completed with no errors and two warnings for unbuilt vignette outputs/inst/doc.
Code/docs/install/all examples passed; tests/manual/vignettes were skipped in that
limited check (the full source suite is recorded separately). No site render or
release qualification follows.

Implemented: projection-native maps/scaling/PVE/PSD/draws; measured-axis lookup
and explicit unavailable uncertainty; external-offset serialization; map validation;
metadata/bootstrap tuple identity; stability recipe guards; typed empty summaries;
audited prior families; opt-in feature intercepts; affine new-feature projection;
`fit_seminonnegative_ebmf()`; `fit_initialized_ebmf()` with one-side anchoring,
mandatory backfits, stable IDs/roles, compact input hashes and selected-only safe
nullchecking. Direct public-native comparisons cover both orientations, new-factor
addition/removal, moments/draws, raw/LDF reconstruction with D on either side,
intercept-only cases and RNG restoration on native success/failure. Intercept and
initialized selectors preserve the native map; the protected intercept is hidden
from semantic factor IDs. Manuals, exports, provenance and author credits match.

The full versioned per-factor support schema and remaining zero-rank refusal/
automatic intercept-bootstrap boundaries remain open. Native parity is software
qualification, not calibration, predictive benefit or consensus recovery GO.

Published real outputs were frozen under legacy-v2; latest legacy-v5 replay passed
at 1e-8 for reference/signed coordinates, rankings, group contrasts and two-factor
query projection, all 14,528 genes and 48/24 cells. The data have no independent
donors or affine examples. Earlier legacy-v1 transient basis IDs were corrected
in v2 and all attempted receipts remain.

SP-v1 claims/gates/seed reservations are frozen in protocol-v1. fixtures-v1 adds
independent hand/source/query truths and masks. native-v2 verifies native hand
reconstruction, intercept ones, protected lifecycle and semi-nonnegative capability.
SCI-001/002/003 remain partial; stress/profile/real/residual cases and all core
claim-specific SCI-004a/b gates are pending. No scientific GO, release or push.

Next: follow `development/major-improvements-progress.md`; finish FND-003/005
schema/zero-rank/refit boundaries, then GB-006/005 and excluded fixed-B/profile
candidate gates before SC integration. Optional GBCD requires installed exact-pin
qualification. Do not narrow the full 1–3 goal to completed fitting wrappers.

Read-only scheduler inspection of 181 distinct submitted G2 study IDs found no
queued/running jobs. Its evaluator uses the live checkout with frozen source hash
checks. No terminal scientific verdict/accounting audit was run. Older G2 counts
below are historical; do not restart/resubmit or rewrite frozen study roots.

---

## Historical G1/G2 handoff (read-only study context)

Last checked: 2026-10-10 12:12 UTC (14:12 Europe/Amsterdam). G1 artifacts are finalized; read them without overwriting.
Workspace: `/exports/para-lipg-hpc/mdmanurung/flashier-utils`.

## 1. Goal and immediate next action

**G1 is complete; G2 approved and all180 primary folds submitted:** all30 pilot jobs completed/verified, summary and accounting finalized, final evaluation SHA256 manifest passed, and independent summary/accounting audit passed. Read `docs/plans/2026-10-09-consensus-ebmf-g2-ceiling-proposal.md`. It recommends60 fresh primary-regime datasets/180folds within a separate550CPU-hour ceiling, using a single exact95% baseline-shortfall bound to check whether the required5pp recovery gain is unattainable. It permits NO-GO for this generator's recovery gate or INCONCLUSIVE, never GO. All60 primary inputs and60300initializer/bootstrap seeds are frozen. First fold25742299 completed/verified and the controller released allremaining179folds after held-allocation checks. Current snapshot:36completed/verified,2running,142pending;11complete independent datasets. No primary job/method failures among verified folds. The eventual scientific endpoint is GO, NO-GO or INCONCLUSIVE; public consensus integration requires GO plus real-donor evidence and foundation repairs. G1 completion does not establish benefit or clear package integration.

**First action for a zero-context agent:** read the G2 ceiling proposal and its planning self-check at `validation/consensus-ebmf/g2-design-20261009/arithmetic.log`; G2 is approved; check the live G2 controller (exec session44555) and log `validation/consensus-ebmf/g2-run-20261010T0715Z.log`, then follow section8. Do not duplicate/resubmit jobs or ask for G2 approval again. A durable proposal copy is in that planning archive. Do not run the status helper, finalizer, summarizer or per-job verifier against the finalized evaluation archive: they overwrite frozen files. Existing public-package defects remain a separate workstream; G1 did not authorize those package edits.

Latest user requests: approve, then resume (2026-10-10). Approval covers the exact60-dataset/550CPU-hour ceiling proposal. Separate G2 source/input/protocol/manifest frozen under current `validation/consensus-ebmf/ce04-g2-ceiling-v2/`; failed v1 preflight is preserved/frozen; approval.json binds the proposal hash. Guard7mockedtests, allocation refusal checks, exact decision/failure/denominator checks, helper parity and Snakemake dry-run passed. No new user permission is needed within G2. Infrastructure25742287 failed an overly strict whole-native-object identity check after writing its checkpoint; numerical coordinates reproduced exactly. Its1830CPU-second reservation and13allocatedCPU-seconds remain charged. Corrected v2 infrastructure25742297 passed cross-job numerical checkpoint reuse with0native compute calls and changed-configuration refusal (12allocatedCPU-seconds). All30 G1 tasks completed; no resubmission/retry was needed. Do not ask the user to approve G1 again. G2 ceiling-check execution is approved; full benefit/harm and G3 real-data designs remain separate. Earlier parallel subagent work was authorized; no subagent was used in this planning turn.

## 2. Workspace, preferences and authorization

- Read root `AGENTS.md`; activate installed `i-have-adhd:i-have-adhd` at a new session. Keep concise next actions and bounded work. Ponytail full is active. Never suggest changing models unless the user asks.
- Use **R4_51 through `scripts/run-r.sh`**. The wrapper already invokes Rscript: `scripts/run-r.sh script.R` or `scripts/run-r.sh -e '...'`; do not put `Rscript` after the wrapper. Load checkout with `pkgload::load_all(".")`. The validation-library installed package was stale; plain `library(flashier.utils)` is insufficient.
- Use `/usr/bin/python3.11`; default `python3` is3.6. Snakemake executable: `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/snakemake/bin/snakemake`.
- Branch `main`, HEAD `75a1709e9979f66fa46beaee994443b804f9ca64`. Source version `1.1.0.9000`; local release tags `v1.0.0`, `v1.1.0`. Git status: modified `HANDOFF.md`, pre-existing untracked `flashier-gbcd/` and `flashier-starcat/`. No public source/API edits in this session. Do not rewrite tags, discard these folders or commit local scientific artifacts.
- Scientific scripts/plans/receipts live in ignored `benchmarks/`, `docs/plans/`, `validation/`, excluded from builds. Preserve credits, `inst/NOTICE`, licenses and `inst/provenance.csv`; keep CREDITS/vignette credits consistent when attribution changes.
- Writable roots: repository and `/tmp`; `.git` read-only. SLURM socket access requires approved sandbox escalation. Use the exact approved G2 scope (60datasets,550CPU-hour ceiling) in justification, not another user permission request. Auto-review has approved these read/submit/release actions; no rejection remains.
- Parallel work explicitly authorized. One subagent completed read-only evaluator/precision/package reviews; older agents hit usage limits. Their status does not substitute for receipts. Give explicit ownership and do not edit frozen source concurrently.
- G1 approval: user “i approve”,2026-10-09, **100 allocated CPU-hours**, including infrastructure, seed155 smoke,30 pilot tasks, failure/cancellation/retry grace. G2 primary precision/budget, G3 real-data cohort/preprocessing, G4 optional early package repairs, G5 verdict adoption and G6 publication remain separate decisions. Prepare concrete designs before requesting those decisions.

## 3. Package status: implemented versus verified versus pending

| Workstream | Current evidence | Still open |
|---|---|---|
| Existing package |50 exports and 50 provenance rows;19 tutorial Rmds; release tags through1.1.0; development1.1.0.9000. |Existing correctness defects and build/site assurance below. |
| CE-02 native multistart |`ce02-development-v3` frozen, completed. |Public multistart API is not implemented; excluded helper only. |
| CE-03 candidates/refits |Development-v3 passed15 checks; independent verification-v1 passed8; reused30 saved v1 bootstraps. |General edge cases/public integration remain distinct; no scientific benefit claim. |
| CE-04 development |All30 tasks verified;150primary/210sensitivity rows; final archive SHA and independent summary/accounting audit PASS. |Scientific precision and primary/real-data gates remain open. |
| CE-04 G2 ceiling check |Approved60datasets/180folds/550CPU-hour ceiling; v2 checkpoint preflight passed; first primary fold running. |First full verification, remaining179folds, final ceiling report; no GO permitted. |
| CE-04 full benefit and real donors |G1 descriptive evidence and design reviews exist. |Full GO precision/reliability and real-data protocol remain unresolved. |
| CE-05 consensus API |Not implemented. |Requires scientific GO, real-donor evidence and FND prerequisites. |
| GBCD/STARCAT extensions |Audited implementation plan and original bundles exist. |Package features not started. Within-dataset CE consensus is different from cross-dataset SC catalogs. |

Public-package source was re-reviewed today and remains unchanged. The October8 audit's native defect probes and full-suite receipt are in `validation/feature-completeness-review-2026-10-08/` and `docs/plans/2026-10-08-feature-completeness-review.md`. Its old CE-03/CE-04 progress paragraphs are superseded by this handoff/current receipts; its package defects remain open:

1. **FND-001:** projection inherits source native indices/scaling/PVE; posterior extraction uses positional maps; retained views and external-offset/intersection uncertainty handling need repair. See `R/projection.R`, `R/extract.R`, `R/representation.R`. Independent fresh-native mean/SD/PVE/draw oracles are required.
2. **FND-002:** `R/operations.R:129` passes `.select_ids()` integer positions as native `fixed_idx`, rather than full-axis logical masks. Whole-side/mean tests do not cover this defect; verify partial fixing in both orientations.
3. Additional audited defects: metadata tuple collisions, valid zero-factor consumer crashes, ordinary refit recipes inherited by edited/projected fits, insufficient map validation. FND-003/004 and relevant FND-005/006 contracts remain open.
4. Local generated site is stale (1.0.0/13 articles versus current19 tutorial sources). `scripts/check-site.py` hardcodes34 exports/13 articles and is not called by deployment workflow. `.Rbuildignore` lacks exclusions for HANDOFF and the two local planning bundles. This is a future-build risk, not proof published archives contain them.

Package-work priority is SCI-001 freeze → FND-001/FND-002 → remaining contracts/docs/build fixes, under the relevant plan/scope. G1 does not authorize those substantial public-package edits. Audited extension authority: `docs/plans/2026-10-07-gbcd-starcat-audited-implementation-plan.md` (46 tasks: audit complete, implementation not started). Do not present the plan as implemented code.

Historical full checkout suite: October8 exit0, no failures/reported skips,11 existing native iteration-limit warnings. That build attempt was interrupted and produced no archive. **No fresh full-suite, R CMD check, site render, live CI/release verification or new defect numerical reproduction was run in this status/handoff turn.** Public-site state may differ from stale local `docs/`.

## 4. Completed G1 execution and budget

**30 completed/verified;0running/0pending.** Postprocessing finished17:18:30UTC. Worker log ends with archival completion. No pilot retry/cancellation occurred; the earlier intentional infrastructure cancellation remains in accounting. All10 independent datasets have all3folds.

| Regime/task | Jobs | State at snapshot |
|---|---|---|
|separated|Exact IDs in submissions.jsonl|All6 completed/verified. |
|identifiable_donor_heterogeneity|Exact IDs in submissions.jsonl|All6 completed/verified; first task reused155 seed-smoke fits. |
|near_collinear|Exact IDs in submissions.jsonl|All6 completed/verified. |
|rare_donor_associated|Exact IDs in submissions.jsonl|All6 completed/verified. |
|null|Exact IDs in submissions.jsonl|All6 completed/verified; inexpensive_best_5 failed conditioning in null-t2-f3, preserved as a method failure, not a failed job. |

Exact mapping for **all 30** tasks/IDs, accepted reservations and releases: `validation/consensus-ebmf/ce04-development-v1/submissions.jsonl`; frozen manifest `tasks.json` alongside it. Never infer the full mapping from an ID range. Final `ce04-evaluation-v2/live-status.json` is hashed and must not be refreshed in place.

- Reserved 339390CPU-seconds = **94.275 CPU-hours**; remaining 20610seconds = **5.725 hours**. Reservations retained for every attempt; no automatic refunds/requeues/retries.
- Final actual allocation:41913CPU-seconds = **11.6425CPU-hours** across all33 accepted infrastructure/smoke/pilot jobs. `budget-status.json`, completion marker and final `accounting.tsv` agree. Independent exact-ID accounting audit passed. Older provisional snapshots are retained.
- Interim report `validation/consensus-ebmf/g1-interim-20261009T152552Z/README.md` preserves the earlier17-task snapshot, with85primary/119sensitivity rows, five complete datasets and zero method failures in allfive arms among verified tasks. Receipt/protocol hash assertions passed. Fresh-fold timing excludes reused smoke; rare/null timing/outcomes are unavailable in that snapshot. This report is outside the final evaluation archive; do not rewrite it when later jobs complete.
- Final report `validation/consensus-ebmf/g1-final-report-20261009/README.md`:both consensus arms0method failures; nested6valid no-program null folds; inexpensive baseline1conditioning failure retained with infinite loss/NAFDP. Both consensus arms gained0pp primary recovery against best155 (100% recovery in both pilot datasets); seed/nested primary RMSE changed+0.29%/+0.63%. These are descriptive two-dataset comparisons, not a verdict.
- Each pilot job:1CPU,4G request,medium,3h max,no-requeue; reserve10830CPU-seconds including30s grace. Both planned3h retry slots cannot fit remaining ceiling: at most one additional such reservation fits. Retry wrappers are not prepared; guard refuses duplicates and over-budget reservations.
- Never cancel unrelated jobs or lower frozen resource requests merely to bypass quota. If a task fails, inspect its exact job/log/checkpoints and preserve the attempt. A retry needs a distinct manifest slot and unchanged method/checkpoint configuration. Current worker/summary assumes original30 job receipts; a retry requires a separate completion report retaining original failures, not silent substitution or rewriting frozen receipts.

## 5. Resume commands and verifier lifecycle

**G1 finalization is complete. Do not execute the restart/manual verification commands below against this archive.** They are historical recovery guidance. Final evaluation manifest passed `sha256sum --check --quiet SHA256SUMS` from its archive root. Independent R4_51 `g1-final-report-20261009/audit.R` passed recalculation of30receipts/50dataset rows/30paired rows/failure denominators/33accounting IDs. Its first relative-tolerance CSV roundoff failure is retained; final paired comparison uses absolute1e-12 tolerance. No frozen output changed.

User requested a check in two hours on2026-10-09 at15:31UTC. The detached checker (formerPID799550) stopped without producing its scheduled17:31:15UTC snapshot. Process absence and empty log confirmed; ran the unchanged read-only checker directly at17:51:59UTC. `validation/consensus-ebmf/status-check-20261009T173115Z.json` records the scheduled and actual times, all30ended/receipts present and completed postprocessing; the recovery log has the same stem plus`-recovery.log`. It changed no jobs or frozen outputs. User was told the scheduled background check stopped; do not claim it executed on time. No automatic chat notification is available. Receipt counts indicate presence, not a new verification pass.

Run from repository root. First check completion without modifying archived outputs:

```bash
cd /exports/para-lipg-hpc/mdmanurung/flashier-utils
ls validation/consensus-ebmf/ce04-evaluation-v2/POSTPROCESS_COMPLETE.json
tail -n 12 validation/consensus-ebmf/ce04-postprocess-resume-20261009T1309.log
```

A missing completion marker is normal while jobs run. If it **exists**, read marker and `summary/receipt.rds`, `dataset-metrics.csv`, `paired-development-differences.csv`, `failure-denominators.csv`, `compute-ledger.csv`, `accounting.tsv`. Verify the final checksum manifest from its archive directory. Do not rerun helpers that overwrite `live-status.json`, summary or receipts after archival freeze.

If not complete, refresh exact-ID status (SLURM escalation required):

```bash
/usr/bin/python3.11 validation/consensus-ebmf/ce04-evaluation-v2/launcher/status-ce04-pilot.py
```

**Automatic worker completed:** `benchmarks/finish-consensus-ebmf-pilot.py`, former exec session75637, final log above. Restarted at13:09UTC after confirming the original process was absent; verified the remaining receipts and finalized all outputs at17:18UTC. Original `ce04-evaluation-v2/postprocess.log` is preserved. No restart is needed or safe against the finalized archive.

If the worker really stopped and archive is incomplete, inspect its error first. Restart only this no-submit/no-retry command, preserving evidence/logs:

```bash
/usr/bin/python3.11 benchmarks/finish-consensus-ebmf-pilot.py
```

Keep it in a resumable exec session with output captured to a **new** log outside `ce04-evaluation-v2/`; never overwrite prior error evidence. The finalizer excludes only the original `postprocess.log`, so a new open log inside that archive would invalidate its final checksums. The `--once` check is only safe while incomplete; after all tasks verified it enters finalization and rewrites outputs. Avoid invoking it after completion.

Manual per-job verification, only while outputs remain unfrozen and no worker is verifying the same job:

```bash
scripts/run-r.sh benchmarks/verify-consensus-ebmf-ce04.R   validation/consensus-ebmf/ce04-evaluation-v2/<task>/job-<job_id>
```

Worker failures are deliberate for ended jobs without receipts, mismatched hashes or invalid accounting. It submits/retries/cancels nothing. Summarizer `benchmarks/summarize-consensus-ebmf-ce04.R` requires30 independently verified original task receipts; it cannot fabricate completion for a missing/failed job. Preserve failures and issue a qualified partial report if completion cannot be established.

Final archive manifest excludes runtime `cache/`, `.snakemake/` and still-open `postprocess.log`; scientific receipts/source/ledgers are included. Verify checksums before changing any final artifact. This handoff and accounting snapshot are outside that frozen method archive.

## 6. Frozen algorithm and artifact map: do not change during jobs

Method authority: `validation/consensus-ebmf/ce04-evaluation-v2/protocol.rds` plus `protocol.sha256`, source snapshots and source/artifact hashes. IDs and generator streams are fixed. Inputs bind all pilot fixtures, generator receipts/protocol, verified smoke receipt and155 checkpoints. Source changes can make pending tasks fail preflight; no public R or excluded method changes while jobs depend on those bytes.

- Pilot:2 generated trials ×5 regimes ×3 donor folds =30 tasks. Generator has12 donors/72 samples/48 features, variable donor cell counts and repeated visits. Bootstrap complete donor occurrences within strata; no donor crosses train/test folds. Primary and pilot generator seed streams are disjoint; no primary draws generated.
- Five arms:ordinary_native(seed666),best_155,seed_consensus_155,nested_30_by_5,inexpensive_best_5. Shared306 discovery attempts/fold =155 full-training starts +30×5 bootstrap starts +1 ordinary fit. Logical charges1/155/155/155/5. First primary-regime trial/fold reuses the seed-smoke pool and needs151 new attempts, not another155 full-training fits.
- **Both consensus ranks/anchors come from best eligible first requested 5 starts**, not best155. Registered seed q=min(N−1,46) when N>1 and support>=124 distinct starts; one vote/start/cluster. Basin votes diagnostic only. Do not alter voting to rescue support or tune using held-out results.
- Frozen native model:rows,point-exponential priors both sides,max_factors4,var_type0,S=NULL,nullcheckFALSE,greedy maxiter500/backfit2000,no warmstart/extrapolation; explicit seed passed to native greedy initializer. Eligible highest finite ELBO, smallest start-ID tie; valid zero-rank winner remains valid. Same-solution rule:equal rank,abs ΔELBO<=1e-4,max fitted Δ<=1e-3. ELBOs are never compared across different bootstrap datasets.
- Identified anchored regimes use strict congruence>0.90 one-to-one truth recovery/FDP. Near-collinear is **nonidentifiable**:subspace diagnostics only, no family recall/FDP claim. Null has no structured truth families but positive censored-noise background; all estimated nonzero calls are false. Censored-Gaussian observations deliberately mismatch native Gaussian residual likelihood.
- Masked scoring:odd inference/even scoring features;macro-donor RMSE. Every supported arm checks scoring-only+37 perturbation leaves actual masked inference/activity/tau unchanged. Conditioning:inference rank relativeSVD1e-8,rcond>=1e-6. Zero-program predictor has no estimated held-out native noise. Failed candidates never become zero/baseline; failed recall0 for identified nonnull families, FDP unavailable, RMSEInf.
- Seven registered one-at-a-time sensitivities:congruence0.80/0.95,recovery0.60/0.90,density0.25/0.75/Inf around primary0.90/0.80/0.50. Reuse discovery fits; never choose a winner from these scores. Dataset is Monte Carlo unit:average3folds within each trial before paired SD across2trials.
- Evaluator uses source-only local clones to repair empty-reference alignment and empty-match aggregation; public package code unchanged. `ce03_profiles()` sorts stable member IDs, but generic direct seed-consensus input-order invariance is **not established**; October8 raw-helper counterexample remains a public-integration hardening task. Passed fixture checks are not universal proofs.
- Snakemake owns dependency/staleness detection via `benchmarks/workflow/consensus-ebmf-ce04.smk`, one per-job directory/lock and writable XDG cache. Frozen `.sbatch` and submit guard must remain unchanged.

| Archive | Meaning |
|---|---|
|`ce01-development-v2`, `ce02-development-v3`|Baseline protocol/software/native multistart evidence; retained immutable. |
|`ce03-development-v1`, `ce03-development-v2`|Failed runner/finalizer attempts; v1 contains all 30 reusable donor checkpoints. Do not erase or rerun blindly. |
|`ce03-development-v3`, `ce03-verification-v1`|Passed15+8 fixture checks; both-orientation conditional refits; software evidence only. |
|`ce04-development-v1`|Approved fixed-task manifest, append-only reservations/submissions/releases and budget guard tests. |
|`ce04-infra-v1`|Create25735341 deliberately cancelled after saved checkpoint; resume25735347 completed with zero checkpoint-compute calls. Matching hashes;125 total allocatedCPU-seconds. |
|`ce04-generators-v1`|Unused nonidentifiable generator proposal preserved; no method fits used it. |
|`ce04-generators-v2`|Current frozen10 pilot fixtures/oracles; pure A/B anchors establish minimal noiseless rank-K identifiability in identified regimes, not noisy/higher-rank recovery guarantees. |
|`ce04-seed-smoke-v1`|Job25736334 completed344s;155 eligible fits; first-five anchor rank4;517 profiles survive;3 clusters155votes retained, fourth52 excluded. Independently verified155 hashes/reuse/refusals/reconstruction, SHA frozen. |
|`ce04-evaluation-v1`|Unused failed Snakemake preflight (read-only default cache); no evaluation fitting; preserved/frozen. |
|`ce04-evaluation-v2`|All30 task receipts verified; final summary/accounting and archival SHA manifest complete/checked. |

All archive shorthand is under `validation/consensus-ebmf/`. Helpers: `benchmarks/consensus-ebmf-seed-smoke.R`, `benchmarks/consensus-ebmf-evaluate.R`, `benchmarks/consensus-ebmf-callback.R`, `benchmarks/consensus-ebmf-generators.R`. CE02/CE03 native/candidate functions are extracted from frozen source via parse into dedicated environments; do not source whole historical runners and trigger old fits.

## 7. After G1: concrete G2 review and remaining limits

Concrete proposal prepared2026-10-09: `docs/plans/2026-10-09-consensus-ebmf-g2-ceiling-proposal.md`, mirrored in `validation/consensus-ebmf/g2-design-20261009/`. Exact one-sided Clopper-Pearson ceiling check:60 fresh independent datasets; zero baseline shortfall events gives4.870291% upper probability, bounding both candidates' maximum possible mean gain below5pp. One event makes this screen INCONCLUSIVE. Other regimes/real data have zero new runs. Fixed scope retains allfive arms/seven sensitivities; only the single negative recovery-ceiling claim is inferential. Named reservations548.5333CPU-hours including180folds,2retry and2infrastructure slots; requested ceiling550. Pilot allocation extrapolation54.5CPU-hours from onlyfive fresh primary folds, not a guarantee. Planning self-check passed; full GO design remains unresolved. Existing G1 evaluator/guard are pilot/100-hour-specific and must not be reused unchanged or edited in place. Approval would authorize separately versioned excluded orchestration; no package changes/publication included.

Read the final G1 report, `docs/plans/2026-10-09-consensus-ebmf-g2-design-review.md` and live tracker `docs/plans/2026-10-08-consensus-ebmf-plan.md`. Review the primary recovery ceiling first:best155 already recovered100% in both pilot datasets, so no observed5pp benefit exists. This does not prove improvement impossible in the population. Two pilot datasets/regime are insufficient for the registered scientific precision. Primary freeze still needs:complete arm/endpoint/regime simultaneous-bound family;nonzero variance safeguard/minimum n instead of unreliable1.2×two-trialSD;regime n caps/operating power;exact/bootstrap method;zero-baseline/missing-loss and failure-rate rules;full fresh-fold cost and explicit CPU budget. Final fresh-fold timing is in the report. Do not promote descriptive plug-in n to execution authorization or use folds/cells as independent trials. Any revised generator needs separately versioned development/calibration evidence before primary freeze.

Prospective margins remain recovery improvement>=5pp,relative prediction degradation<=5%,FDP increase<=2pp,plus numerical-failure and real-donor prediction evidence. No scientific verdict yet. G3 dataset candidate/preprocessing is discussed in tracker (Kang18_8vs8/EH2259); no real-donor split/preprocessing/evaluation protocol is frozen; cached or timing data is not primary input authority. Pancreas technology transfer cannot supply missing donors; never invent donor IDs.

No public consensus/GBCD/STARCAT APIs, new release, publication, primary synthetic run or real-data comparison have been performed. Suggested next-session skills:ADHD default;`bioinformatics-scientific-workflow`/`snakemake-workflow-engine` for the active pilot;`review-clear-execute` only when reviewing a concrete next implementation plan. User instructions override skill defaults. Do not update Codex memory unless explicitly requested; updating this handoff is not a request to edit memory.


## 8. Current G2 execution — approved, running (2026-10-10)

Approval: user `approve`, then `resume`, for the exact proposal at
`docs/plans/2026-10-09-consensus-ebmf-g2-ceiling-proposal.md` (durable original
`validation/consensus-ebmf/g2-design-20261009/proposal.md`). Approval receipt binds
its SHA256. No G3/G4/full benefit/publication scope follows.

Current study root: `validation/consensus-ebmf/ce04-g2-ceiling-v2/`.
Failed infrastructure-only v1 root is sealed with SHA256SUMS; do not overwrite it.
The only version change corrects comparison of serialized native fits to compare
L_pm/F_pm/ELBO/PVE instead of entire objects with function environments. All primary
inputs and initializer maps are byte-identical between versions; no primary fit
ran under v1. Original fitting, generator, thresholds and seven sensitivities remain
unchanged. All scientific evaluator helper bodies/formals match G1.

Exact submitted IDs at this snapshot:
- v1 infra-create25742287: FAILED after13allocatedCPU-seconds; checkpoint exists,
  G1 software fixture oracles passed; reservation1830CPU-seconds retained.
- v2 infra-resume25742297: COMPLETED12allocatedCPU-seconds; PASS receipt,
  zero native compute calls, numerical checkpoint identity/configuration refusal.
- v2 first primary `identifiable_donor_heterogeneity-t1-f1`25742299: RUNNING.

Guard/manifest:184 total named attempts across versions:1prior failed infrastructure
slot +1v2infra slot +180original folds +2unused retry slots. All named reservations
548.5333CPU-hours <=550. The v2 guard mandates1830prior reservedCPU-seconds and
retains all new reservations. Current prior+all180submitted reservations1953060seconds
(542.516667CPU-hours); all allocations are charged across versions. Live quota for medium
is restricted (8G/4CPUs), so two4Gtasks may run concurrently; other user jobs remain
untouched. About54.5allocatedCPU-hours was extrapolated from onlyfiveG1fresh primary
folds, not a runtime guarantee; queue/wall time may exceed a day.

Foreground resumable automatic controller: **exec session44555**, command:

```bash
/usr/bin/python3.11 validation/consensus-ebmf/ce04-g2-ceiling-v2/launcher/run-approved-study.py 25742299
```

Its external log is `validation/consensus-ebmf/g2-run-20261010T0715Z.log`.
Controller waits for first job COMPLETED and its full verification/receipt hash,
then held-submit/inspect/releases only remaining unreserved original slots and
calls `benchmarks/finish-consensus-ebmf-g2.py`. Per-job Snakemake runs both evaluate
and verify inside the allocated job. Controller/worker locks prevent duplicates.
No retry or unrelated scheduler action is automatic. The controller is additional
orchestration with its own hash/self-check under `launcher/`; it does not mutate
frozen method inputs/source. No automatic chat notification is promised.

First resume action: inspect session44555, the log and completion marker. If the
controller is alive, do not restart it. If it ended, read its exact error before
retrying orchestration. Its resume logic skips accepted/released originals and
refuses ambiguous reservations/held jobs. Two retry slots are reserved in the
approved design but no retry mapping/wrapper is prepared; inspect any failure and
preserve all attempts before preparing a separate unchanged-config resume mapping.
The launcher refuses generic retry slots, so nothing retries silently.

```bash
tail -n 12 validation/consensus-ebmf/g2-run-20261010T0715Z.log
ls validation/consensus-ebmf/ce04-g2-ceiling-v2/POSTPROCESS_COMPLETE.json
```

No completion marker is expected yet. At completion, read `summary/README.md`,
`summary/decision.json`, all60`ceiling-events.csv` rows, all-attempt/descriptive
metrics, failures, exact-ID accounting/budget and finalSHA256SUMS. Missing/failed
folds count shortfall=1, never disappear. Any source/protocol/ledger-integrity
problem blocks a verdict. Zero dataset shortfalls bounds gain below5pp; one event
makes this fixed60trial screen INCONCLUSIVE. It cannot yield GO.

New excluded files: `benchmarks/*g2*`, `benchmarks/workflow/consensus-ebmf-g2.smk`,
the frozen v2 source/input/protocol/manifest, the separately hashed controller,
and preflight logs under `validation/consensus-ebmf/g2-preflight-20261010/`.
Seven mocked budget tests, allocation/refusal checks, fixed-denominator/integer-truth
failure checks, unchanged helper parity, primary generator oracles,60300seed audit,
Snakemake dry-run and checkpoint preflight passed. A preliminary summary parse error
and failed v1 native-object identity preflight are retained in logs, not hidden.
No fresh public-package full suite/build/site/CI/release validation was run; public
package source/API and G1 frozen archives remain unchanged.


### G2 progress audit,2026-10-10 12:12UTC

Controller session44555 is alive and monitoring all180submitted primary jobs.
First-fold marker is PASS. Snapshot36verified/2running/142pending;11complete
independent datasets. No primary job failures, ended jobs without receipts, or
method failures across allfive arms in the36verified folds. Read-only independent
R audit checked frozen source/artifact identities,36receipt hashes, exact manifest
and job/task identities, integer family count and all-attempt accounting.
Evidence: `validation/consensus-ebmf/g2-status-20261010T1212Z/` (README,audit.R/log,
exact-ID accounting/queue snapshots,receipt.json,SHA256SUMS). Terminal allocation
29803CPU-seconds/8.278611CPU-hours includes both infrastructure attempts; active
allocation is not final. Reservation1953060seconds/542.516667CPU-hours within550.
Two retry slots remain unused; no retries or method/source/schedule changes.

The controller has passed its first-fold launch stage; it is now in the final
monitor/report worker. Do not launch jobs again or rerun finalization while it is
alive. No completion marker/final scientific verdict yet. Next bounded action:
check controller/log for new completions; if complete, independently audit all60
dataset events, denominator, exact Clopper-Pearson decision and final manifests.
Public package tests/build/site/CI were not rerun during this monitoring audit.
