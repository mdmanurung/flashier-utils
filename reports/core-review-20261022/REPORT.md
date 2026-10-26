# Core review of flashier.utils

Baseline: `63f197843cd3115c072fbc8038bf4027f69a141f`, initially clean `main`.
Runtime: R 4.5.1 through `scripts/run-r.sh`; flashier 1.0.60 at
`fd44811ab22d61ee905e1964175571c452fe7683`; ebnm 1.1-42 at
`825edf1e308af9f3e594a5482300b9fb99793129`.

## Delivery status

**Both approved API changes are implemented and qualified in development version
2.0.0.9000.** The user approved the export reduction and all nine signature
reorders on 2026-10-23 ("approve both"). The package exposes 21 public functions;
34 retained implementations and manuals are internal. The site promotes the 21
public functions and four tutorials. Scientific calculations and defaults were
not changed by this migration.

Source and final installed suites each pass **1,475 assertions, 0 failures/errors/
skips, 11 known native projection iteration warnings**. The final tarball passes
`R CMD check --no-manual`: **Status OK**, including public examples and all four
rebuilt vignettes with cached pancreas input. The direct-native matrix passes
**514/514**. Pre/post migration and source/isolated-installed smoke outputs are
exactly identical for activity, effects, offsets, factor IDs and coordinate manifests.

Final evidence: `validation/core-api-migration-20261023/integrity-receipt.json`.
Final tarball: `validation/core-api-migration-20261023/flashier.utils_2.0.0.9000.tar.gz`.
The receipt binds 149 package source/staging hashes, the tarball, full tests,
native/caller checks, isolated installation, site checks and preserved attribution.
It verifies that all 120 function bodies and defaults match the pre-migration
source, and all 111 existing direct caller bindings retain their meaning. Forwarded
dots are omitted from that static binding comparison and exercised by regression tests.

The first complete 2.0 candidate/check is preserved under `first-candidate/`.
Its only difference from the final source is a README changelog link that the site
checker caught. The final tarball was rebuilt, fully rechecked and reinstalled.
Earlier 1.1 qualification remains at `validation/core-review-20261023-final/`.

The package changes are committed in 219375d, e8046f1 and dcb5c47 and pushed to
origin/main. Live git ls-remote verified dcb5c471f1624035f6e11e095927f10fd544f7db.
Review reports and validation archives remain local. No release, tag rewriting,
scheduler study or external reviewer call was performed. PDF manual was not checked locally. Hosted R CMD check and pkgdown subsequently
passed for the pushed SHA; see the delivery receipt. Native iteration warnings do not establish convergence.

## Reproduced defects and repairs

References identify the repaired implementation; original behavior is recoverable
at the baseline commit. Regression cases are in `tests/testthat/test-core-review.R`
unless indicated otherwise.

| Severity | Defect | Repair / evidence |
| --- | --- | --- |
| High | `standardize_factors()` could not declare bare native-fit measurement units, so comparisons and projections rejected an otherwise valid fit. | Append validated `feature_scale`, `preprocessing`, `independent_unit` declarations; refuse conflicts with known recorded declarations. Native objects remain unmodified. `R/representation.R`, both orientations tested. |
| High | Held-out preprocessing could reorder scoring features after activities were inferred, silently scoring against misaligned truth. | Validate scoring dimensions and IDs just as for inference; record a fold failure instead of metrics. `R/holdout.R`. |
| High | limma silently ignored formula offsets while lm honored them, changing coefficients by exactly the supplied offset. | Refuse offsets for non-lm engines until a native equivalent is explicitly qualified. `R/association.R`. |
| High | Grouping columns named `basis_id`, `factor`, or summary-statistic names could overwrite labels or mix programs during centering. | Reject reserved output-name collisions before grouping. `R/metadata.R`. |
| Medium | `fix_factors()` required orientation metadata but offered no argument for bare native fits. | Append `sample_side`, resolve it explicitly, retain it on the returned fit; compare with public native fixing on both axes/orientations. `R/operations.R`. |
| Medium | A zero-factor reference crashed program alignment, although a zero-factor target worked. | Construct genuinely empty typed match columns; test empty reference/target/both. `R/alignment.R`. |
| Medium | Draw-based backprojection silently ignored requested term/contrast/condition selectors when their columns were absent. | Reject missing selector columns consistently with coefficient mode. `R/backprojection.R`. |
| Medium | `validate_backprojection(feature_scale=NULL)` bypassed the declared-unit boundary. | Reuse strict measurement-scale validation for NULL, NA, empty, unknown, numeric and vector inputs. `R/validation.R`. |
| Low | An empty global factor selection emitted max warnings and an opaque row-count error. | Fail early with a factor-selection error. `R/features.R`. |
| Low | Metadata ranking parsed nonsyntactic column names as formula code. | Build a formula from literal name symbols; verify the same result as a simply named column. `R/metadata.R`. |

An additional late compatibility repair lets `factor_stability()` callbacks return
an explicitly declared `factor_representation` from a directly fitted native model.
Both orientations pass focused regressions, without users editing private attributes.
The final full-package qualification includes this repair and its six assertions.

Verification repairs:

- The old rotation fixture used the wrong inverse transformation on effects;
  correct it to the inverse transpose and assert reconstruction equality
  (`tests/testthat/test-fu1-stability.R`).
- The baseline suite had 1 failure, 1,390 passes and 11 warnings. Its failure
  compared extraction provenance timestamps across calls. Test input immutability
  directly instead (`tests/testthat/test-fu2a-fu2c.R`).
- Replace site checker constants (34 exports / 13 articles) with current
  export/registry/manual/vignette inventories, enforce at most four tutorials,
  and check stale links, attribution and local links. Six independent fixtures
  exercise the checker (`scripts/test-check-site.py`).
- Exclude tracked scientific `reports/` from package builds, as already required
  for local archives. `inst/NOTICE`, licenses and `inst/provenance.csv` are retained.

## API inventory and applied migration

`exports.csv` preserves the 55 pre-migration exports, their task, source, origin, upstream
symbol, source callers, test references, vignette references and proposed disposition.
`arguments.csv` records every current formal argument and its position.
`PROPOSAL.md` lists the approved 21 retained functions and every internalization.
These inventories report source references, not a claim of complete branch coverage.

The recommended core consists of native fitting/standardization, matched extraction,
feature selection and two essential plots, grouped summaries/contrasts and basis-safe
association/backprojection, uncertainty, replication/stability, held-out evaluation,
projection, and provenance. No new comparison export is needed: existing tools plus
one concise case study cover the demonstrated use case.

Specialist diagnostics/gates/networks, Seurat/UMAP, generic statistical adapters,
and low-level ranking/display helpers are candidates for internalization. Native
editing and specialized prior fitting can be shown through public flashier calls.
Do not delete internal implementations still called by core functions or discard
upstream notices. Export removal requires migration documentation and updating
namespace/provenance tests, examples, references and pkgdown together.

Applied argument convention:

1. Primary fit/data first; comparisons use reference/target as their first pair.
2. Metadata and required scientific declarations precede optional controls.
3. Use consistent factor, sample/feature, group, contrast and ID selectors where
   concepts match; singular `factor` remains valid for genuinely one-factor tasks.
4. Keep computation/display controls last. Use named scientific arguments in
   examples. Do not silently reinterpret old positional calls.

The approved table in proposed-arguments.csv now matches all 21 public signatures.
The nine reorders align sample-side/representation declarations, feature/plot
selection, explicit seeds, biological units and measurement units. No scientific
default or return calculation changed. Named external calls retain their meaning;
old positional calls and removed package exports can break. MIGRATION.md records
the external-caller transition. Retained internal manuals keep their examples
conditional rather than presenting them as supported public workflows.

## Native compatibility and fidelity

Entry point: `standardize_factors(native_fit, sample_side=..., feature_scale=...)`.
Most interpretation consumers accept its matched representation. Keep the original
native fit and pass orientation plus representation explicitly for posterior
uncertainty/sampling. Bare callback fits in holdout need explicit projection priors;
stability refits need a callback when a replayable wrapper recipe is absent.

| Area | Evidence / interpretation |
| --- | --- |
| Native numeric parity | `validation/core-api-migration-20261023/native-matrix.csv`: 514 assertions passed across 28 native fixtures and two expected unsupported-input cases. |
| Input combinations | Both orientations; zero/one/multiple factors; signed/semi/nonnegative support; supported dense/sparse/noise inputs; fixing/intercepts/reordering/removal/serialization. |
| Coordinates | Raw/LDF/RMS reconstructions versus native fitted values; offsets included where declared; posterior SD/LFSR transformations, PVE and direct seeded sampler summaries. |
| Fitting and editing wrappers | Existing direct-native oracles for ordinary fitting, initialization, multistart, fixing/refitting and projection pass in the complete suite. No fitting algorithm was replaced. |
| Independent utilities | Hand calculations, explicit information-boundary failures, literal-name/grouping fixtures, and the repaired rotation oracle; native parity is not claimed for nonexistent upstream utilities. |
| Installed identity | Isolated install of the built package; identical activity/effect/offset/factor/manifest outputs against the source-loaded smoke fixture. |

Limits: one pinned flashier revision on Linux/R 4.5.1, not all past/future versions.
Unnamed fitted objects are intentionally rejected. Bare native intercepts reconstruct
correctly but remain ordinary factors without role metadata. Twenty-four initialized
matrix fixtures reach their bounded 100-iteration limit: numeric parity at attained
fits is established, convergence is not. The suite's 11 native projection iteration
warnings remain visible.

Version-sensitive object fields remain documented in
`validation/core-review-20261022/native-field-dependencies.md`. Public mean accessors
match values but lose dimnames for greedy fits; a blind replacement would lose IDs.
PVE, LFSR and sampler output have no dedicated public accessors in the pinned version.
This dependency is a compatibility boundary, not evidence of numerical drift.

## Different priors: prediction is distinct from program recovery

The reproducible local study runs three truth regimes × two data seeds × two prior
settings: 12 full-data fits plus 24 donor-fold masked refits. All completed without
recorded model/fold failures or warning conditions. Native numerical-zero messages
are retained. It records ranks, runtime, objectives, failures/warnings, training and
known-signal metrics, shared-fold prediction, program matching, support diagnostics,
and SVD subspace distance. Native objectives are not used to rank different priors.

| Example | Observation | Interpretation |
| --- | --- | --- |
| Signed sparse | Normal/normal signal RMSE .1087/.1156; normal/point-Laplace .0928/.1119. Masked RMSE ranges overlap. | A small illustrative contrast, not calibrated superiority. |
| Nonnegative | Constrained signal RMSE .0934/.0999 versus normal .1075/.1202; constrained masked RMSE can reach .5203, normal .3747. | Better full-data signal reconstruction did not guarantee better held-out prediction; some constrained fold fits dropped rank. |
| Correlated/overlapping weak programs | Both settings match only 1 of 3 truth programs at .8 congruence despite masked R² .7889–.9238. | Good prediction does not establish complete factor recovery. |
| Exact rotation | Signal difference about 1.0e-17; effect-subspace distance about 2.6e-16; pairwise congruences .707, below .8. | Individual-factor disagreement can preserve the same reconstructed signal. |
| Unequal-rank split | Two truth programs versus three target programs; exact same signal, one duplicate program unmatched. | One-to-one unmatched/spurious counts need split/merge interpretation and the complete similarity matrix. |

Two seeds, identity preprocessing, zero offsets and .3 Gaussian noise make this an
illustration, not calibration. Truth recovery is available only for simulations;
real-data agreement/transfer does not establish biological conservation or truth.
See `prior-comparison-runs.csv`, `prior-comparison-masked.csv`, split matrices and
`prior-comparison.rds` in the validation archive for metrics and loaded-source hashes.

## Documentation consolidation

19 tutorials / 3,209 lines become four / 379 lines (about 88% fewer source lines):

| Tutorial | Purpose |
| --- | --- |
| first-fit | Direct native GTEx fit, extraction, defining features, one plot. |
| advanced-use | Explicit simulated data, coordinate preservation, conditional uncertainty, one contrast/backprojection. |
| compare-priors | Common data/folds, reconstruction versus masked prediction, matching and subspaces. |
| pancreas-programs | One attributed published-data interpretation/transfer example. |

`vignette-dispositions.json` records every original heading/chunk and disposition.
Full original sources are archived under `validation/core-review-20261022/original-vignettes/`.
The removed credits vignette body exactly duplicated retained `CREDITS.md`; authors,
licenses and notices are preserved, and all current references link to the full credits.
README/navigation were updated; full argument tours remain in function references.
The real-data example executes with the retained pinned pancreas cache; it makes no
donor-level inference claims. Package builds never download that data.

## Validation receipts

- Final source and installed suites: **1,475 passed each, 0 failed/errors/skips; 11 native warnings each**.
- Independent native matrix: **514/514 passed**, warnings and unsupported cases recorded.
- Site checker fixtures: **6 passed**.
- Build with all four vignettes and real cached pancreas input: **passed**.
- Isolated installation and exact source/installed smoke comparison: **passed**.
- R CMD check: **Status OK**; generated-site/link validation: **21 exports / four articles pass**, with current input hashes verified.
- PDF manual and cross-platform/current-upstream CI are not performed locally.

All executable study scripts, receipts, failed attempts, runtime information and
source hashes are in `validation/core-review-20261022/`, excluded from package builds.
A source snapshot and tarball manifest bind the tested package. The approved API migration is complete; any commit/push/release remains a separate delivery action.


## Final site qualification (2.0.0.9000)

The new site passes checks for all 21 public references, four article pages,
attribution, source/manual coverage and local links; six checker fixtures pass.
Internal manuals/pages retain authorship but are absent from promoted navigation.
All 55 provenance records, CREDITS.md, NOTICE and LICENSE are preserved. The
provenance generator now retains audited rows added or corrected after the old
plan; a disposable regeneration preserves all 55 records byte-for-byte.

The local site uses integrity-verified cached CDN assets and the supported
`home.sidebar=FALSE`, `news.cran_dates=FALSE` overrides because external metadata
lookups are unavailable. External CRAN metadata validation and deployment are not
claimed. See site-receipt.json, site-build.log, site-home-refresh.log and the asset
integrity CSV in validation/core-api-migration-20261023/.

## Historical pre-migration closeout (2026-10-23)

The source, staging and final tarball agree on all tested implementation/manual/test
payloads. Existing native 514/514 and site input hashes match current source; the
site checker was rerun successfully against the retained site. All 55 public exports,
LICENSE, inst/NOTICE and inst/provenance.csv remain unchanged from baseline.
Final archive: validation/core-review-20261023-final/.

Two verification-launch/reporting mistakes were corrected without changing package
code: R CMD check needed its output directory created first (check-launch-attempt.log);
CSV serialization needed the nested test result column omitted (the complete saved
RDS is retained, and test-summary.log/final-tests-counts.csv recover all counts).
The initial strict byte comparison of DESCRIPTION correctly detected build-generated
metadata/reformatting; verify-integrity.py explicitly validates every original DCF
field and the exact allowed added fields. No test rerun was needed for CSV recovery.

No commit, push, release, API pruning or new scientific study was performed.
At this historical checkpoint the API decision was pending. The user subsequently
approved both changes; the current migration and receipt are recorded above.
PDF manual and cross-platform CI remain unchecked. Native iteration warnings do not
establish convergence; the bounded prior comparison remains illustrative.

## Commit and push delivery (2026-10-23)

Three focused commits deliver the correctness repairs, tutorial consolidation and
approved API migration. All 149 qualified package inputs match the final committed
tree byte-for-byte. The tracked tree is clean; reports/core-review-20261022/ remains
untracked as a local review archive. Delivery receipt:
validation/core-api-migration-20261023/delivery/receipt.json.
Hosted R CMD check (39937386891) and pkgdown (39937387117) completed successfully
for dcb5c471f1624035f6e11e095927f10fd544f7db. The connector returned an empty list;
the direct GitHub CLI query confirmed both successful workflows.
