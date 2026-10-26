# Adversarial review of the EBMF–cNMF benchmarks

Date: 2026-10-25. Scope: the local PBMC pilots, saved program comparisons,
matched-rank cNMF extension, and existing stability evidence. This review checks
the comparison design and disposable scripts. It is not a whole-package audit.

**Verdict: the work supports feasibility and conditional descriptions of saved
fits. It does not establish method superiority or novel biological programs.**
One reproducible sign-handling defect affects the constrained loading metric.
Matched-rank results substantially weaken the earlier uniqueness examples.
Original fits and frozen reports remain unchanged.

## New counterevidence from the authorized next step

Two native cNMF K12 consensus runs used five starts each, seeds14/15. They used
the same cells, same 300 genes, same counts, and same variance scaling as EBMF.
Prepared matrices differed from the TSV reference by at most 1.78e-15.

| Method | Actual factors | Relative training SSE, two starts/runs |
|---|---:|---:|
| cNMF K4, original five-start consensus | 4 | 0.742413 |
| cNMF K12, five-start consensus per run | 12 | 0.624363 / 0.621611 |
| EBMF normal scores / point-Laplace effects | 12 | 0.645691 / 0.645691 |
| EBMF point-exponential scores / point-Laplace effects | 12 | 0.655769 / 0.656732 |

These are training errors under different penalties and optimizer procedures.
They cannot rank predictive performance. The cNMF receipts contain unclosed-file
ResourceWarnings; no optimizer iteration-limit warning was caught.

| Previously highlighted EBMF direction | Outside cNMF K4 span | Outside cNMF K12 span, seeds14/15 |
|---|---:|---:|
| Normal-score F4, HLA/CD74-rich contrast | 87.9% | 18.2% / 19.9% |
| Nonnegative-score F9, HLA/CD74-rich contrast | 86.3% | 15.7% / 17.6% |
| Normal-score F3, NKG7-rich versus IL7R/LTB contrast | 79.6% | 27.5% / 9.8% |

The ordinary-span numbers above are unaffected by the sign-handling defect.
They describe squared loading residuals in variance-scaled count coordinates.
They are not percentages of biological expression variance or discovery
probabilities. Remaining differences could reflect constraints, regularization,
rank allocation, optimizer choices, rare cells or noise.

Evidence: [rank sensitivity examples](../validation/ebmf-cnmf-adversarial-review-20261025/logs/rank-sensitivity-examples.csv),
[K12 cNMF fit log](../validation/prior-cnmf-rank12-20261025/logs/fit.log),
[combined loading tables](../validation/prior-cnmf-rank12-20261025/logs/coverage.csv).

## Findings that block stronger scientific claims

### 1. A valid source sign flip changes constrained coverage — confirmed code defect

The comparison chooses both source poles only when scores cross zero. Otherwise
it always chooses the positive loading pole. A signed factor whose scores are
all negative needs the opposite realized pole. Flipping both scores and
loadings is a valid reorientation here and preserves reconstruction exactly.

In the independent synthetic oracle, the reported residual changes from 0 to 1
after that valid flip. In real signed F2, it changes from 0.21628 to 1.0;
the same failure occurs in seed15. Current ordinary-span results remain valid.
The original cone metric is not invariant to a permitted representation change.

**Fix before reuse:** define whether the endpoint means supported capability
or observed directions. For observed directions, enumerate positive and negative
score poles explicitly. For signed capability, assess both poles regardless of
observed signs. Add the sign-flip oracle to the comparison's runnable checks.
The independent corrected observed-pole probe passes. The original disposable
scripts were not silently edited during this review.

Code: `validation/prior-cnmf-greedy12-programs-20261024/analyze.py:58` and inherited
copies. Evidence: [sign-flip probe](../validation/ebmf-cnmf-adversarial-review-20261025/logs/sign-flip-probe.csv).

### 2. Unequal rank generated much of the apparent uniqueness — confirmed design confound

The original comparison gave EBMF 12 loading vectors and cNMF four. A larger
basis can cover smaller-basis directions while having many directions outside
that smaller span. The K12 extension reduces the HLA/CD74 gap from 86–88% to
16–20%. The NKG7/IL7R gap also drops sharply. Treat the earlier examples as
rank-conditional descriptions, not method-specific discoveries.

**Fix:** compare a prespecified matched-rank frontier and recommended native
rank-selection pipelines separately. Do not select the most favorable rank
after examining program labels.

### 3. No method selected its preferred rank

cNMF K4 and K12 were supplied, not selected using its K diagnostics. Every
greedy EBMF run hit Kmax12. Backfit and nullcheck retaining 12 does not show
that the data favor exactly 12. It shows that the ceiling was reached.

**Fix:** freeze the rank-selection rule and ceiling policy first. Report actual
rank, active rank, whether the ceiling binds, and rank sensitivity separately.

### 4. Training reconstruction cannot establish predictive superiority

The same entries supply both fitting and all reported reconstruction errors.
More factors and weaker penalties can improve training loss while fitting
outliers. The newer cNMF SSE is lower than EBMF's at K12, but that does not
establish generalization. No entry-level pseudo-replication or significance test
was performed, which is appropriate at this stage.

**Fix:** discover both models inside paired training folds. Estimate query
scores using permitted entries and evaluate disjoint held-out entries. Include
a simple prediction baseline and keep preprocessing inside the training folds.

### 5. There is no recovery truth or calibrated absent-program endpoint

PBMC gene names and reference marker familiarity do not say which latent
programs should exist in this particular subset. A loading outside another
model's span can be noise. We have no planted shared/private-program benchmark,
null recovery rate, split/merge oracle at the fitted-model level, or independent
validation of the highlighted candidates.

**Fix:** use planted identifiable and deliberately nonidentifiable scenarios.
Include omitted programs, overlapping genes, signed bipolar programs, weak
signals, collinearity and single-cell outliers. Report false candidates as well
as successful recovery. Keep real-data interpretation a separate endpoint.

### 6. cNMF's consensus robustness settings were scaled down or disabled

Both K4 and K12 use five starts and density2.0, so no spectra are excluded.
The linked tutorial uses 20 starts and investigates lower density thresholds;
its text discusses more starts for real data. We also disable final TPM-based
usage refitting. The latter choice preserves a valid count-space reconstruction
pair, but it changes the native workflow being evaluated.

**Fix:** retain the controlled common-matrix arm. Add an explicitly separate,
budgeted native cNMF arm with its diagnostics, sufficient starts and declared
density policy. Preserve paired units if final usage refitting is enabled.

Code: `validation/prior-cnmf-rank12-20261025/fit-cnmf.py:40–49` and upstream
`consensus()` in the frozen cNMF source.

### 7. Initialization agreement is not data stability

Two EBMF seeds are insufficient to characterize local-optimum sensitivity.
Native power-style initialization can converge toward nearly identical
directions despite distinct hashes. Two cNMF consensus seeds do not resolve
that limitation either. The input cells, feature selection and preprocessing
are unchanged in all these comparisons.

**Fix:** record initializer diversity, then vary cells or biological units in
separate refits. Nest multiple starts within each refit. Do not pool seed-by-
resample fits as independent biological replicates.

### 8. Extra directions often concentrate in a few cells

One cell accounts for 61.7% or 67.8% of squared HLA/CD74 scores in the two EBMF
arms. The top five account for 71.9% or 78.9%. Signed F3's top five account for
73.6%. Several other EBMF directions are almost single-cell effects. Original
cNMF factors3/4 are also concentrated; this issue applies to both methods.

**Fix:** use leave-one-cell-out and repeated cell-subsampling diagnostics.
Inspect those cells' QC and sample identity. A rare real state and an outlier
need different interpretations; neither can be inferred from concentration alone.

Evidence: `validation/prior-cnmf-greedy12-programs-20261024/logs/score-concentration.csv`.

### 9. The model arms change several assumptions at once

cNMF requires nonnegative spectra and scores. Both current EBMF arms allow
signed point-Laplace effects, and one also allows signed scores. Priors,
support constraints, regularization, optimizer and consensus all differ.
Equal rank does not isolate a prior effect or give equal model classes.

**Fix:** include an EBMF nonnegative/nonnegative control at the same ranks.
Keep feature priors fixed when studying score-prior changes. Attribute results
to the complete specified model, not to one changed assumption.

### 10. Observation SD1 fixes EBMF's shrinkage environment

The common matrix is scaled to gene variance one, but that is not a measurement
noise estimate. EBMF fixes observation SD1. Count data remain heterogeneous,
and changing noise assumptions can change shrinkage and retained rank.
Unpenalized NMF optimizes a different criterion even when its squared-error
data term matches Gaussian EBMF's data term.

**Fix:** prespecify noise-model sensitivity and inspect residuals. Compare
native objectives only within the same model/settings. Do not rank ELBO against
NMF loss or interpret lower training SSE as better empirical Bayes inference.

Code: `validation/prior-cnmf-greedy12-20261024/run.R:24`.

## Other important limits and implementation risks

### 11. Five-start consensus and individual EBMF fits have unequal compute

Each cNMF output pools five NMF fits and then refits usages. Each EBMF output
is one greedy/backfit pipeline. The current records do not provide an
equal-compute comparison. Merely counting outer seeds conceals that difference.

**Fix:** report candidate-fit count, stages, convergence, CPU time and memory.
Use a separate fixed-compute sensitivity if computational efficiency is a claim.

### 12. Existing timing values cover different stages

cNMF timers include preparation and file work; EBMF timers cover separate
native optimization stages. The two new cNMF runs take about 8.11s and 1.65s,
so cold-start/cache effects alone are substantial. No speed winner follows.

**Fix:** separate end-to-end and optimizer timing. Repeat cold and warm cases
under the same thread policy, environment and input size.

### 13. The sample is too narrow for cohort-level conclusions

The pilot samples 180 cells across library-size deciles from one public PBMC
specimen. It is not cell-type-stratified, does not establish representative rare
states, and contains no independent donor replication. Balanced library-size
sampling does not resolve that limitation.

**Fix:** expand to the full specimen for a technical sensitivity, then use
multi-donor data for biological stability. Cell resampling remains conditional
on this specimen and must be labeled accordingly.

### 14. Only 300 selected genes define the apparent program universe

Filtering and HVG selection use the full pilot before these descriptive
comparisons. That is acceptable for a conditional saved-fit description.
It would leak information if reused unchanged as a held-out preprocessing
pipeline. A program outside these genes cannot be declared absent.

**Fix:** fit selection/scaling within training folds for prediction. Record
feature overlap and gene-universe coverage. Add compatible all-gene profiles
for interpretation without confusing them with reconstruction weights.

### 15. Span coverage does not show that the target used the same program

The projection optimizes combinations of target loading vectors without using
target sample scores. A direction can be representable while the learned
target usages never express it in the relevant cells. Raw pairwise usage
correlations are saved but were not the basis of the earlier uniqueness claims.

**Fix:** report loading geometry and aligned observed usage separately.
Test candidate augmentation on held-out entries with an equal-extra-capacity
control. Improved training fit cannot satisfy that predictive endpoint.

### 16. A large residual has no built-in novelty significance

With 300 dimensions and a four-dimensional target span, an independent
isotropic random direction has expected uncovered fraction 1-4/300, about
98.67%. The executable illustration gives 98.61%. This is not a calibrated
null for correlated genes; it demonstrates why an 80–90% residual is not by
itself evidence of a recovered biological program.

**Fix:** calibrate the chosen endpoint using appropriate synthetic/null controls
and validate candidates independently. Keep exploratory gene examples separate
from a formal program-discovery family.

### 17. Support is inferred from method display names

`name.startswith("signed")` and the corresponding target-name check determine
whether coefficients can be negative. Renaming a model can therefore change
the feasible set without changing its fitted priors. Current fixed names are
consistent, but the prototype is unsafe to generalize by renaming or adding arms.

**Fix:** read declared score/effect support from validated model metadata.
Add a label-renaming invariance check before public API integration.

Code: `validation/prior-cnmf-greedy12-programs-20261024/analyze.py:58–60`.
Evidence: `validation/ebmf-cnmf-adversarial-review-20261025/logs/label-support-probe.csv`.

### 18. Nearest counterparts can reuse one target factor

Independent maximum cosines are not a one-to-one correspondence. They can make
several source directions appear to have the same stable target. The independent
assignment probe enforces one-to-one matching and preserves weak matches rather
than discarding them. For nonnegative-score EBMF, the unmatched direction's
assigned cosine is essentially zero; cNMF K12 also has two very weak assignments.

**Fix:** report all-pairs similarities, explicit one-to-one matching and joint
span coverage together. Do not turn weak matches into a recovered-program count
without a declared rule or handle splits by assignment alone.

### 19. Summary residuals are not explained variance

The plots average loading-direction residuals. They do not weight by realized
component usage, and signed factors with two examined poles contribute twice.
Overlapping components also prevent interpreting component norms as additive
PVE. A high or low mean residual is not global expression reconstruction error.

**Fix:** label the denominator and averaging unit. Keep per-factor span,
per-pole cone, actual usage, component magnitude and prediction in separate
tables. Do not sum component norms as explained variance.

### 20. Gene-weight interpretation depends on preprocessing units

The comparison uses variance-scaled count-space weights. These are distinct
from cNMF TPM spectra and regression-based gene scores. Rare genes and large
library sizes can have strong effects in these coordinates. Listing top genes
does not validate a pathway, cell identity, dense prior or sparse support.

**Fix:** keep weight units explicit, inspect library/QC associations, and show
positive and negative loadings. Assess posterior inclusion separately if needed;
point-Laplace posterior means need not contain exact zeros. Add independently
supported biology only after technical validation.

### 21. The sequence is exploratory and adapts after outcomes

Normal-score random starts collapsed; positive starts were then added. Priors
and ranks were extended after inspecting earlier results. Those steps were
requested, recorded and retained, which supports an honest exploration. They
do not constitute one unchanged confirmatory benchmark protocol.

**Fix:** retain the history, freeze the next comparative endpoint/budget before
running it, and keep all failed or degenerate fits. Do not exclude inconvenient
starts or choose a model globally by ELBO across different priors.

### 22. Reproducibility and software checks do not establish scientific validity

Artifacts preserve source hashes, input IDs, many warnings, checkpoints and
runtime versions. Input equality and reconstruction checks passed. However,
the R benchmark calls a package-private initializer helper from the checkout;
its source hash and checkout commit are not self-contained in the original
study receipt. Top-level Python requirements are pinned, but the temporary
environment and all code dependencies still matter for replay.

The latest recorded package receipt has 1,577 assertions, zero failures/skips
and 11 known projection warnings. The historical synthetic stability benchmark
has 20 donor-group refits per scenario and a weak collinear factor recovering
in 7/20. Neither proves real PBMC biological stability. The new fit scripts
mostly call native flashier/EBNM; model behavior is not a new flashier-utils
factorization algorithm.

**Fix:** capture checkout/helper hashes and full resolved environments in future
receipts. Keep software correctness, optimizer behavior, sampling stability,
prediction and biological recovery as separate qualification statements.

Evidence: `validation/plot-gallery-20261023/receipt.json`,
`validation/benchmarks/stability/full/metrics.csv`, and `benchmarks/stability.R`.

## Stability evidence after adversarial one-to-one alignment

These are two-start/run comparisons on the same fixed data, not resampling
confidence intervals. Cosines use valid sign reorientation only for fully
signed EBMF factors.

| Model | Minimum assigned cosine | Median assigned cosine | What the evidence supports |
|---|---:|---:|---|
| Normal-score EBMF K12 | 0.999999994 | Approximately 1 | Near-identical loading axes for these two starts |
| Nonnegative-score EBMF K12 | Approximately 0 | Approximately 1 | Most axes recur; one direction lacks a stable counterpart |
| cNMF K12, two five-start consensuses | 0.0525 | 0.9884 | Most axes recur; two assignments are very weak |

No general stability winner follows. The procedures, consensus budgets and
hyperparameters differ. In particular density filtering is disabled for cNMF.
The strongest new cNMF directions outside signed EBMF are also start-sensitive:
seed14 factors7/8 have assigned cross-seed cosines about 0.060/0.052.

## Checks that passed and prevent false findings

- Prepared matrices and IDs match across methods; there is no observed
  mismatched-input defect in this controlled comparison.
- Raw cNMF usage is correctly paired with count-space consensus spectra when
  final TPM usage refitting is disabled. Normalized usage or gene scores were
  not incorrectly used for reconstruction.
- Native checkpoint reconstructions reproduce the reported errors.
- Gene permutation, positive factor rescaling, zero-source and duplicated-basis
  oracles pass. NNLS first-order conditions pass on these inputs.
- Normalized loading bases have condition numbers about 1.40–2.48. The current
  empirical gaps are not explained by a grossly singular target basis.
- Original failures, collapsed starts and capped fits remain visible. Existing
  frozen archives were not overwritten to improve the narrative.

## Smallest defensible next benchmark

1. Repair and qualify comparison metrics first: sign/pole contract, support
   metadata, label invariance and one-to-one versus split/merge reporting.
2. Freeze a small factorial design: common matrix, matched ranks, a fully
   nonnegative EBMF control, declared noise policy and separate native-pipeline
   sensitivity. State compute and convergence policies before outcomes.
3. Add planted shared/private/bipolar/collinear/noise/outlier controls. Evaluate
   both false candidates and missed programs with independent analytic truth.
4. Run paired train-only preprocessing and masked prediction. Repeat starts
   within cell subsamples; use donors as resampling units when donor data exist.
5. Confirm candidate programs using observed usages and independent cells/data.
   Keep gene-level descriptions provisional until that evidence exists.

This review did not launch those larger studies, change the package API, rerun
package-wide checks, commit, push or publish. The authorized K12 cNMF extension
is the only new model fitting performed here.

## GBCD addendum

The user subsequently requested GBCD. A native-default pilot returned 13 active
components at initialization kmax12, after a 90s attempt timed out and an
unchanged-default retry ran for246.359s. It emitted two iteration-limit warnings.
This is a provisional fit, not a converged method benchmark.

Three additional comparison contracts matter: GBCD expects natural-log input;
its kmax is not a final-rank bound; and component activity support remains
unknown after upstream processing. The new pilot respects all three. It uses
shared log-expression/score association profiles for cross-method inspection,
not incompatible native weights or a cross-scale SSE ranking. It does not
reuse the defective constrained-coverage endpoint for GBCD.

GBCD has corresponding HLA/CD74 and NKG7-rich association patterns. For the
normal-score EBMF examples, absolute association-profile cosines are0.821 and
0.934. Their raw activity correlations are0.627 and-0.798, respectively, so
profile resemblance is not identical usage. This further weakens any claim
that those biological programs are uniquely recovered by EBMF.

[GBCD report](../validation/prior-cnmf-gbcd-20261025/REPORT.md).
The dedicated StarCAT catalog/reference pipeline remains planned; the core
projection/alignment primitives and shipped GBCD adapter do not complete it.
[Current StarCAT status](starcat-status-20261025.md).

## Review artifacts and attribution

[Executable probes](../validation/ebmf-cnmf-adversarial-review-20261025/probes.py),
[probe results](../validation/ebmf-cnmf-adversarial-review-20261025/logs/probes.json),
[one-to-one matches](../validation/ebmf-cnmf-adversarial-review-20261025/logs/one-to-one-start-matches.csv),
[review journal](../validation/ebmf-cnmf-adversarial-review-20261025/JOURNAL.md).

Native FLASH/flashier: Jason Willwerscheid, Peter Carbonetto, Wei Wang and
Matthew Stephens. EBNM: Jason Willwerscheid, Matthew Stephens and Peter
Carbonetto. cNMF: Dylan Kotliar and contributors; its MIT notice remains with
copied source. PBMC3k: 10x Genomics. NumPy/SciPy and R Core supply the solvers.
