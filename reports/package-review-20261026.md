# Package review — scientific correctness and programming (2026-10-26)

Scope: all 29 files under R/ read in full; probes run with scripts/run-r.sh;
full testthat suite rerun (164 test blocks, 0 failures, 11 known
"Maximum number of iterations reached" warnings in test-projection.R).

## Verdict
No correctness bug found in any traced path. Main risks are defaults and one
missing capability.

## Scientific findings (ranked)
1. Projection re-estimates the activity prior g and residual variance on the
   query (`project_factor_activity()`, `.native_projection()`); documented at
   R/projection.R:238. Probe: same cell in two query sets differs ~1e-3.
   Recommend `prior_policy = "fixed_reference"` using `fit$L_ghat`/`F_ghat`
   with `ebnm_*(g_init = g, fix_g = TRUE)`; makes holdout "fixed_programs" truly fixed.
2. Residual-variance defaults disagree: `project_factor_activity()` var_type=0,
   `project_factor_features()` feature-wise (hard-coded), `factor_holdout()`
   inherits training fit. Recommend activity projection default to reference
   `configuration$var_type`.
3. `factor_stability()` resamples donors with replacement; duplicated rows bias
   EB prior/noise estimates. Expose `resample = "subsample"` (fraction arg).
4. `test_factors(engine="limma")` equals lm for small K (probe: identical SE/df/p
   at K=2). Flag exists; document that moderation needs K in the tens.
5. Minor: `min_congruence` strict `>` (R/alignment.R:36,42) vs "acceptance
   threshold" wording; fgsea on point-mass effects ranks many exact zeros;
   variational posterior SDs/intervals are mean-field underdispersed (add a sentence).

## Programming findings
- Bare-fit default coordinates are LDF: `factor_activity(fit) != fit$L_pm`;
  `factor_uncertainty(fit)`/`factor_intervals(fit)` likewise. Say so in function docs.
- `.metadata_activity()` never passes `extra="drop"`: cohort-wide metadata tables
  fail in every grouping/testing function.
- `factor_intervals()` calls `fit$sampler(64)` per batch, drawing full L and F;
  ~1.5 GB/batch at 100k cells x 30 factors. 256 MiB guard counts retained values only.
- `basis_id` hashes `.engine_metadata()` (R/representation.R:306-317), contradicting
  R/runtime.R:1 comment; harmless (revalidation uses stored manifest).
- Dense no-space formatting in projection.R, programs.R, program-markers.R,
  program-diagnostics.R, seurat.R; run styler once.
- `project_factor_features()` uses `identical()` for fixed values,
  `.native_projection()` uses `all.equal(tol=1e-12)`; `.projection_S()` skips
  positivity check for matrix S.

## Verified correct
LDF/raw/rms reconstruction and multiplier recovery; anchor sign into intervals
and prob_positive; lm contrast SE diag(C'VC) and t p-values; back-projection
diag(B Sigma B'); adaptive nullcheck indices with protected intercept; sample-side
prior reversal in fit_ebmf/fit_seminonnegative_ebmf/fit_initialized_ebmf; masked
holdout information boundary; ImmGenT membership filter equivalence.
