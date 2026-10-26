# API simplification and argument stress-check (2026-10-26)

Scope: the 24 exports and 8 shared validators (`.validate_matrix`, `.validate_ids`,
`.match_ids`, `.resolve_sample_side`, `.resolve_view`, `.validate_view`,
`.validate_basis`, `.check_engine`). Baseline: main @ d43fe33.

## Verdict
The package refuses almost every bad input with a clear message. 322 probes ran
across six edge axes: shape, class, selection, orientation, metadata and numeric
bounds. Probes and results: `edge-matrix.csv`; harness:
`validation/api-review-20261026/probe.R`.

| Verdict | Probes |
|---|---|
| ok | 142 |
| clear refusal | 164 |
| crash, now a clear refusal | 1 |
| invalid input accepted, now refused | 4 |
| unclear error, now clear | 7 |
| needless refusal, now accepted | 2 |
| reserved argument, now removed | 2 |

No probe returned a wrong numeric result.

Each shared validator was reached through the exports: `.validate_matrix` via
`fit_ebmf()` and projection, `.match_ids` via the metadata functions,
`.validate_basis` via `backproject_contrast()`, and `.resolve_sample_side` and
`.resolve_view` via `standardize_factors()` and every extractor. The ID selector
`.select_ids()` was the only shared helper that needed a change.

## Defects fixed
- **One-sample query crashed projection.** Native flashier drops a length-1
  dimension in its wrapup. `project_factor_activity()` now refuses fewer than two
  query samples.
- **Non-positive `S` was accepted.** Scalar `S <= 0` passed `fit_ebmf()` and both
  projections; a zero matrix passed `fit_ebmf()`. The shared `.check_S()` now
  refuses non-finite or non-positive `S` everywhere.
- **Unknown `var_type` leaked a flashier error.** `fit_ebmf()` now names the
  allowed values.
- **`standardize_factors()` on a view or matrix asked for `sample_side`.** It now
  names the wrong input.
- **Position or logical selections gave an ID-format error.** `.select_ids()` now
  says selections must be character IDs. Factor IDs are coerced to character.
- **`factor_intervals()` required `side`.** It now defaults to activity, like
  `factor_uncertainty()`.

## API simplified (breaking)
Twelve arguments accepted only one value or were reserved. Each only raised an
error, so they are removed. The refusal reasons moved into the function docs.

| Function | Removed |
|---|---|
| `factor_replicability()` | `alignment`, `bootstrap`, `refit_inputs`, `control` |
| `factor_stability()` | `resample`, `alignment` |
| `factor_holdout()` | `freeze_noise`, `workers` |
| `test_factors()` | `factor_uncertainty` |
| `backproject_contrast()` | `uncertainty` |
| `plot_factor_activity()`, `plot_factor_features()` | `engine` |

Trailing arguments now follow the controls-last rule: `fit_ebmf()` puts
`intercept` before `control`; `factor_stability()` puts `scheme, fraction` after
`B`; `project_factor_activity()` puts `prior_policy` after `ebnm_fn`;
`plot_factor_features()` takes `...` last. Exact before/after orders:
`proposed-arguments.csv`. Old-to-new calls: `MIGRATION.md`.

`factor_stability(workers=)` stays because it drives the fork backend. Removing
`alignment` does not change the stability checkpoint hash.

## Checked and left alone
- **Required declarations** (`sample_side`, `feature_scale`, `preprocessing`,
  `unit_col`, `seed`, `engine`, `independence`) are a scientific contract.
- **Locked defaults** (`prior_policy = "query"`, `extra = "drop"`, explicit GoM
  memberships, LDF bare-fit coordinates) are unchanged.
- **`match.arg` vector defaults** were skipped. Forwarded vector defaults would
  change the parameters provenance records.
- **`factor_features(factor)` stays singular** because it selects exactly one factor.
- **Whole `test_factors()` results** are not accepted by `backproject_contrast()`.
  Per-factor models carry no cross-factor covariance; pass `$table`, as the
  vignette does.
- **Open items from package-review-20261026** are already resolved in source:
  strict `min_congruence` wording, small-K limma note, mean-field
  underdispersion and LDF notes.

## Verification
- Full suite: 185 blocks, 0 failures, 0 skips, 11 known projection maxiter warnings.
- Probe rerun: only the 16 intended cells changed; see `edge-matrix.csv`.
- `R CMD check --no-manual`: `validation/api-review-20261026/check.log`.
