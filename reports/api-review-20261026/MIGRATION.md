# Migration: 2026-10-26 argument cleanup (2.0.0.9000)

Named calls that never used a removed argument keep their meaning. Drop any
removed argument from your call; each accepted only its default.

| Old call | New call |
|---|---|
| `factor_replicability(r, t, alignment = "congruence", control = list())` | `factor_replicability(r, t)` |
| `factor_stability(X, fit, md, "donor", 7, 100, "unit", ...)` | `factor_stability(X, fit, md, "donor", 7, 100, ...)` |
| `factor_stability(..., resample = "unit", alignment = "congruence")` | `factor_stability(...)` |
| `factor_holdout(..., freeze_noise = FALSE, workers = 1L)` | `factor_holdout(...)` |
| `test_factors(..., factor_uncertainty = "ignore")` | `test_factors(...)` |
| `backproject_contrast(..., uncertainty = "conditional")` | `backproject_contrast(...)` |
| `plot_factor_activity(..., engine = "ggplot2")` | `plot_factor_activity(...)` |
| `plot_factor_features(..., engine = "ggplot2")` | `plot_factor_features(...)` |

Positional calls past these positions change binding:
- `factor_stability()` position 7 is now `scheme`, then `fraction`, `strata`.
- `fit_ebmf()` positions 14 and 15 are now `intercept`, `control`.
- `project_factor_activity()` position 8 is now `prior_policy`.
- `plot_factor_activity()` position 7 is now `max_samples`.
- `plot_factor_features()` arguments after `style` are now positional; `...` is last.

Name these arguments in calls. Exact orders: `proposed-arguments.csv`.
