# Handoff — flashier.utils post-1.1.0

**Date:** 2026-10-06 · **Branch:** main (pushed) · **Status:** 1.1.0 released (tag v1.1.0, GitHub release); dev version 1.1.0.9000

## Goal
Interpret EBMF factors with tutorials on real data. 1.1.0 is out; only optional follow-ups remain.

## Next action
Check CI on the latest main commit: `gh run list --limit 3`.

## State
- Local `rcmdcheck` on 1.1.0.9000: 0 errors, 0 warnings, 0 notes
- Done after release: `factor_gate_alignment(direction = c("high","low"))` (default unchanged, summary gains `direction`), `eligible_above` documented as assay-specific (manual + `vignettes/protein-interpretation.Rmd`, no behaviour change)
- Do not rewrite tags v1.0.0 or v1.1.0

## Locked decisions
- Real data in tutorials where the design supports it (flashier `gtex`, Jamboree pancreas); upstream data never shipped; pancreas chunks gated on `FLASHIER_UTILS_PANCREAS_DIR`
- Keep released contracts (AGENTS.md): `.native_projection` uses `all.equal(tol=1e-12)`; projection metric is `fit_r2` (in-sample)
- `eligible_above` default stays 0.5 (documented, not required) to preserve the public API
- `stability-replicability` keeps the GTEx SNP-gene-pair resampling version

## Dead ends — do not redo
- ImmGen-T real data: upstream `data/` untracked, so `protein-interpretation` stays simulated
- Pancreas data has no donor IDs: never claim donor independence
- Fixing both "loadings" and "factors" in flashier 1.0.60 errors; fix one mode
- Do not split program and projection code across commits (shared helpers)

## Read first
- `R/program-markers.R:146` — `factor_gate_alignment`
- `AGENTS.md` — credit and archive rules

## Open items
- Stale Codex memory (says release pending, coverage below nominal) — correct only if asked
- Deferred by design: SuSiE, brms sampling, CorShrink dosym, sparse plus explicit S — need a design decision
