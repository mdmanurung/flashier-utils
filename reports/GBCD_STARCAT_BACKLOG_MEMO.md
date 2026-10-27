# GBCD/STARCAT backlog memo (flashier-utils)

Date: 2026-10-16. The full GBCD and STARCAT builds remain descoped and are
tracked in-repo: `flashier-gbcd/` and `flashier-starcat/` (trackers at 0%).
This memo anchors the backlog and states the entry conditions.

**Reconciled 2026-10-26.** The "0%" trackers were stale, not a record of zero
progress. GBCD-related fitting shipped under the audited plan's GB items:
semi-nonnegative, intercept, initialized and multistart fitting (GB-001/002/003/006,
2026-10-13/14) and the optional exact-pin GBCD adapter `fit_gbcd_programs()`
(GB-004, 70a30d6, 2026-10-23). GB-005b display ordering was completed on
2026-10-27 (canonical MDS sign, declared `group_order`, dense-distance refusal). STARCAT remains unstarted apart from the
`fixed_reference` projection policy (3922cd4). No release GO is recorded.

## Scope summary

- **GBCD** (`flashier-gbcd/`): gene-batch-condition decomposition build per the
  in-repo spec — GB-00x item series. GB-005 (activity prevalence threshold +
  bounded within-group MDS ordering in `summarize_factor_activity()` /
  `plot_factor_activity()`) is DELIVERED and tested (2026-10-16). The
  remaining GB items are unstarted. (Superseded 2026-10-26: see the note above.)
- **STARCAT** (`flashier-starcat/`): sample-side catalog build. Unstarted.

## Entry conditions (what unblocks the build)

1. flashier-utils release GO (see reports/GO_MEMO_SCI004.md) — the builds
   extend released primitives; building on an unreleased base risks rework.
2. A registered rule for the first GBCD milestone (repo law: no R/ change
   without a pre-registered rule).
3. Slurm access for the main runs (sandbox cannot execute them; staged scripts
   + cost estimates only).

## Cost posture

Prior full-build estimates were "weeks of gated work" (tracker note). No
re-estimate is made here — the trackers own the numbers once work starts.

## Disposition

Tracked backlog; not gating the flashier-utils release; revisit after Stage B
of cellprograms produces its first real-data operating experience (the
immune-aging PBMC analysis is the first candidate consumer).
