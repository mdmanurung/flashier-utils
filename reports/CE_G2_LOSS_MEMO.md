# CE G2 loss memo + rerun decision doc (flashier-utils)

Date: 2026-10-16. Context: the CE G2 benchmark artifacts were session-local and
are lost (recorded in HANDOFF.md reality section). This memo records what was
lost, what survives, and the rerun decision with cost.

## What was lost

- CE G2 result tables and the rendered report for the second consensus-EBMF
  evaluation round. The run consumed approximately 550 CPU-hours on the user's
  Slurm cluster; no frozen roots, result RDS files, or logs survive in the repo
  or shared storage.

## What survives

- The benchmark design and any committed harness code under `benchmarks/`
  (consensus-ebmf helpers were session-local and are ALSO lost — see HANDOFF).
  The surviving basis for reconstruction is the HANDOFF description, the
  flashier-gbcd/flashier-starcat trackers, and the package's current tested
  behavior (1230-test suite, FND-001/002 oracles).

## Rerun decision

Recommendation: DO NOT rerun CE G2 as-is. Reasons:
1. The original harness code is lost; reconstruction is itself a multi-day
   engineering task with no registered rule to gate the result.
2. Stage A of cellprograms has since superseded the question CE G2 addressed
   (consensus-EBMF operating point) with a registered, benchmarked operating
   point (model_selection factorial, ~324 CPU-h, RESULTS.md committed).
3. The same compute budget is better spent on the preregistered Stage B arms
   (batch-adjustment benchmark, anchored refinement) which have live decision
   rules.

If the user wants CE G2 reconstructed anyway: estimate ~40-80 CPU-h to rebuild
the harness + ~550 CPU-h for the run, requiring Slurm access (unreachable from
the sandbox — staged scripts only). Decision owner: user. Default: track as
backlog, revisit only if a Stage B result contradicts the consensus-EBMF
assumptions.
