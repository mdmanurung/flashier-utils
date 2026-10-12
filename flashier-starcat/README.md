# flashier-utils Consensus Reference Development Bundle — Revised

This bundle is the implementation-ready plan for adding cross-dataset consensus factor references and fixed-reference EBMF projection to `flashier-utils`.

It supersedes the first bundle.

## Main changes from v1

- separates inferential factor coordinates from cross-dataset comparison profiles;
- uses `factor_reference` rather than `flash_reference`;
- specifies STARCAT's feature/profile/graph logic more precisely;
- explicitly prevents LSAP assignment from pruning consensus graph edges;
- defines a common-feature reference policy and zero-offset schema-v1 contract;
- removes a production NNLS engine from core scope;
- moves frozen prior transport to an experimental research track;
- renames "novel program" discovery to **residual-program discovery** and gates it behind native `flashier` feasibility tests;
- replaces additive known/novel PVE language with reconstruction-improvement metrics;
- reduces the initial public API surface;
- adds a source audit, feasibility probes, deterministic-ID contract and tracker validation.

## Files

- `00_REVIEW_NOTES.md` — iterative review findings and corrections.
- `01_IMPLEMENTATION_SPEC.md` — revised architecture and implementation requirements.
- `02_MILESTONES.md` — dependency-ordered core, conditional and experimental milestones.
- `03_EXECUTION_TRACKER.md` — readable task backlog.
- `03_EXECUTION_TRACKER.csv` — machine-readable task backlog.
- `04_API_AND_DATA_CONTRACTS.md` — normative schemas and error policies.
- `05_TEST_BENCHMARK_PLAN.md` — exact invariants and scientific validation plan.
- `06_FLASHIER_FEASIBILITY_PROBES.md` — go/no-go probes for residual programs and experimental priors.
- `07_SOURCE_AUDIT_AND_DECISIONS.md` — audited sources and explicit architectural decisions.
- `08_AGENT_HANDOFF.md` — instructions for coding agents.
- `09_BUNDLE_VALIDATION.md` — static consistency checks for the final bundle.
- `milestones.json` — machine-readable release tracks and milestone dependencies.
- `BUNDLE_MANIFEST.json` — file hashes, task counts and audit snapshot.

## Core path

```text
M0 -> M1 -> M2 -> M3 -> M4 -> M5 -> M6 -> M10 -> M12
```

## Conditional residual-program path

```text
M6 -> M7 -> GO -> M8 -> residual validation/docs
```

## Experimental prior track

```text
M6 -> M11
```

The experimental prior track does not block the core release.
