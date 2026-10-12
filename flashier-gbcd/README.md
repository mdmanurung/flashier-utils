# flashier-utils extension plan — v3 adversarially reviewed bundle

This bundle supersedes v1 and v2.

It is based on a fresh audit of the current `flashier-utils` code and the
upstream APIs actually pinned/used by the package. The objective is not to add
the largest number of features. It is to add semi-NMF, initialized fitting,
GBCD interoperability, multistart diagnostics, and useful display summaries
without weakening the package's existing orientation, representation,
projection, uncertainty, provenance, or resampling guarantees.

## Audited baseline

- `mdmanurung/flashier-utils`
  - commit: `75a1709e9979f66fa46beaee994443b804f9ca64`
- `willwerscheid/flashier`
  - version: `1.0.60`
  - commit: `fd44811ab22d61ee905e1964175571c452fe7683`
- `stephenslab/ebnm`
  - version: `1.1-42`
  - current audited commit: `825edf1e308af9f3e594a5482300b9fb99793129`
- `stephenslab/gbcd`
  - version: `0.2-10`
  - commit: `670e3afb50236303ddf99db33a3394abd4a711bd`
- `stephenslab/gbcd-workflow`
  - commit: `e20227a1633b8b2d9d5660f88f3bfc2d7e9c001b`
- `jhmarcus/drift-workflow`
  - commit: `da05b47f313f183f994155384ca6790f7b865d6e`

## Most important v3 corrections

1. **Projection manifests need repair before intercept support.**
   Projection representations currently inherit source native indices, PVE,
   scaling, and multipliers even though their returned native FLASH fits use a
   new coordinate system. v3 makes projection-native mapping explicit.

2. **Refit safety is broader than the new functions.**
   Edited or projected native fits can currently retain an old fitting
   configuration that does not reproduce the reference object. v3 adds an
   explicit refit policy to existing fit-producing operations.

3. **FLASH intercepts become representation offsets.**
   The fixed all-ones native intercept is removed from semantic factor axes and
   mapped to `factor_representation$offset`.

4. **New-feature projection becomes offset-aware.**
   When a source representation has a feature offset,
   `project_factor_features()` must estimate new-feature offsets rather than
   silently returning zero offsets.

5. **GBCD `Kmax` is not a strict final maximum.**
   The wrapper uses `kmax`, not `max_factors`; upstream final rank can be as high
   as `2*kmax - 1`.

6. **GBCD has mixed activity constraints.**
   The returned `Baseline` originates from a signed point-Laplace component,
   whereas GEP memberships use a nonnegative prior. The global representation
   therefore cannot honestly claim that every activity factor is nonnegative.

7. **Multistart must preserve sign constraints.**
   Replacing the native greedy initializer merely to change its seed otherwise
   loses FLASH's automatically inferred sign constraints.

8. **Initialized-factor nullchecking must preserve IDs.**
   New factors are nullchecked from highest native index downward so user
   supplied factor names remain aligned when factors are removed.

## Files

| File | Purpose |
|---|---|
| `IMPLEMENTATION_SPEC.md` | normative implementation contract |
| `PATCH_MAP.md` | exact package files/functions to change |
| `MILESTONES.md` | dependency-ordered execution plan |
| `TRACKER.md` | implementation checklist |
| `ACCEPTANCE_TESTS.md` | mandatory validation matrix |
| `ARCHITECTURE_DECISIONS.md` | architectural constraints |
| `REVIEW_NOTES.md` | iterative audit findings and corrections |
| `UPSTREAM_AUDIT.md` | source pins, APIs, licensing/reuse notes |

## Recommended implementation sequence

1. Foundation invariants: support metadata, refit policy, projection manifest repair
2. Semi-nonnegative EBMF
3. Native intercept → semantic offset
4. Initialized/anchored EBMF
5. GBCD adapter
6. Group prevalence summaries
7. Multistart optimization sensitivity
8. STRUCTURE-style within-group ordering

GBCD-style split-signed initialization remains a benchmark-only idea until the
simpler native semi-NMF implementation is shown to be inadequate.
