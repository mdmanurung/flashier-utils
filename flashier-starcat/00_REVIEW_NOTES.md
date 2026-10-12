# Iterative Review Notes

This revision was produced after three review passes over the first implementation bundle and an audit of the current repositories and public APIs.

## Audit snapshot

- `mdmanurung/flashier-utils` main: `75a1709e9979f66fa46beaee994443b804f9ca64`
- `immunogenomics/starCAT` main: `5399183aa6093ff972bb605d2e33d64ce41b944c`
- `flashier-utils` currently pins `flashier` at `fd44811ab22d61ee905e1964175571c452fe7683` in `DESCRIPTION`.
- STARCAT repository license: MIT.

The implementation agent must refresh these snapshots before coding if repository HEAD has moved materially.

---

# Review pass 1 — architecture and scientific contracts

The first bundle had the right high-level direction but conflated three different objects:

1. the inferential coordinates of a fitted EBMF model;
2. normalized profiles used only to compare programs across datasets;
3. the final reference programs used for query projection.

The revised design separates them explicitly.

### Changes

- Do **not** turn comparison normalization into a new inferential `factor_representation` by default.
- Introduce an internal **comparison-profile layer** derived from existing `factor_representation` objects.
- Keep the original source basis immutable and identified by `basis_id`.
- Use a separate `factor_reference` class for consensus programs.
- Require schema-v1 consensus references to have zero offsets; consensus intercept/baseline estimation is deferred.
- Make feature-universe policies explicit for pairwise comparison and final reference construction.
- Add source-factor eligibility records instead of assuming every fitted factor is usable.

This avoids changing inferential activity units merely to compare factor shapes.

---

# Review pass 2 — STARCAT fidelity and `flashier` feasibility

The first bundle simplified STARCAT's catalog construction too aggressively.

STARCAT's published catalog construction uses:

- normalized GEP spectra;
- Pearson correlation on the union of the two datasets' overdispersed genes;
- cross-dataset edges only;
- reciprocal top-ranked neighbors (top seven for TCAT);
- correlation threshold `R > 0.5`;
- descending-edge greedy merging;
- a two-thirds graph-density rule;
- final consensus spectra formed by element-wise averaging on a common feature universe.

The revised specification therefore distinguishes:

- `mode = "generic"`: suitable for ordinary `flashier` factors;
- `mode = "starcat2025"`: a specialized profile/graph contract requiring the statistics needed to reproduce the published STARCAT normalization logic.

The revised bundle also records a critical implementation rule: **the consensus graph must use all-pairs program similarities, not the one-to-one LSAP assignment used by the existing `factor_replicability()` API.** Existing alignment code can provide feature matching, sign rules and similarity primitives, but its assignment step must not determine graph edges.

Public `flashier` already exposes the main primitives needed for the core projection layer (`flash_factors_init`, `flash_factors_fix`, `flash_backfit`, `flash_greedy`, `flash_fit_get_g`, `flash_ebnm`). Therefore no upstream `flashier` change is required for consensus construction or ordinary fixed-program EBMF projection.

---

# Review pass 3 — scope, validation and agent executability

Several features from the first bundle were demoted from core functionality.

### Frozen reference priors are experimental, not core

A consensus program is assembled from factors learned in multiple source datasets, so it has no unique fitted source activity prior. Moreover, activity distributions depend on source composition. The package should not silently choose or pool source priors.

The revised plan therefore:

- keeps current query-side prior re-estimation as the supported core behavior;
- permits explicit user-supplied projection priors;
- moves learned/frozen reference-prior transport into an experimental milestone;
- adds feasibility probes using `flash_fit_get_g()` and `flash_ebnm(g_init=..., fix_g=TRUE)` before any public API is designed.

### A second production NNLS engine is removed from core scope

STARCAT parity is important for validation, but `flashier-utils` does not need to become another STARCAT implementation. Production projection remains native EBMF. STARCAT/NNLS comparisons belong in external benchmark scripts and archived validation receipts.

### "Novel programs" are renamed residual programs

A factor discovered after fitting a reference is not automatically biologically novel. The revised API is `discover_residual_programs()`. Results are described as query-residual programs until independently replicated or otherwise validated.

### Additive factor PVE is not used

Known-versus-residual contribution is reported using reconstruction quantities such as:

- known-reference fit R2;
- combined fit R2;
- improvement in SSE / fit R2 after adding residual programs.

The bundle does not claim additive factor-wise PVE.

### Public API is smaller

Core new exports are limited to the high-level workflow and accessors. Similarity, graph construction and consensus estimation are internal implementation layers unless a later use case justifies exposing them.

### Projection readiness now includes activity-support homogeneity

The current native projection path initializes the fixed reference programs jointly with one activity-side prior family. Therefore a schema-v1 reference may contain programs with heterogeneous activity support for discovery/inspection, but it is marked `projection_ready=FALSE` unless the selected projection programs share one compatible activity-support contract. Per-factor projection priors are deferred to a separate extension rather than approximated silently.

---

# Net result

The revised architecture is:

```text
independent factor_representation objects
        |
        v
comparison profiles (internal, non-inferential)
        |
        v
all-pairs cross-dataset similarity
        |
        v
reciprocal/density-constrained program graph
        |
        v
factor_reference
        |
        +--> native fixed-program EBMF projection
        |
        +--> optional residual-program discovery after feasibility probe
```

Core release work is intentionally narrower than the first bundle: **consensus reference construction + native EBMF projection first; residual discovery second; prior transport experimental.**
