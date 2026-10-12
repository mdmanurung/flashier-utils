# Source Audit and Design Decisions

This document records what is borrowed conceptually, what already exists in `flashier-utils`, and where the revised design intentionally differs from STARCAT.

---

# 1. Audited sources

## flashier-utils

Repository:

https://github.com/mdmanurung/flashier-utils

Audit snapshot:

```text
75a1709e9979f66fa46beaee994443b804f9ca64
```

Relevant existing files:

```text
R/representation.R
R/alignment.R
R/projection.R
R/holdout.R
R/stability.R
R/fit.R
R/provenance.R
```

Existing capabilities confirmed in the snapshot:

- named `factor_representation` with basis hashing;
- raw/LDF coordinate extraction;
- constraint-aware pairwise congruence alignment;
- fixed-program projection using public `flashier` calls;
- exact verification that fixed feature programs do not change;
- feature mapping, coverage, rank and memory guards;
- biological-unit stability/refitting machinery;
- provenance registry requirements.

---

## STARCAT paper

Kotliar D, Curtis M, Agnew R, et al. Reproducible single-cell annotation of programs underlying T cell subsets, activation states and functions. *Nature Methods*. 2025;22:1964-1980.

DOI:

https://doi.org/10.1038/s41592-025-02793-1

Published catalog-construction elements used conceptually:

1. independently discovered GEPs from multiple reference datasets;
2. normalize program spectra before cross-dataset comparison;
3. Pearson correlation on pair-specific informative genes;
4. cross-dataset edges only;
5. correlation threshold > 0.5;
6. reciprocal top-seven neighbors for TCAT;
7. descending-edge greedy merging;
8. at least two-thirds connectivity for merged sets;
9. common final feature universe;
10. element-wise average to form consensus programs;
11. fixed-reference usage inference in query data.

The paper additionally uses gene-wise SD-normalized counts and row-normalized query usages in STARCAT. Those details are **not automatically imposed** on generic `flashier` factors.

---

## STARCAT repository

Repository:

https://github.com/immunogenomics/starCAT

Audit snapshot:

```text
5399183aa6093ff972bb605d2e33d64ce41b944c
```

License: MIT.

Relevant source:

```text
src/starcat/build_consensus_reference.py
src/starcat/starcat.py
```

The planned R implementation is independently written. No code copy is required.

---

## flashier

Documentation:

https://willwerscheid.github.io/flashier/

Current `flashier-utils` dependency pin in the audited `DESCRIPTION`:

```text
fd44811ab22d61ee905e1964175571c452fe7683
```

Public capabilities relevant here:

```text
flash_init
flash_factors_init
flash_factors_fix
flash_backfit
flash_greedy
flash_fit
flash_fit_get_g
flash_ebnm
```

---

## ebnm

Documentation:

https://stephenslab.github.io/ebnm/

Relevant supported pattern:

```text
g_init + fix_g = TRUE
```

This establishes technical feasibility for fixed-prior experiments, not scientific transportability.

---

# 2. Decision log

## D001 — Keep consensus logic in `flashier-utils`

**Decision:** do not add cross-dataset catalog concepts to `flashier` core.

**Reason:** consensus programs, datasets, feature mappings and reference catalogs are application-layer concerns. `flashier` already exposes the needed inference primitives.

---

## D002 — Preserve source inferential coordinates

**Decision:** comparison normalization uses internal profiles instead of overwriting/repurposing source `factor_representation` objects.

**Reason:** comparison scale and inferential activity scale serve different purposes.

---

## D003 — All-pairs similarity, not LSAP, drives consensus graph

**Decision:** reuse existing alignment primitives before assignment but do not use `clue::solve_LSAP()` to decide consensus edges.

**Reason:** consensus clustering requires a many-dataset graph. One-to-one pairwise assignment would discard legitimate graph structure.

---

## D004 — One factor per dataset per consensus

**Decision:** prohibit duplicate source datasets inside a consensus group.

**Reason:** prevents one source dataset's factor splitting from inflating replication support and matches the intended independent-replication logic.

---

## D005 — Common reference feature basis in schema v1

**Decision:** no union-with-missing/implicit-zero consensus reference in the first schema.

**Reason:** averaging absent features as zero would conflate unmeasured with zero effect; missing-value consensus complicates projection semantics.

---

## D006 — Source offsets must be zero in schema v1

**Decision:** reject source representations with nonzero offsets for consensus construction.

**Reason:** averaging dataset-specific baselines/intercepts has no general invariant meaning. Existing ordinary source fits use zero offsets.

---

## D007 — Generic versus STARCAT-2025 modes

**Decision:** provide a generic EBMF mode and a specialized STARCAT-2025 profile contract.

**Reason:** STARCAT's variance-normalized nonnegative spectra assumptions do not hold automatically for signed or log-scale EBMF factors.

---

## D008 — No universal generic similarity threshold

**Decision:** generic mode does not silently declare 0.5 to be a biologically universal factor-similarity cutoff; `build_factor_reference(mode="generic")` requires an explicit `min_similarity`.

**Reason:** STARCAT's 0.5 threshold was calibrated in its program representation/context. Users may choose it, but generic signed factors require benchmarking.

---

## D009 — Singletons are retained diagnostically, excluded from consensus reference by default

**Decision:** default `min_datasets=2` for projection-ready consensus programs.

**Reason:** cross-dataset recurrence is the evidence defining a consensus. STARCAT manually retained a few biologically interpretable singletons; that human interpretation is not generalized automatically.

---

## D010 — Core projection stays EBMF

**Decision:** do not add a new production NNLS solver solely to mimic STARCAT.

**Reason:** the package's purpose is to extend `flashier`; existing native EBMF projection is already implemented and tested. STARCAT/NNLS remains an external validation comparator.

---

## D011 — Query prior family may be specified; learned consensus prior is not inferred automatically

**Decision:** core reference projection uses an explicit prior family with query-side hyperparameter estimation.

**Reason:** a consensus program has no unique source prior, and activity distributions are composition-dependent.

---

## D012 — Fixed-prior transport is experimental

**Decision:** `flash_fit_get_g()` and `flash_ebnm(... fix_g=TRUE)` are evaluated only in the experimental track.

**Reason:** technical availability does not establish a valid prior-combination strategy across datasets.

---

## D013 — Residual, not novel, program terminology

**Decision:** public extension is named `discover_residual_programs()`.

**Reason:** unexplained query structure may reflect batch effects, preprocessing mismatch, reference misspecification or biology. Novelty requires replication/validation.

---

## D014 — No additive known/residual PVE claim

**Decision:** compare known-only and combined reconstruction through SSE/fit-R2 changes.

**Reason:** factor contributions in nonorthogonal matrix factorizations are not generally additive.

---

## D015 — Residual discovery requires a native feasibility gate

**Decision:** `flash_greedy()`/backfit behavior with fixed reference programs is tested before the API is implemented.

**Reason:** avoids designing around assumptions about fixed-flag persistence or nullcheck behavior.

---

## D016 — Benchmarks stay out of package builds

**Decision:** downloaded real-data and heavy scientific benchmark receipts remain in the external/local validation archive.

**Reason:** this follows the repository's `AGENTS.md` instructions and keeps package checks reproducible/lightweight.

## D017 — Projection-ready schema v1 requires one compatible activity-support contract

**Decision:** a `factor_reference` may contain heterogeneous program activity supports for discovery and inspection, but schema-v1 native projection is enabled only when all selected programs share one compatible activity support.

**Reason:** the current public fixed-program projection initializes the programs jointly with one activity-side EBNM prior family. Assigning different activity priors to individual fixed programs is not part of the audited core path and must not be approximated silently.

