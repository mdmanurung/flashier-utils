# Milestones

## Dependency graph

```text
M0 Baseline freeze
 |
 v
M1 Foundation invariants
 |  - prior support
 |  - refit policy
 |  - projection manifest repair
 |
 +------------------+
 |                  |
 v                  v
M2 Semi-NMF     M3 Intercept -> offset
                    |
                    v
              M4 Initialized EBMF
                    |
          +---------+---------+
          |                   |
          v                   v
     M5 GBCD adapter      M6 Group prevalence
          |                   |
          +---------+---------+
                    |
                    v
              M7 Multistart
                    |
                    v
              M8 Plot ordering
```

Split-signed initialization is deferred research.

---

## M0 — Freeze current baseline

### Deliverables

- record current repository commit;
- record package version;
- save current exports;
- confirm current CI status;
- snapshot ordinary fit/projection/stability behavior;
- add any missing orientation regression fixtures.

### Exit criteria

A reproducible pre-change test baseline exists.

---

## M1 — Foundation invariants

### Scope

- prior-support correction;
- fitting-recipe/refit policy;
- projection manifest/native-index repair;
- external feature-offset metadata support.

### Deliverables

- corrected `.prior_support()`;
- `fit_recipe`;
- `refit_requires_callback`;
- stability guard;
- `_native_changed()` refit-policy update;
- projection-native metadata rebuild;
- source-vs-projection index separation;
- projection PVE correction;
- `feature_offset` metadata support in `standardize_factors()`.

### Exit criteria

- ordinary fit stability auto-refit still works;
- edited/projected fits cannot silently auto-refit as ordinary EBMF;
- projection representation PVE/native indices correspond to projection fit;
- factor uncertainty on projection representations uses correct multipliers;
- old ordinary-fit basis IDs remain unchanged.

---

## M2 — Semi-nonnegative EBMF

### Deliverables

- `fit_seminonnegative_ebmf()`;
- orientation-safe prior mapping;
- support validation;
- provenance/docs/tests.

### Exit criteria

- exact native public-call equivalence;
- nonnegative semantic activity;
- signed semantic effects;
- both matrix orientations pass.

---

## M3 — Intercept → offset

### Deliverables

- `fit_ebmf(intercept=...)`;
- intercept-only branch;
- adaptive-only nullcheck;
- native factor-role metadata;
- raw offset extraction;
- adaptive native-index manifest;
- uncertainty/interval index fixes;
- semantic PVE/rank behavior;
- editing guards;
- offset-aware new-feature projection.

### Exit criteria

- raw and LDF reconstruction identity;
- native intercept absent from semantic factor axes;
- `factor_pve(fit)` and `factor_pve(view)` agree on adaptive factors;
- uncertainty/interval extraction maps correct native columns;
- activity projection retains known offset;
- feature projection estimates new offsets;
- `max_factors=0` intercept-only fit works.

---

## M4 — Initialized / anchored EBMF

### Deliverables

- `fit_initialized_ebmf()`;
- strict basis validation;
- prior/support consistency checks;
- semantic/native mapping;
- mandatory initial backfit;
- optional factor addition;
- descending new-factor nullcheck;
- stable user IDs;
- callback-required refit policy.

### Exit criteria

- public FLASH pipeline equivalence;
- fixed side unchanged;
- initialized IDs preserved;
- nullchecked new IDs remain aligned;
- both orientations pass;
- automatic stability refit refuses.

---

## M5 — GBCD adapter

### Deliverables

- optional pinned GBCD dependency;
- `fit_gbcd_programs()`;
- no-missing/preprocessing validation;
- `kmax` semantics;
- memory preflight;
- scale conversion;
- mixed component support metadata;
- signature-statistics table;
- provenance dependency updates.

### Exit criteria

- natural-log reconstruction matches GBCD documented scale;
- Baseline retained as factor;
- global activity support is conservative;
- GEP component support recorded separately;
- conditional LM+ash statistics are not mislabeled as EB-SNMF uncertainty;
- package works without GBCD installed.

---

## M6 — Group prevalence

### Deliverables

Extend `summarize_factor_activity()` with:

```text
active_n
active_proportion
active_weighted_proportion
```

when thresholding is requested.

### Exit criteria

Manual calculations agree exactly; no inference or hidden display scaling is
introduced.

---

## M7 — Multistart optimization sensitivity

### Deliverables

- `fit_multistart_ebmf()`;
- deterministic initializer seeds;
- sign-constraint-preserving native random starts;
- safe softImpute/IRLBA restrictions;
- failure-aware candidate table;
- public ELBO extraction;
- explicit selection/refit policy.

### Exit criteria

- actual native initializer seeds differ;
- constrained priors retain required initialization signs;
- repeated configuration is reproducible;
- failed candidates are recorded;
- ELBO selection is correct;
- selected fit cannot silently enter ordinary stability refitting.

---

## M8 — STRUCTURE ordering

### Deliverables

- optional grouped one-dimensional MDS ordering;
- deterministic fallback/ties;
- display metadata.

### Exit criteria

- default order unchanged;
- groups remain contiguous only when MDS ordering requested;
- sample values/set unchanged;
- zero-distance/singleton groups do not fail.

---

# Release gates

## Gate A — Native/semantic identity

After M1/M3:

```text
native indices
representation multipliers
native PVE
offset
```

all map correctly.

## Gate B — Refitting integrity

After M4:

Every native reference fit is either:

```text
automatically reproducible
```

or:

```text
explicitly callback-required
```

with no silent fallback.

## Gate C — External-engine integrity

After M5:

GBCD outputs retain original statistical meaning and are not presented as native
FLASH posterior uncertainty.

## Gate D — Optimization integrity

After M7:

Multistart is genuinely multistart and preserves model constraints.

## Gate E — Release

- current CI matrix green;
- R devel green;
- pkgdown builds;
- provenance registry complete;
- credits reviewed;
- NEWS records projection-metadata correctness change;
- no private legacy FLASH workflow APIs introduced.
