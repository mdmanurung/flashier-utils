# Upstream audit and source pins

## 1. flashier-utils baseline

```text
repository: mdmanurung/flashier-utils
commit: 75a1709e9979f66fa46beaee994443b804f9ca64
package version: 1.1.0.9000
license: MIT + file LICENSE
```

Relevant current files:

```text
R/fit.R
R/operations.R
R/representation.R
R/extract.R
R/intervals.R
R/projection.R
R/stability.R
R/metadata.R
R/plots.R
R/runtime.R
inst/provenance.csv
```

---

## 2. flashier

```text
repository: willwerscheid/flashier
version: 1.0.60
commit: fd44811ab22d61ee905e1964175571c452fe7683
license: BSD 3-Clause + file LICENSE
```

Public APIs relied on by this bundle:

```text
flash
flash_init
flash_add_intercept
flash_greedy
flash_greedy_init_default
flash_greedy_init_softImpute
flash_greedy_init_irlba
flash_backfit
flash_nullcheck
flash_factors_init
flash_factors_fix
flash_factors_unfix
flash_factors_remove
flash_factors_reorder
flash_fit
flash_fit_get_elbo
ldf
```

Important audited behaviors:

### `flash_add_intercept()`

- fixed all-ones loadings when `rowwise=TRUE`;
- fixed all-ones factor when `rowwise=FALSE`;
- initializes the opposite side from residual means;
- calls `flash_factors_init()`;
- fixes one native side;
- backfits the intercept.

### `flash_greedy()`

If `init_fn=NULL`, FLASH:

1. inspects EBNM priors;
2. infers sign constraints;
3. passes them into `flash_greedy_init_default()`.

If a custom `init_fn` is supplied, the caller owns constraint-aware
initialization.

### `flash_greedy_init_default()`

Has its own seed argument with a fixed default. Varying only an outer R seed is
not a reliable multistart mechanism.

### `flash_greedy_init_softImpute()` / `flash_greedy_init_irlba()`

Do not expose sign-constraint arguments.

### `flash_nullcheck()`

Supports `kset` and may remove zeroed factors, changing native indices.

### `ldf()`

Normalizes every native factor/loadings column and returns a separate D.
Therefore an all-ones intercept side is no longer necessarily all ones in the
normalized LDF coordinates.

---

## 3. ebnm

```text
repository: stephenslab/ebnm
version: 1.1-42
audited current commit: 825edf1e308af9f3e594a5482300b9fb99793129
license: GPL >= 3
```

Relevant current exported priors include:

```text
ebnm_point_normal
ebnm_point_laplace
ebnm_point_exponential
ebnm_normal
ebnm_normal_scale_mixture
ebnm_unimodal
ebnm_unimodal_symmetric
ebnm_unimodal_nonnegative
ebnm_unimodal_nonpositive
ebnm_generalized_binary
ebnm_npmle
ebnm_deconvolver
ebnm_flat
ebnm_horseshoe
ebnm_ash
```

`ebnm_generalized_binary()` has a default scale and can be used directly as an
exported nonnegative prior.

Parameterized wrappers/closures remain explicit-callback workflows in
`flashier-utils`.

---

## 4. GBCD package

```text
repository: stephenslab/gbcd
version: 0.2-10
commit: 670e3afb50236303ddf99db33a3394abd4a711bd
license: MIT + file LICENSE
authors: Yusha Liu, Peter Carbonetto, Jason Willwerscheid, Matthew Stephens
```

Public API:

```text
fit_gbcd
fit_gbcd_control_default
gep_annotation
```

### `fit_gbcd()` input

```text
Y: cell × gene normalized/log-transformed expression
Kmax >= 2
```

The documentation states final factor count is often near Kmax and never above:

```text
2*Kmax - 1
```

### Covariance implementation

Pinned source chooses between:

- low-rank-plus-sparse representation; or
- explicit cell × cell covariance matrix.

This motivates the wrapper memory preflight.

### Returned feature scale

`F$lfc` is on base-2 log-fold-change scale.

Documented reconstruction:

```r
tcrossprod(L, log(2) * F$lfc)
```

### Returned statistics

Pinned `fit_ebmf_to_Y()`:

- returns `F$lfc` from EB-SNMF posterior means;
- separately runs gene-wise `lm()` conditional on estimated L;
- sends coefficient estimates/SEs to `ashr::ash`;
- returns posterior-mean/posterior-SD as `z_score`;
- returns ashr local false sign rates as `lfsr`.

These are related outputs but not one single posterior uncertainty object.

### Baseline

Returned columns are named:

```text
Baseline, GEP1, GEP2, ...
```

The Baseline is not a constant all-ones sample activity offset.

---

## 5. GBCD workflow repository

```text
repository: stephenslab/gbcd-workflow
commit: e20227a1633b8b2d9d5660f88f3bfc2d7e9c001b
license file: present
```

Use only as historical/workflow context where the maintained `gbcd` package now
provides the numerical implementation.

Do not copy workflow source when package delegation is available.

---

## 6. DRIFT workflow

```text
repository: jhmarcus/drift-workflow
commit: da05b47f313f183f994155384ca6790f7b865d6e
default branch: master
license file at repository root: not observed in the audited root listing
```

Useful conceptual ideas:

```text
explicit known-factor initialization
fixed baseline factors
initialization sensitivity / multiple starts
STRUCTURE-style sample ordering
```

The code uses an older FLASH API.

Policy:

```text
conceptual attribution only
independent implementation against current public flashier APIs
no source copying unless licensing is independently clarified
```

---

## 7. Provenance classification summary

| New capability | Numerical source | Reuse classification |
|---|---|---|
| semi-NMF wrapper | flashier | call-only wrapper/adaptation |
| native intercept | flashier | call-only wrapper |
| initialized EBMF | flashier | call-only orchestration |
| multistart | flashier; DRIFT conceptual motivation | new orchestration, no DRIFT code |
| GBCD adapter | gbcd | call-only wrapper |
| group prevalence | package-owned | new QoL/descriptive |
| MDS ordering | base stats; DRIFT conceptual motivation | independent display implementation |

No upstream implementation source needs to be redistributed for these features.
