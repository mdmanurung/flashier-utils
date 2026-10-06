# flashier.utils 1.0.0.9000

- Add signed program networks, annotation comparisons, descriptive diagnostics,
  native ROC specificity, protein thresholds and equal-population gates.
- Add native fixed-activity and fixed-program EBMF projection with explicit
  units, preprocessing, coverage, rank and allocation checks.
- Add workflow tutorials and a published-data pancreas example adapted from
  Peter Carbonetto's single-cell Jamboree analysis. Tutorials now use real data
  (flashier's GTEx matrix; pancreas caches from `.github/scripts/prepare-pancreas.R`,
  built by the pkgdown workflow and gated off during package checks). These are
  development changes, not a release.

# flashier.utils 1.0.0

* Added signed decoupleR source/target/mor networks through pathways, with native
  ULM/MLM scoring, no default centering, explicit one-to-one feature mapping,
  network hashes and retained annotation losses. Duplicate edges and invalid
  weights are refused; weighted networks cannot use fgsea. Scores annotate
  estimated feature programs and condition on the EBMF fit.
* Implemented the 34 planned exports around public flashier and optional engines.
* Enforced named coordinates, biological units and fixed-basis backprojection.
* Added donor resampling/checkpoints, independent-fit matching and leakage-safe
  projection/masked scoring with retained failures and explicit denominators.
* Added Seurat/uwot and display-only ECDF; original inferential scores stay signed.
* Added source receipts, complete helper notices and actual runtime records.
* Explicitly narrowed unsupported NNLM/cycling, SuSiE, CorShrink dosym/bootstrap,
  frozen-noise projection and brms/uncertain-response modes. No fallback engines.
* Confirmed package identity flashier.utils and maintainer Mikhael Manurung.
  GitHub distribution with Linux release/devel, macOS and Windows checks.
* Full-engine CI verifies every suggested dependency loads; installs GSL for mashr.
