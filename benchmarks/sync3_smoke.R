## SYNC-3 real-data smoke: refit one Stephenson cell type through the new
## flashier.utils backend and confirm parity with the committed masked fit.
.libPaths("/workspace/.Rlib")
pkgload::load_all("/workspace/repo2", quiet = TRUE)
library(flashier.utils)

f_old <- readRDS("/workspace/stephenson_K30_masked_fit.rds")
ct <- "CD4_naive"
Y <- f_old$pseudobulk[[ct]]
meta <- f_old$sample_metadata
genes <- f_old$preprocessing$feature_sets[[ct]]

# refit just this cell type through the migrated ordinary path
x <- as_cell_program_data(setNames(list(Y), ct), sample_metadata = meta)
t0 <- Sys.time()
fit_new <- canonicalize_programs(fit_celltype_programs(
  x, loading_prior = "point_laplace", var_type = 1,
  features = list(CD4_naive = genes), max_factors = 30L,
  covariates = c("institute", "sex"), covariate_mode = "residualize", seed = 1L))
cat("refit time:", format(Sys.time() - t0), "\n")

new <- fit_new$fits[[ct]]$flash
old <- f_old$fits[[ct]]$flash
cat("backend used flashier.utils:", !is.null(attr(new, "flashbridge_metadata")), "\n")
cat("old n_factors:", old$n_factors, " new n_factors:", new$n_factors, "\n")
cat("old ELBO:", old$elbo, " new ELBO:", new$elbo, "\n")
if (old$n_factors == new$n_factors) {
  cat("L_pm max abs diff:", max(abs(as.matrix(old$L_pm) - as.matrix(new$L_pm))), "\n")
  cat("F_pm max abs diff:", max(abs(as.matrix(old$F_pm) - as.matrix(new$F_pm))), "\n")
  cat("scores max abs diff:", max(abs(f_old$scores[[ct]] - fit_new$scores[[ct]])), "\n")
}
cat("fallback:", fit_new$fits[[ct]]$fallback, "\n")
cat("SMOKE DONE\n")
