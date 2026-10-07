restore_programs <- function(ldf,signed_effects=FALSE) {
  Y <- cache$reference_expression
  factor_ids <- paste0("F",seq_len(ncol(ldf$F)))
  A <- sweep(ldf$L[rownames(Y),,drop=FALSE],2,ldf$D,"*")
  B <- ldf$F
  colnames(A) <- colnames(B) <- factor_ids
  priors <- list(ebnm::ebnm_point_exponential,
                 if(signed_effects) ebnm::ebnm_point_normal
                 else ebnm::ebnm_point_exponential)
  object <- fit_ebmf(Y,"rows",ebnm_fn=priors,max_factors=0,var_type=2,
                     feature_scale="shifted_log_reference_library_normalized_measurements",
                     verbose=0)
  metadata <- attr(object,"flashbridge_metadata")
  metadata$factor_ids <- factor_ids
  metadata$preprocessing <- list(scope="external_reference",revision=cache$revision,
    transformation="log1p(filtered measurement / library size * reference mean)",
    reference_library_mean=cache$reference_library_mean,
    source_sha256=cache$source_sha256)
  object <- flashier::flash_factors_init(object,list(A,B),ebnm_fn=priors)
  object <- flashier::flash_factors_fix(object,seq_along(factor_ids),"factors")
  # Native calls replace the object; preserve the wrapper declarations.
  attr(object,"flashbridge_metadata") <- metadata
  view <- standardize_factors(object,sample_side="rows",scaling="raw",
    prior_support=c(activity="nonnegative",
                    effects=if(signed_effects) "signed" else "nonnegative"))
  stopifnot(identical(unname(view$activity),unname(A)),
            identical(unname(view$effects),unname(B)))
  view
}
