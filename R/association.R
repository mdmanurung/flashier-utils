# Public engine delegation; continuous factor activities are responses.
.check_engine <- function(package, symbols = character()) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Requested engine ", package, " is unavailable; install it explicitly in your R library")
  if (any(!symbols %in% getNamespaceExports(package))) stop("Installed ", package, " lacks the required public API; inspect supported versions")
  invisible(TRUE)
}

.association_design <- function(fit, metadata, formula, engine, id_col, unit_col, independence, factors, na_action) {
  if (!inherits(formula, "formula") || length(formula) != 2L) stop("formula must be one-sided, modelling activity as the response")
  if (missing(unit_col) || is.null(unit_col) || length(unit_col) != 1L) stop("unit_col must explicitly identify the independent biological unit")
  independence <- match.arg(independence, c("independent", "repeated"))
  random <- grepl("|", paste(deparse(formula), collapse = ""), fixed = TRUE)
  if (random && engine != "dream") stop("Random effects require an explicitly supported repeated-subject engine")
  data <- .metadata_activity(fit, metadata, c(all.vars(formula), unit_col), id_col, na_action)
  .select_ids(factors, colnames(data$A), "factor")
  data$A <- data$A[, .select_ids(factors, colnames(data$A), "factor"), drop = FALSE]
  if (!ncol(data$A)) stop("No factor responses selected")
  units <- as.character(data$metadata[[unit_col]])
  if (anyNA(units) || any(!nzchar(trimws(units)))) stop("Independent unit IDs must be nonempty")
  repeats <- anyDuplicated(units) > 0L
  if (engine %in% c("lm", "limma") && (independence != "independent" || repeats)) stop("Independent engine cannot test repeated records; aggregate explicitly or use dream with a donor random effect")
  if (engine == "dream") {
    .check_engine("reformulas", c("findbars", "nobars"))
    if (independence != "repeated" || !random) stop("dream requires repeated independence and a random-effect formula")
    bars <- reformulas::findbars(formula)
    if (!any(vapply(bars, function(bar) unit_col %in% all.vars(bar[[3]]), logical(1)))) stop("dream random effects must include the declared biological unit")
    fixed <- reformulas::nobars(formula)
  } else fixed <- formula
  X <- stats::model.matrix(fixed, data$metadata)
  if (qr(X)$rank != ncol(X)) stop("Design is rank deficient; requested effects are not estimable")
  if (nrow(X) <= ncol(X)) stop("No residual degrees of freedom")
  rownames(data$metadata) <- rownames(data$A)
  data$design <- X; data$formula <- formula; data$fixed_formula <- fixed
  data$n_units <- length(unique(units)); data$independence <- independence; data$unit_col <- unit_col
  data
}

.contrast_matrix <- function(contrasts, design, include_intercept) {
  terms <- colnames(design)
  if (is.null(contrasts)) {
    index <- if (include_intercept) seq_along(terms) else which(terms != "(Intercept)")
    C <- diag(length(terms))[, index, drop = FALSE]
    dimnames(C) <- list(terms, terms[index])
    return(C)
  }
  if (is.numeric(contrasts) && is.null(dim(contrasts))) contrasts <- matrix(contrasts, ncol = 1L, dimnames = list(names(contrasts), "contrast"))
  if (!is.matrix(contrasts) || !is.numeric(contrasts) || any(!is.finite(contrasts))) stop("contrasts must be a finite named term by contrast matrix")
  .validate_ids(rownames(contrasts), "Contrast design term")
  .validate_ids(colnames(contrasts), "Contrast")
  if (any(!rownames(contrasts) %in% terms)) stop("Unknown contrast design terms")
  C <- matrix(0, length(terms), ncol(contrasts), dimnames = list(terms, colnames(contrasts)))
  C[rownames(contrasts), ] <- contrasts
  if (!ncol(C) || any(colSums(abs(C)) == 0)) stop("No nonzero contrasts to test")
  C
}

#' Test continuous factor activities using an explicit inferential engine
#'
#' Conditional on the fitted activity basis. Independent engines reject repeated
#' units. dream requires the declared unit in random effects. Native warnings and
#' errors remain visible; failures are recorded. No count preprocessing or voom.
#' @param fit Native fit or inferential view (never display percentiles).
#' @param metadata Explicit sample and biological unit metadata.
#' @param formula One-sided response design formula.
#' @param engine lm, limma or dream, chosen explicitly.
#' @param id_col Sample ID column.
#' @param unit_col Required independent biological unit column.
#' @param independence Required independent or repeated design declaration.
#' @param factors Optional unique response factors.
#' @param contrasts Named term by contrast matrix; omitted terms have zero weight.
#' @param include_intercept Include intercept tests when no contrasts supplied.
#' @param na_action fail or explicit joint omit.
#' @param adjust stats::p.adjust method.
#' @param adjust_scope all_tests or by_term.
#' @param factor_uncertainty Only ignore is supported; posterior SD is not SE.
#' @param control Engine settings: lm joint_coefficient_covariance, limma/dream ebayes,
#'   and dream native settings (dropped-method changes are forbidden).
#' @param keep_models Retain native fits for diagnostics.
#' @return table, models, design, contrasts, diagnostics, coefficient_covariance
#'   where available, and analysis_metadata identifying conditional inference.
#' @export
#' @author Mikhael Manurung; R Core Team (stats::lm); Gordon Smyth and limma contributors; Gabriel Hoffman (variancePartition::dream).
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: NEW_ORCHESTRATION / call_only.
#'   Credits stats::lm; limma; dream; later brms; use provenance("test_factors")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' metadata <- data.frame(sample=rownames(X),donor=rownames(X),
#'   group=rep(c("A","B"),each=12))
#' test_factors(view,metadata,~group,"lm",unit_col="donor",independence="independent")$table
test_factors <- function(fit, metadata, formula, engine, id_col = "sample", unit_col,
                          independence, factors = NULL, contrasts = NULL, include_intercept = FALSE,
                          na_action = "fail", adjust = "BH", adjust_scope = "all_tests",
                          factor_uncertainty = "ignore", control = list(), keep_models = TRUE) {
  if (missing(engine) || missing(unit_col) || missing(independence)) stop("Choose engine, unit_col and independence explicitly")
  engine <- match.arg(engine, c("lm", "limma", "dream"))
  if (factor_uncertainty != "ignore") stop("Response uncertainty modes are not validated; native posterior SD is not sampling SE")
  adjust_scope <- match.arg(adjust_scope, c("all_tests", "by_term"))
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control))))) stop("control must be a uniquely named list")
  data <- .association_design(fit, metadata, formula, engine, id_col, unit_col, independence, factors, na_action)
  C <- .contrast_matrix(contrasts, data$design, include_intercept)
  if (!ncol(C)) stop("No non-intercept terms to test")
  ids <- colnames(data$A)
  warnings <- character(); errors <- list(); covariance <- list()
  capture <- function(code) withCallingHandlers(force(code), warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
  table <- expand.grid(factor = ids, contrast = colnames(C), stringsAsFactors = FALSE)
  table$term <- table$contrast
  for (name in c("estimate", "std.error", "statistic", "df", "p.value")) table[[name]] <- NA_real_
  table$status <- "failed"; table$reason <- NA_character_
  models <- NULL
  if (engine == "lm") {
    if (any(!names(control) %in% "joint_coefficient_covariance")) stop("Unsupported lm control")
    models <- vector("list", length(ids)); names(models) <- ids
    response_name <- ".flashbridge_response"
    while (response_name %in% names(data$metadata)) response_name <- paste0(response_name, "_")
    model_formula <- stats::as.formula(call("~", as.name(response_name), formula[[2L]]), env = environment(formula))
    for (k in seq_along(ids)) {
      md <- data$metadata; md[[response_name]] <- data$A[, k]
      model <- tryCatch(capture(stats::lm(model_formula, data = md)), error = function(e) e)
      rows <- which(table$factor == ids[k])
      if (inherits(model, "error")) { table$reason[rows] <- conditionMessage(model); errors[[ids[k]]] <- conditionMessage(model); next }
      models[[k]] <- model
      beta <- stats::coef(model); V <- stats::vcov(model)
      estimates <- drop(crossprod(C, beta))
      variances <- colSums(C * (V %*% C))
      se <- sqrt(pmax(variances, 0)); df <- stats::df.residual(model)
      statistic <- estimates / se
      table[rows, c("estimate", "std.error", "statistic", "df", "p.value")] <- list(estimates, se, statistic, rep(df, length(rows)), 2 * stats::pt(-abs(statistic), df))
      table$status[rows] <- ifelse(is.finite(statistic), "supported", "non_estimable")
      table$reason[rows] <- ifelse(is.finite(statistic), NA_character_, "zero residual variation")
    }
    if (isTRUE(control$joint_coefficient_covariance)) {
      md <- data$metadata; md[[response_name]] <- I(data$A)
      joint <- capture(stats::lm(model_formula, data = md))
      V <- stats::vcov(joint)
      expected <- unlist(lapply(ids, function(id) paste(id, rownames(C), sep = ":")), use.names = FALSE)
      if (!identical(rownames(V), expected) || !identical(colnames(V), expected)) stop("Unprobed native mlm covariance ordering")
      for (j in seq_len(ncol(C))) {
        mapping <- matrix(0, nrow(V), length(ids), dimnames = list(rownames(V), ids))
        for (k in seq_along(ids)) mapping[((k - 1L) * nrow(C) + 1L):(k * nrow(C)), k] <- C[, j]
        covariance[[colnames(C)[j]]] <- crossprod(mapping, V %*% mapping)
      }
      models$joint <- joint
    }
  } else {
    .check_engine("limma", c("lmFit", "contrasts.fit", "eBayes"))
    allowed <- if (engine == "dream") c("ebayes", "dream") else "ebayes"
    if (any(!names(control) %in% allowed)) stop("Unsupported engine control")
    ebayes <- control$ebayes
    if (is.null(ebayes)) ebayes <- list()
    if (!is.list(ebayes) || (length(ebayes) && is.null(names(ebayes))) || any(!names(ebayes) %in% c("proportion", "stdev.coef.lim", "trend", "robust", "winsor.tail.p"))) stop("Unsupported moderation setting")
    native_contrast_names <- colnames(C)
    if (engine == "limma") {
      model <- capture(limma::lmFit(t(data$A), data$design))
      model <- limma::contrasts.fit(model, C)
      model <- capture(do.call(limma::eBayes, c(list(fit = model), ebayes)))
    } else {
      .check_engine("variancePartition", c("dream", "eBayes"))
      settings <- control$dream
      if (is.null(settings)) settings <- list()
      if (!is.list(settings) || (length(settings) && is.null(names(settings))) || any(!names(settings) %in% c("ddf", "control", "BPPARAM"))) stop("Unsupported dream setting")
      native_C <- C
      native_contrast_names <- utils::tail(make.unique(c(colnames(data$design), colnames(C))), ncol(C))
      colnames(native_C) <- native_contrast_names
      trivial <- apply(C, 2L, function(column) sum(column != 0) == 1L && sum(column) == 1)
      if (all(trivial)) {
        native_contrast_names <- rownames(C)[apply(C, 2L, which.max)]
        native_L <- list()
      } else native_L <- list(L = native_C)
      model <- capture(do.call(variancePartition::dream, c(list(exprObj = t(data$A), formula = formula, data = data$metadata), native_L, settings)))
      errors <- attr(model, "errors")
      model <- capture(do.call(variancePartition::eBayes, c(list(fit = model), ebayes)))
    }
    for (i in seq_len(nrow(table))) {
      k <- match(table$factor[i], rownames(model$coefficients)); j <- match(native_contrast_names[match(table$contrast[i], colnames(C))], colnames(model$coefficients))
      if (is.na(k) || is.na(j)) { table$reason[i] <- "native fit omitted failed factor/contrast"; next }
      se <- model$stdev.unscaled[k, j] * sqrt(model$s2.post[k])
      df <- if (is.matrix(model$df.total)) model$df.total[k, j] else model$df.total[k]
      table[i, c("estimate", "std.error", "statistic", "df", "p.value")] <- list(model$coefficients[k, j], se, model$t[k, j], df, model$p.value[k, j])
      table$status[i] <- if (is.finite(table$p.value[i])) "supported" else "non_estimable"
    }
    models <- model
  }
  table$adj.p.value <- if (adjust_scope == "all_tests") stats::p.adjust(table$p.value, method = adjust) else stats::ave(table$p.value, table$term, FUN = function(p) stats::p.adjust(p, method = adjust))
  table$engine <- engine; table$n_samples <- nrow(data$A); table$n_units <- data$n_units
  table$basis_id <- data$view$manifest$basis_id; table$inference <- "conditional_on_factor_fit"
  metadata <- c(data$view$manifest, list(engine = engine, engine_version = if (engine == "lm") as.character(getRversion()) else as.character(utils::packageVersion(if (engine == "dream") "variancePartition" else "limma")),
    formula = formula, unit_col = unit_col, independence = independence, matching = data$matching, exclusions = data$exclusions,
    adjustment = adjust, adjust_scope = adjust_scope, tested_family = table[, c("factor", "contrast")],
    moderation = control$ebayes, small_factor_moderation_caveat = engine != "lm", inference = "conditional_on_factor_fit"))
  attr(table, "analysis_metadata") <- data$view$manifest
  .attach_provenance(list(table = table, models = if (keep_models) models else NULL, design = data$design, contrasts = C,
       coefficient_covariance = covariance, diagnostics = list(warnings = warnings, errors = errors, native_attributes = if (engine == "dream") attributes(models) else NULL), analysis_metadata = metadata), "test_factors", match.call(), .resolved_parameters("test_factors", environment()))
}
