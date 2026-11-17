# Public engine delegation; continuous factor activities are responses.
.check_engine <- function(package, symbols = character()) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Requested engine ", package, " is unavailable; install it explicitly in your R library")
  }
  if (any(!symbols %in% getNamespaceExports(package))) {
    stop("Installed ", package, " lacks the required public API; inspect supported versions")
  }
  invisible(TRUE)
}

.association_design <- function(fit, metadata, formula, engine, id_col, unit_col, independence, factors, na_action) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("formula must be one-sided, modelling activity as the response")
  }
  if (missing(unit_col) || is.null(unit_col) || length(unit_col) != 1L) {
    stop("unit_col must explicitly identify the independent biological unit")
  }
  independence <- match.arg(independence, c("independent", "repeated"))
  random <- grepl("|", paste(deparse(formula), collapse = ""), fixed = TRUE)
  if (random && engine != "dream") {
    stop("Random effects require an explicitly supported repeated-subject engine")
  }
  data <- .metadata_activity(fit, metadata, c(all.vars(formula), unit_col), id_col, na_action)
  .require_biological_factors(colnames(data$A))
  .select_ids(factors, colnames(data$A), "factor")
  data$A <- data$A[, .select_ids(factors, colnames(data$A), "factor"), drop = FALSE]
  if (!ncol(data$A)) {
    stop("No factor responses selected")
  }
  if (identical(unit_col, id_col)) {
    message(
      "unit_col is the sample ID, so every sample is treated as its own biological unit. ",
      "Name the donor or subject column if samples share a unit."
    )
  }
  units <- as.character(data$metadata[[unit_col]])
  if (anyNA(units) || any(!nzchar(trimws(units)))) {
    stop("Independent unit IDs must be nonempty")
  }
  repeats <- anyDuplicated(units) > 0L
  if (engine %in% c("lm", "limma") && (independence != "independent" || repeats)) {
    stop(
      "Independent engine cannot test repeated records; aggregate explicitly or use dream with a donor random effect"
    )
  }
  if (engine == "dream") {
    .check_engine("reformulas", c("findbars", "nobars"))
    if (independence != "repeated" || !random) {
      stop("dream requires repeated independence and a random-effect formula")
    }
    bars <- reformulas::findbars(formula)
    # the bar must group by the biological unit itself: a finer grouping (donor:visit) gives each
    # repeat its own intercept, so the donor correlation is not modelled and tests are anti-conservative
    has_unit_intercept <- function(bar) {
      if (attr(stats::terms(stats::as.formula(call("~", bar[[2]]))), "intercept") != 1L) {
        return(FALSE)
      }
      vars <- all.vars(bar[[3]])
      if (!length(vars) || !all(vars %in% names(data$metadata))) {
        return(FALSE)
      }
      group <- interaction(data$metadata[vars], drop = TRUE)
      n_pairs <- nrow(unique(data.frame(group = group, unit = units)))
      # one-to-one: each level lies within one unit and each unit within one level
      n_pairs == nlevels(group) && n_pairs == length(unique(units))
    }
    if (!any(vapply(bars, has_unit_intercept, logical(1)))) {
      stop(
        "dream random effects must include a random intercept for the declared biological unit; ",
        "the grouping must be the biological unit (one level per unit_col value), not a finer one such as donor:visit"
      )
    }
    fixed <- reformulas::nobars(formula)
  } else {
    fixed <- formula
  }
  if (engine != "lm" && length(attr(stats::terms(fixed), "offset"))) {
    stop("Formula offsets are supported only by engine = 'lm'")
  }
  X <- stats::model.matrix(fixed, data$metadata)
  if (qr(X)$rank != ncol(X)) {
    stop("Design is rank deficient; requested effects are not estimable")
  }
  if (nrow(X) <= ncol(X)) {
    stop("No residual degrees of freedom")
  }
  rownames(data$metadata) <- rownames(data$A)
  data$design <- X
  data$formula <- formula
  data$fixed_formula <- fixed
  data$n_units <- length(unique(units))
  data$independence <- independence
  data$unit_col <- unit_col
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
  if (is.numeric(contrasts) && is.null(dim(contrasts))) {
    contrasts <- matrix(contrasts, ncol = 1L, dimnames = list(names(contrasts), "contrast"))
  }
  if (!is.matrix(contrasts) || !is.numeric(contrasts) || any(!is.finite(contrasts))) {
    stop("contrasts must be a finite named term by contrast matrix")
  }
  .validate_ids(rownames(contrasts), "Contrast design term")
  .validate_ids(colnames(contrasts), "Contrast")
  if (any(!rownames(contrasts) %in% terms)) {
    stop("Unknown contrast design terms")
  }
  C <- matrix(0, length(terms), ncol(contrasts), dimnames = list(terms, colnames(contrasts)))
  C[rownames(contrasts), ] <- contrasts
  if (!ncol(C) || any(colSums(abs(C)) == 0)) {
    stop("No nonzero contrasts to test")
  }
  C
}

#' Test continuous factor activities using an explicit inferential engine
#'
#' Regress each factor's activity on a sample-level design and return one row
#' per factor and term, with p-values. Conditional on the fitted activity basis. Independent engines reject repeated
#' units. dream requires the declared unit in random effects. Native warnings and
#' errors remain visible; failures are recorded. No count preprocessing or voom.
#' @param fit Native fit or a view from [standardize_factors()] (never display percentiles).
#' @param metadata Data frame with one row per sample, including the design
#'   columns, the sample ID and the biological unit.
#' @param formula One-sided design formula over `metadata` columns, e.g. `~ group + age`.
#' @param engine `"lm"` (stats), `"limma"` or `"dream"` (variancePartition),
#'   chosen explicitly; the last two need their packages installed. limma always
#'   applies empirical Bayes moderation; with few factors it reproduces lm (identical
#'   SE, df and p at K = 2) and needs K in the tens to borrow strength. Moderation
#'   assumes the factors' residual variances are on a similar scale, so use `lm` for
#'   views in raw scale (a warning appears when they differ by more than tenfold).
#'   dream reports its own Satterthwaite t, df and p unless `moderate = TRUE`.
#' @param id_col Name of the sample ID column (default `"sample"`).
#' @param unit_col Name of the biological unit column, e.g. `"donor"`. Required.
#' @param independence `"independent"` (one sample per unit) or `"repeated"`
#'   (several samples per unit; needs `engine = "dream"` with the unit as a
#'   random effect). Required.
#' @param factors Optional character factor IDs to test. Default all.
#' @param contrasts Optional numeric matrix, rows named by design terms and
#'   columns by contrast; omitted terms have zero weight.
#' @param include_intercept Also test the intercept when no contrasts are
#'   supplied (default `FALSE`).
#' @param na_action `"fail"` (default) or `"omit"` rows with NA jointly.
#' @param adjust Multiple-testing method name for [stats::p.adjust()] (default `"BH"`).
#' @param adjust_scope `"all_tests"` (default) or `"by_term"`.
#' @param control Optional named list of engine settings: lm `joint_coefficient_covariance`,
#'   limma/dream `ebayes`, and dream native settings (dropped-method changes are forbidden).
#' @param keep_models Retain native model fits for diagnostics (default `FALSE`).
#' @param moderate `engine = "dream"` only: apply empirical Bayes moderation
#'   ([variancePartition::eBayes()]) after the fit (default `FALSE`). The moderated
#'   df pools residual df across factors (about 8000 in a simulation where the
#'   between-donor df is 8), which makes p-values for terms that vary between donors
#'   anti-conservative once there are three or more factors. limma always
#'   moderates and `lm` never does, so the argument is ignored for them.
#' @return List with `$table` (one row per factor and term: `estimate`,
#'   `std.error`, `statistic`, `df`, `p.value`, `adj.p.value`), `$models`, `$design`, `$contrasts`, `$diagnostics`,
#'   `$coefficient_covariance` where available, and `$analysis_metadata`
#'   identifying conditional inference. With `adjust_scope = "all_tests"` the
#'   Benjamini-Hochberg family is the non-intercept tests; intercept tests (from
#'   `include_intercept = TRUE`) are adjusted within their own family.
#'   `adj.p.value` is adjusted only within this call.
#' @seealso [backproject_contrast()] to express a term in feature units
#' @export
#' @author Mikhael Manurung; R Core Team (stats::lm); Gordon Smyth and limma contributors; Gabriel Hoffman (variancePartition::dream).
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: NEW_ORCHESTRATION / call_only.
#'   Credits stats::lm; limma; dream; later brms; use factor_provenance("test_factor_association")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48), 24, 2), matrix(rnorm(32), 16, 2)) +
#'   matrix(rnorm(384, sd = 0.2), 24, 16)
#' dimnames(X) <- list(paste0("s", 1:24), paste0("g", 1:16))
#' fit <- fit_ebmf(X, "rows",
#'   max_factors = 2, seed = 4, verbose = 0,
#'   ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
#'   nullcheck = FALSE, backfit = TRUE, feature_scale = "simulated_centered_intensity"
#' )
#' view <- standardize_factors(fit)
#' # one sample per donor, so the sample ID can serve as the unit; name the donor column otherwise
#' metadata <- data.frame(
#'   sample = rownames(X), donor = rownames(X),
#'   group = rep(c("A", "B"), each = 12)
#' )
#' test_factor_association(view, metadata = metadata, formula = ~group, engine = "lm", unit_col = "donor", independence = "independent")$table
test_factor_association <- function(
  fit,
  metadata,
  formula,
  engine,
  unit_col,
  independence,
  id_col = "sample",
  factors = NULL,
  contrasts = NULL,
  include_intercept = FALSE,
  na_action = "fail",
  adjust = "BH",
  adjust_scope = "all_tests",
  control = list(),
  keep_models = TRUE,
  moderate = FALSE
) {
  if (!.is_flag(moderate)) {
    stop("moderate must be TRUE or FALSE")
  }
  if (missing(engine) || missing(unit_col) || missing(independence)) {
    stop("Choose engine, unit_col and independence explicitly")
  }
  engine <- match.arg(engine, c("lm", "limma", "dream"))
  adjust_scope <- match.arg(adjust_scope, c("all_tests", "by_term"))
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  if (!.is_named_list(control)) {
    stop("control must be a uniquely named list")
  }
  data <- .association_design(fit, metadata, formula, engine, id_col, unit_col, independence, factors, na_action)
  C <- .contrast_matrix(contrasts, data$design, include_intercept)
  if (!ncol(C)) {
    stop("No non-intercept terms to test")
  }
  ids <- colnames(data$A)
  noted <- new.env()
  noted$warnings <- character()
  errors <- list()
  covariance <- list()
  capture <- function(code) {
    withCallingHandlers(force(code), warning = function(w) noted$warnings <- c(noted$warnings, conditionMessage(w)))
  }
  table <- expand.grid(factor = ids, contrast = colnames(C), stringsAsFactors = FALSE)
  table$term <- table$contrast
  for (name in c("estimate", "std.error", "statistic", "df", "p.value")) {
    table[[name]] <- NA_real_
  }
  table$status <- "failed"
  table$reason <- NA_character_
  models <- NULL
  if (engine == "lm") {
    if (any(!names(control) %in% "joint_coefficient_covariance")) {
      stop("Unsupported lm control")
    }
    models <- vector("list", length(ids))
    names(models) <- ids
    response_name <- ".flashbridge_response"
    while (response_name %in% names(data$metadata)) {
      response_name <- paste0(response_name, "_")
    }
    model_formula <- stats::as.formula(call("~", as.name(response_name), formula[[2L]]), env = environment(formula))
    for (k in seq_along(ids)) {
      md <- data$metadata
      md[[response_name]] <- data$A[, k]
      model <- tryCatch(capture(stats::lm(model_formula, data = md)), error = function(e) e)
      rows <- which(table$factor == ids[k])
      if (inherits(model, "error")) {
        table$reason[rows] <- conditionMessage(model)
        errors[[ids[k]]] <- conditionMessage(model)
        next
      }
      models[[k]] <- model
      beta <- stats::coef(model)
      V <- stats::vcov(model)
      estimates <- drop(crossprod(C, beta))
      variances <- colSums(C * (V %*% C))
      se <- sqrt(pmax(variances, 0))
      df <- stats::df.residual(model)
      statistic <- estimates / se
      table[rows, c("estimate", "std.error", "statistic", "df", "p.value")] <- list(
        estimates,
        se,
        statistic,
        rep(df, length(rows)),
        2 * stats::pt(-abs(statistic), df)
      )
      table$status[rows] <- ifelse(is.finite(statistic), "supported", "non_estimable")
      table$reason[rows] <- ifelse(is.finite(statistic), NA_character_, "zero residual variation")
    }
    if (isTRUE(control$joint_coefficient_covariance)) {
      md <- data$metadata
      md[[response_name]] <- I(data$A)
      joint <- capture(stats::lm(model_formula, data = md))
      V <- stats::vcov(joint)
      expected <- unlist(lapply(ids, function(id) paste(id, rownames(C), sep = ":")), use.names = FALSE)
      if (!identical(rownames(V), expected) || !identical(colnames(V), expected)) {
        stop("Unprobed native mlm covariance ordering")
      }
      for (j in seq_len(ncol(C))) {
        mapping <- matrix(0, nrow(V), length(ids), dimnames = list(rownames(V), ids))
        for (k in seq_along(ids)) {
          mapping[((k - 1L) * nrow(C) + 1L):(k * nrow(C)), k] <- C[, j]
        }
        covariance[[colnames(C)[j]]] <- crossprod(mapping, V %*% mapping)
      }
      models$joint <- joint
    }
  } else {
    .check_engine("limma", c("lmFit", "contrasts.fit", "eBayes"))
    allowed <- if (engine == "dream") c("ebayes", "dream") else "ebayes"
    if (any(!names(control) %in% allowed)) {
      stop("Unsupported engine control")
    }
    ebayes <- control$ebayes
    if (is.null(ebayes)) {
      ebayes <- list()
    }
    if (engine == "dream" && !moderate && length(ebayes)) {
      stop("control$ebayes needs moderate = TRUE for engine = 'dream'")
    }
    if (
      !is.list(ebayes) ||
        (length(ebayes) && is.null(names(ebayes))) ||
        any(!names(ebayes) %in% c("proportion", "stdev.coef.lim", "trend", "robust", "winsor.tail.p"))
    ) {
      stop("Unsupported moderation setting")
    }
    native_contrast_names <- colnames(C)
    if (engine == "limma") {
      model <- capture(limma::lmFit(t(data$A), data$design))
      model <- limma::contrasts.fit(model, C)
      model <- capture(do.call(limma::eBayes, c(list(fit = model), ebayes)))
    } else {
      .check_engine("variancePartition", c("dream", "eBayes"))
      settings <- control$dream
      if (is.null(settings)) {
        settings <- list()
      }
      if (
        !is.list(settings) ||
          (length(settings) && is.null(names(settings))) ||
          any(!names(settings) %in% c("ddf", "control", "BPPARAM"))
      ) {
        stop("Unsupported dream setting")
      }
      native_C <- C
      native_contrast_names <- utils::tail(make.unique(c(colnames(data$design), colnames(C))), ncol(C))
      colnames(native_C) <- native_contrast_names
      trivial <- apply(C, 2L, function(column) sum(column != 0) == 1L && sum(column) == 1)
      if (all(trivial)) {
        native_contrast_names <- rownames(C)[apply(C, 2L, which.max)]
        native_L <- list()
      } else {
        native_L <- list(L = native_C)
      }
      model <- capture(do.call(
        variancePartition::dream,
        c(list(exprObj = t(data$A), formula = formula, data = data$metadata), native_L, settings)
      ))
      errors <- attr(model, "errors")
      if (moderate) {
        model <- capture(do.call(variancePartition::eBayes, c(list(fit = model), ebayes)))
      }
    }
    # a single-factor fit comes back without row names; the rows are in factor order
    if (is.null(rownames(model$coefficients)) && nrow(model$coefficients) == length(ids)) {
      rownames(model$coefficients) <- ids
    }
    if ((engine == "limma" || moderate) && length(ids) > 1L) {
      residual_variance <- model$sigma^2
      positive <- residual_variance[is.finite(residual_variance) & residual_variance > 0]
      # a zero-variance factor must not silence the check: say so, and compare the rest
      has_zero <- any(is.finite(residual_variance) & residual_variance <= 0)
      too_wide <- length(positive) > 1L && max(positive) / min(positive) > 10
      if (has_zero || too_wide || length(positive) < 2L) {
        msg <- paste0(
          if (has_zero) {
            "At least one factor has zero residual variance, so the variance ratio is undefined for it. "
          } else {
            ""
          },
          if (length(positive) < 2L) "Fewer than two factors have positive residual variance. " else "",
          if (too_wide) {
            "Factor residual variances differ more than tenfold; moderation assumes a similar scale. "
          } else {
            ""
          },
          "Consider engine = 'lm' or standardised factors."
        )
        noted$warnings <- c(noted$warnings, msg)
        warning(msg, call. = FALSE)
      }
    }
    for (i in seq_len(nrow(table))) {
      k <- match(table$factor[i], rownames(model$coefficients))
      j <- match(native_contrast_names[match(table$contrast[i], colnames(C))], colnames(model$coefficients))
      if (is.na(k) || is.na(j)) {
        table$reason[i] <- "native fit omitted failed factor/contrast"
        next
      }
      # moderated fits carry s2.post and df.total; an unmoderated dream fit has sigma and per-term df
      moderated <- !is.null(model$s2.post)
      se <- model$stdev.unscaled[k, j] * sqrt(if (moderated) model$s2.post[k] else model$sigma[k]^2)
      df <- if (!moderated) {
        model$df.residual[k, j]
      } else if (is.matrix(model$df.total)) {
        model$df.total[k, j]
      } else {
        model$df.total[k]
      }
      table[i, c("estimate", "std.error", "statistic", "df", "p.value")] <- list(
        model$coefficients[k, j],
        se,
        model$t[k, j],
        df,
        model$p.value[k, j]
      )
      table$status[i] <- if (is.finite(table$p.value[i])) "supported" else "non_estimable"
    }
    models <- model
  }
  # a zero standard error gives p = 0; a non-estimable test has no p-value, so it stays out of the BH family
  table$p.value[table$status == "non_estimable"] <- NA_real_
  is_intercept <- table$term == "(Intercept)"
  table$adj.p.value <- if (adjust_scope == "all_tests") {
    # intercept tests (include_intercept = TRUE) form their own family, so they do not change the other adjusted p-values
    stats::ave(table$p.value, is_intercept, FUN = function(p) stats::p.adjust(p, method = adjust))
  } else {
    stats::ave(table$p.value, table$term, FUN = function(p) stats::p.adjust(p, method = adjust))
  }
  table$engine <- engine
  table$n_samples <- nrow(data$A)
  table$n_units <- data$n_units
  table$basis_id <- data$view$manifest$basis_id
  table$inference <- "conditional_on_factor_fit"
  metadata <- c(
    data$view$manifest,
    list(
      engine = engine,
      engine_version = if (engine == "lm") {
        as.character(getRversion())
      } else {
        as.character(utils::packageVersion(if (engine == "dream") "variancePartition" else "limma"))
      },
      formula = formula,
      unit_col = unit_col,
      independence = independence,
      matching = data$matching,
      exclusions = data$exclusions,
      adjustment = adjust,
      adjust_scope = adjust_scope,
      tested_family = table[, c("factor", "contrast")],
      moderation = if (engine == "lm") {
        "none"
      } else if (engine == "limma" || moderate) {
        "ebayes"
      } else {
        "none"
      },
      moderation_settings = control$ebayes,
      small_factor_moderation_caveat = engine == "limma" || (engine == "dream" && moderate),
      inference = "conditional_on_factor_fit"
    )
  )
  attr(table, "analysis_metadata") <- data$view$manifest
  .attach_provenance(
    list(
      table = table,
      models = if (keep_models) models else NULL,
      design = data$design,
      contrasts = C,
      coefficient_covariance = covariance,
      diagnostics = list(
        warnings = noted$warnings,
        errors = errors,
        native_attributes = if (engine == "dream") attributes(models) else NULL
      ),
      analysis_metadata = metadata
    ),
    "test_factor_association",
    match.call(),
    .resolved_parameters("test_factor_association", environment())
  )
}
