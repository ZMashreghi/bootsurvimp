#' Pseudo-Population Bootstrap for Imputed Survey Data
#'
#' @description
#' Pseudo-population bootstrap variance estimation for imputed survey data
#' under stratified sampling (simple random sampling without replacement or
#' conditional Poisson sampling within strata). Missing values of the study
#' variable are imputed by doubly robust (DR) or multiply robust (MR)
#' imputation, either deterministic regression imputation or random hot-deck
#' imputation within imputation classes. Bootstrap pseudo-populations are
#' generated under the nonresponse-model (`"NRM"`) or imputation-model
#' (`"IM"`) approach, or both combined in the multiply robust (`"MR"`) scheme.
#' A pseudo-population version of the Shao and Sitter (1996) reimputation
#' bootstrap (`"Shao.Sitter"`) is also available with mean, ratio, hot-deck or
#' regression imputation.
#'
#' @param data A matrix or data frame containing the sample. The study-variable
#' column is identified through the `study.variable` argument (by default, the
#' column named `study.variable`) and, under stratified sampling, the stratum
#' column through the `stratum` argument. The study variable may contain
#' missing values (`NA`) for the nonrespondents, to be imputed within the
#' function, or it may be already completed, with a response-indicator column
#' (see `resp.ind`) identifying the respondents. A sample of respondents only
#' can be used with the `"NRM"` scheme if a response-probability column (see
#' `resp.prob`) is supplied. Optional columns: inclusion probabilities (see
#' `incl.prob`) and `ksi`, the size variable used to compute the inclusion
#' probabilities in the bootstrap pseudo-populations (1 for every unit if
#' absent).
#'
#' @param population.size A vector of strata population sizes. If it is a named
#' vector, its names must match the stratum labels; if it is unnamed, its
#' elements must follow the sorted order of the stratum labels (as in the
#' `bootsurv` package). For an unstratified design, a single population size.
#'
#' @param R.pop The number of bootstrap replicates to create bootstrap
#' pseudo-populations. The default is `10`.
#'
#' @param R.samp The number of bootstrap replicates to draw bootstrap samples
#' from each bootstrap pseudo-population. The default is `100`.
#'
#' @param study.variable A character string giving the name of the column in
#' `data` that holds the study-variable values. The default is
#' `"study.variable"`.
#'
#' @param stratum A character string giving the name of the column in `data`
#' that identifies the strata. The default is `NULL`, in which case a single
#' stratum is assumed. For a stratified design, set `stratum` to the name of
#' the column in `data` that identifies the strata.
#'
#' @param resp.ind A character string giving the name of the column in `data`
#' that holds the response indicators (1 = respondent, 0 = nonrespondent). The
#' default is `"resp.ind"`. If the column is absent, respondents are identified
#' by the non-missing values of the study variable.
#'
#' @param resp.prob A character string giving the name of the column in `data`
#' that holds known response probabilities. The default is `"resp.prob"`. If
#' the column is present, the `"NRM"` scheme uses these probabilities instead
#' of estimating them; it is required for a sample of respondents only.
#'
#' @param incl.prob A character string giving the name of the column in `data`
#' that holds the inclusion probabilities, which must lie in (0, 1]. The default
#' is `"pi"`. If the column is absent, the inclusion probabilities are computed
#' as `nh/Nh` within each stratum.
#'
#' @param sample.size A vector of the full-sample sizes within strata (named
#' by the stratum labels, or in their sorted order), needed only when `data`
#' contains the respondents only (`"NRM"` scheme), since the full-sample sizes
#' cannot be obtained from the respondents. The bootstrap samples are drawn
#' with these sizes, and the default inclusion probabilities are
#' `sample.size/population.size`. The default is `NULL`.
#'
#' @param aux.IM A list of covariate matrices (one row per unit of `data`), one
#' matrix for each candidate imputation model. An intercept is added
#' automatically. The doubly robust methods and the `"Shao.Sitter"` scheme use
#' the first matrix of the list.
#'
#' @param aux.NRM A list of covariate matrices (one row per unit of `data`), one
#' matrix for each candidate nonresponse model. An intercept is added
#' automatically. The doubly robust methods use the first matrix of the list.
#'
#' @param parameter One of the following population parameters can be
#' specified: `"total"` (population total), `"mean"` (population mean),
#' `"quartile.25"` (population 1st quartile), `"quartile.50"` or `"median"`
#' (population median) or `"quartile.75"` (population 3rd quartile). If the
#' parameter of interest is the population mean or total, the HT-estimator is
#' applied. If the parameter of interest is a population quartile, the
#' estimator in Särndal, Swensson, and Wretman (1992, Chapter 5) is applied.
#' The default is the population total.
#'
#' @param imputation.method The imputation method. For the `"NRM"`, `"IM"` and
#' `"MR"` schemes, one of `"MRDRegI"` (multiply robust deterministic
#' regression imputation), `"MRCLRHDI"` (multiply robust random hot-deck
#' imputation within classes), `"DRDRegI"` (doubly robust deterministic
#' regression imputation) or `"DRCLRHDI"` (doubly robust random hot-deck
#' imputation within classes). For the `"Shao.Sitter"` scheme, one of `"mean"`,
#' `"ratio"`, `"hotdeck"` or `"regression"`.
#'
#' @param bootstrap.method One of the following bootstrap schemes can be
#' applied:
#' \itemize{
#'   \item `"NRM"`: nonresponse-model scheme.
#'   \item `"IM"`: imputation-model scheme.
#'   \item `"MR"`: multiply robust scheme, combining the NRM and IM schemes.
#'   \item `"Shao.Sitter"`: pseudo-population version of the Shao and Sitter (1996) scheme.
#' }
#'
#' @param survey.design The sampling design within strata: `"SRSWOR"` (simple
#' random sampling without replacement) or `"CPS"` (conditional Poisson
#' sampling). The default is `"SRSWOR"`.
#'
#' @param n.classes The number of imputation classes for random hot-deck
#' imputation within classes (`"MRCLRHDI"` and `"DRCLRHDI"`). The classes
#' have (approximately) equal sizes.
#'
#' @param R.MR.samp The number of inner bootstrap samples used by the
#' `"MR"` scheme. The default is `10`.
#'
#' @param R.MR.pop The number of inner bootstrap pseudo-populations used by
#' the `"MR"` scheme. The default is `5`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`boot.statistic`}{A vector of bootstrap statistics (`NULL` for the `"MR"` scheme).}
#'   \item{`boot.var`}{The bootstrap variance estimator of the estimator of the parameter of interest.}
#'   \item{`boot.mean`}{The average of the bootstrap estimator of the parameter of interest (`NA` for the `"MR"` scheme).}
#'   \item{`sample.size`}{The number of sampled units in each stratum.}
#'   \item{`parameter`}{The parameter of interest.}
#'   \item{`imputation.method`}{The imputation method.}
#'   \item{`bootstrap.method`}{The bootstrap scheme.}
#'   \item{`number.bootstrap.pop.replicates`}{The number of bootstrap pseudo-populations, `R.pop`.}
#'   \item{`number.bootstrap.samp.replicates`}{The number of bootstrap samples drawn from each pseudo-population, `R.samp`.}
#' }
#'
#' @references
#' Chen, S., Haziza, D., Léger, C. and Mashreghi, Z. (2019). Pseudo-population
#' bootstrap methods for imputed survey data. Biometrika, 106(2), 369-384.
#' \doi{10.1093/biomet/asz001}
#'
#' Chen, S., Haziza, D. and Mashreghi, Z. (2021). Multiply robust bootstrap
#' variance estimation in the presence of singly imputed survey data. Journal
#' of Survey Statistics and Methodology, 9(4), 810-832.
#' \doi{10.1093/jssam/smaa004}
#'
#' Shao, J. and Sitter, R.R. (1996). Bootstrap for Imputed Survey Data.
#' Journal of the American Statistical Association, 91(435), 1278-1288.
#'
#' @seealso [rescale.boot.imp()] for rescaling bootstrap methods under mean,
#'   ratio, regression and random hot-deck imputation.
#'
#' @examples
#' ## Stratified SRS without replacement (overall 20% sampling fraction;
#' ## strata = age groups), uniform nonresponse within strata.
#' ## Covariates are passed as LISTS of matrices: one matrix per candidate
#' ## imputation model (aux.IM) and per candidate nonresponse model (aux.NRM).
#' aux.im  <- list(as.matrix(data_samp_stsrs_unif[, c("x1", "x2", "x3")]))
#' aux.nrm <- list(as.matrix(data_samp_stsrs_unif[, c("x1", "x2")]))
#' ## Small R.pop / R.samp are used here for speed; use more in practice.
#'
#' ## --- Nonresponse-model (NRM) scheme, doubly robust regression imputation ---
#' set.seed(1)
#' fit.nrm <- pseudopop.boot.imp(data_samp_stsrs_unif, pop_size_st_imp,
#'                               R.pop = 2, R.samp = 20, stratum = "stratum",
#'                               aux.IM = aux.im, aux.NRM = aux.nrm,
#'                               parameter = "mean",
#'                               imputation.method = "DRDRegI",
#'                               bootstrap.method = "NRM")
#' fit.nrm$boot.var               # bootstrap variance estimator
#' head(fit.nrm$boot.statistic)   # bootstrap replicate estimates
#'
#' ## How close is the bootstrap mean to the true population mean?
#' c(estimate = fit.nrm$boot.mean, true = mean(data_pop_st_imp$study.variable))
#'
#' ## --- Imputation-model (IM) scheme, doubly robust class hot-deck, median ---
#' set.seed(1)
#' fit.im <- pseudopop.boot.imp(data_samp_stsrs_unif, pop_size_st_imp,
#'                              R.pop = 2, R.samp = 20, stratum = "stratum",
#'                              aux.IM = aux.im, aux.NRM = aux.nrm,
#'                              parameter = "median",
#'                              imputation.method = "DRCLRHDI",
#'                              bootstrap.method = "IM", n.classes = 4)
#' fit.im$boot.var
#'
#' ## --- Multiply robust imputation with several candidate models, MR scheme ---
#' aux.im.mr  <- list(as.matrix(data_samp_stsrs_unif[, c("x1", "x2")]),
#'                    as.matrix(data_samp_stsrs_unif[, "x1", drop = FALSE]))
#' aux.nrm.mr <- list(as.matrix(data_samp_stsrs_unif[, c("x1", "x2")]),
#'                    as.matrix(data_samp_stsrs_unif[, "x2", drop = FALSE]))
#' set.seed(1)
#' fit.mr <- pseudopop.boot.imp(data_samp_stsrs_unif, pop_size_st_imp,
#'                              R.pop = 2, R.samp = 10, stratum = "stratum",
#'                              aux.IM = aux.im.mr, aux.NRM = aux.nrm.mr,
#'                              parameter = "total",
#'                              imputation.method = "MRDRegI",
#'                              bootstrap.method = "MR",
#'                              R.MR.samp = 2, R.MR.pop = 2)
#' fit.mr$boot.var
#'
#' ## --- Pseudo-population Shao-Sitter scheme with regression imputation ---
#' ## (the first matrix of aux.IM is used as the auxiliary information)
#' set.seed(1)
#' fit.ss <- pseudopop.boot.imp(data_samp_stsrs_unif, pop_size_st_imp,
#'                              R.pop = 2, R.samp = 20, stratum = "stratum",
#'                              aux.IM = aux.im, aux.NRM = aux.nrm,
#'                              parameter = "mean",
#'                              imputation.method = "regression",
#'                              bootstrap.method = "Shao.Sitter")
#' fit.ss$boot.var
#'
#' ## --- Data with other column names: give the names through the arguments.
#' ##     An unnamed population.size follows the sorted stratum labels. ---
#' my.data <- data.frame(income    = data_samp_stsrs_unif$study.variable,
#'                       age.group = data_samp_stsrs_unif$stratum,
#'                       responded = data_samp_stsrs_unif$resp.ind)
#' set.seed(1)
#' fit.custom <- pseudopop.boot.imp(my.data, c(400, 500, 600, 500),
#'                                  R.pop = 2, R.samp = 20,
#'                                  study.variable = "income",
#'                                  stratum = "age.group",
#'                                  resp.ind = "responded",
#'                                  aux.IM = aux.im, aux.NRM = aux.nrm,
#'                                  parameter = "mean",
#'                                  imputation.method = "DRDRegI",
#'                                  bootstrap.method = "NRM")
#' fit.custom$boot.var
#'
#' @export
pseudopop.boot.imp <- function(
  data,
  population.size,
  R.pop = 10,
  R.samp = 100,
  study.variable = "study.variable",
  stratum = NULL,
  resp.ind = "resp.ind",
  resp.prob = "resp.prob",
  incl.prob = "pi",
  sample.size = NULL,
  aux.IM,
  aux.NRM,
  parameter = "total",
  imputation.method,
  bootstrap.method,
  survey.design = "SRSWOR",
  n.classes = NULL,
  R.MR.samp = 10,
  R.MR.pop = 5
) {
  ## Check for missing required arguments
  if (
    is.null(data) ||
      is.null(population.size) ||
      is.null(R.pop) ||
      is.null(R.samp) ||
      is.null(study.variable) ||
      is.null(parameter) ||
      is.null(imputation.method) ||
      is.null(bootstrap.method) ||
      is.null(survey.design)
  ) {
    stop(
      "The following arguments must be non-null: data, population.size, R.pop, R.samp, study.variable, parameter, imputation.method, bootstrap.method and survey.design."
    )
  }

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  if (!(study.variable %in% names(data))) {
    stop(
      "The study-variable column '",
      study.variable,
      "' was not found in 'data'."
    )
  }

  ## Validate the response indicators
  has_resp_ind <- !is.null(resp.ind) && resp.ind %in% names(data)
  if (has_resp_ind) {
    if (!all(data[[resp.ind]] %in% c(0, 1))) {
      stop(
        "The response-indicator column '",
        resp.ind,
        "' must contain only 0s and 1s."
      )
    }
    if (any(is.na(data[[study.variable]][data[[resp.ind]] == 1]))) {
      stop(
        "Respondents (units with ",
        resp.ind,
        " = 1) cannot have missing study-variable values."
      )
    }
  }

  respondents_only <- !((has_resp_ind && any(data[[resp.ind]] == 0)) ||
    any(is.na(data[[study.variable]])))

  ## Validate supplied response probabilities (used as-is wherever present)
  has_resp_prob <- !is.null(resp.prob) && resp.prob %in% names(data)
  if (has_resp_prob) {
    if (
      any(is.na(data[[resp.prob]])) ||
        any(data[[resp.prob]] <= 0 | data[[resp.prob]] > 1)
    ) {
      stop(
        "The response-probability column '",
        resp.prob,
        "' must contain response probabilities in (0, 1]."
      )
    }
  }

  if (respondents_only) {
    if (bootstrap.method %in% c("IM", "MR")) {
      stop(
        "The IM and MR bootstrap schemes require full sample data, including nonrespondent rows."
      )
    }
    if (bootstrap.method == "NRM" && !has_resp_prob) {
      stop(
        "The NRM scheme on a respondents-only sample requires known response probabilities (resp.prob = '",
        resp.prob,
        "')."
      )
    }
    if (bootstrap.method == "NRM" && is.null(sample.size)) {
      stop(
        "For a respondents-only sample, the full-sample sizes within strata must be given through 'sample.size'."
      )
    }
  } else if (!is.null(sample.size)) {
    stop(
      "'sample.size' is only used for a respondents-only sample; for a full sample the stratum sample sizes are taken from 'data'."
    )
  }

  ## Validate the parameter of interest
  valid_parameters <- c(
    "total",
    "mean",
    "quartile.25",
    "quartile.50",
    "median",
    "quartile.75"
  )
  if (!(parameter %in% valid_parameters)) {
    stop(
      "The parameter of interest must be one of: 'total', 'mean', 'quartile.25', 'quartile.50' (or 'median') or 'quartile.75'. The default is the population total."
    )
  }
  ## Validate the imputation method
  valid_imputation_methods <- c(
    "MRDRegI",
    "MRCLRHDI",
    "DRDRegI",
    "DRCLRHDI",
    "mean",
    "ratio",
    "hotdeck",
    "regression"
  )
  if (!(imputation.method %in% valid_imputation_methods)) {
    stop(
      "imputation.method must be one of 'MRDRegI', 'MRCLRHDI', 'DRDRegI', 'DRCLRHDI' (for NRM/IM/MR) or 'mean', 'ratio', 'hotdeck', 'regression' (for Shao.Sitter)."
    )
  }

  ## Validate the bootstrap method
  valid_bootstrap_methods <- c(
    "NRM",
    "IM",
    "MR",
    "Shao.Sitter"
  )
  if (!(bootstrap.method %in% valid_bootstrap_methods)) {
    stop(
      "bootstrap.method must be one of: 'NRM' (nonresponse model), 'IM' (imputation model), 'MR' (multiply robust) or 'Shao.Sitter' (Shao-Sitter reimputation bootstrap)."
    )
  }
  shao_imputation_methods <- c("mean", "ratio", "hotdeck", "regression")
  reimputation_methods <- c("MRDRegI", "MRCLRHDI", "DRDRegI", "DRCLRHDI")
  if (bootstrap.method == "Shao.Sitter") {
    if (!(imputation.method %in% shao_imputation_methods)) {
      stop(
        "The Shao.Sitter bootstrap requires imputation.method to be one of: 'mean', 'ratio', 'hotdeck' or 'regression'."
      )
    }
  } else if (!(imputation.method %in% reimputation_methods)) {
    stop(
      "The NRM, IM, and MR bootstrap schemes require imputation.method to be one of: 'MRDRegI', 'MRCLRHDI', 'DRDRegI', or 'DRCLRHDI'."
    )
  }
  if (!(survey.design %in% c("SRSWOR", "CPS"))) {
    stop("survey.design must be one of 'SRSWOR' or 'CPS'.")
  }
  if (
    !is.numeric(population.size) ||
      any(!is.finite(population.size)) ||
      any(population.size <= 0)
  ) {
    stop("Each population size must be a positive numeric value.")
  }
  if (!is.numeric(R.pop) || R.pop %% 1 != 0 || R.pop < 1) {
    stop(
      "R.pop must be a positive integer."
    )
  }
  if (!is.numeric(R.samp) || R.samp %% 1 != 0 || R.samp < 2) {
    stop(
      "R.samp must be an integer >= 2 (a within-population variance needs at least two samples)."
    )
  }
  if (!is.numeric(R.MR.samp) || R.MR.samp %% 1 != 0 || R.MR.samp < 1) {
    stop(
      "R.MR.samp must be a positive integer."
    )
  }
  if (!is.numeric(R.MR.pop) || R.MR.pop %% 1 != 0 || R.MR.pop < 1) {
    stop(
      "R.MR.pop must be a positive integer."
    )
  }

  ## Data Prep ##

  # A single stratum is assumed when no stratum column is given.
  if (is.null(stratum)) {
    if (length(population.size) > 1) {
      stop(
        "The 'stratum' argument must be specified when 'population.size' contains more than one value. Specify the name of the stratum column for a stratified design, or provide a single population size for an unstratified design."
      )
    }
    stratum <- ".stratum"
    data[[stratum]] <- 1
  } else if (!(stratum %in% names(data))) {
    stop(
      "The stratum column '",
      stratum,
      "' was not found in 'data'. Omit the 'stratum' argument for an unstratified design."
    )
  }
  if (!has_resp_ind) {
    resp.ind <- ".resp.ind"
    data[[resp.ind]] <- as.numeric(!is.na(data[[study.variable]]))
  }
  ord <- order(data[[stratum]])
  data <- data[ord, , drop = FALSE]
  if (!is.list(aux.IM) || !is.list(aux.NRM)) {
    stop("aux.IM and aux.NRM must each be a list of covariate matrices.")
  }
  if (
    any(
      vapply(aux.IM, function(a) nrow(as.matrix(a)), integer(1)) != nrow(data)
    ) ||
      any(
        vapply(aux.NRM, function(a) nrow(as.matrix(a)), integer(1)) !=
          nrow(data)
      )
  ) {
    stop(
      "Every matrix in aux.IM and aux.NRM must have one row per unit in 'data' (nrow(data))."
    )
  }
  if (imputation.method == "ratio" && ncol(as.matrix(aux.IM[[1]])) != 1) {
    stop(
      "For ratio imputation, the first matrix of 'aux.IM' must contain a single auxiliary variable."
    )
  }
  aux.IM <- lapply(aux.IM, function(x) as.matrix(x)[ord, , drop = FALSE])
  aux.NRM <- lapply(aux.NRM, function(x) as.matrix(x)[ord, , drop = FALSE])

  if (!"ksi" %in% names(data)) {
    data$ksi <- rep(1, nrow(data))
  }

  strata <- sort(unique(data[[stratum]]))
  L <- length(strata)
  if (L != length(population.size)) {
    stop(
      "The length of 'population.size' must be equal to the number of strata. Each element of 'population.size' is the population size of one stratum."
    )
  }

  # Match the population sizes to the strata: by name for a named vector,
  # otherwise in the sorted order of the stratum labels.
  if (is.null(names(population.size))) {
    names(population.size) <- as.character(strata)
  }
  if (!all(as.character(strata) %in% names(population.size))) {
    stop(
      "The names of 'population.size' must match the stratum labels in the stratum column of the data set."
    )
  }
  strata_population <- as.numeric(population.size[as.character(strata)])
  strata_sample_size <- as.numeric(table(data[[stratum]]))
  if (!is.null(sample.size)) {
    # Full-sample sizes of a respondents-only sample, matched as population.size
    if (
      !is.numeric(sample.size) ||
        length(sample.size) != L ||
        any(sample.size %% 1 != 0)
    ) {
      stop(
        "'sample.size' must contain one integer sample size for each stratum."
      )
    }
    if (is.null(names(sample.size))) {
      names(sample.size) <- as.character(strata)
    }
    if (!all(as.character(strata) %in% names(sample.size))) {
      stop(
        "The names of 'sample.size' must match the stratum labels in the stratum column of the data set."
      )
    }
    respondents_size <- strata_sample_size
    strata_sample_size <- as.numeric(sample.size[as.character(strata)])
    if (any(strata_sample_size < respondents_size)) {
      stop(
        "Each full-sample size in 'sample.size' must be at least the number of respondents in the stratum."
      )
    }
  }
  if (any(strata_population < strata_sample_size)) {
    stop(
      "Each stratum population size must be at least as large as the corresponding sample size."
    )
  }
  respondents_per_stratum <- tapply(
    data[[resp.ind]],
    data[[stratum]],
    function(ri) sum(ri == 1)
  )
  if (any(respondents_per_stratum < 1)) {
    stop("Every stratum must contain at least one respondent.")
  }
  # n.classes is only used by the class hot-deck imputers.
  if (imputation.method %in% c("MRCLRHDI", "DRCLRHDI")) {
    if (
      is.null(n.classes) ||
        !is.numeric(n.classes) ||
        n.classes %% 1 != 0 ||
        n.classes < 1
    ) {
      stop(
        "For class hot-deck imputation, n.classes must be a positive integer."
      )
    }
    if (n.classes > sum(strata_sample_size)) {
      stop(
        "For class hot-deck imputation, n.classes cannot exceed the sample size."
      )
    }
  }
  if (!(incl.prob %in% names((data)))) {
    for (h in seq_along(strata)) {
      data[[incl.prob]][data[[stratum]] == strata[h]] <- strata_sample_size[h] /
        strata_population[h]
    }
  }
  if (
    any(is.na(data[[incl.prob]])) ||
      any(data[[incl.prob]] <= 0 | data[[incl.prob]] > 1)
  ) {
    stop(
      "The inclusion probabilities (incl.prob = '",
      incl.prob,
      "') must be in (0, 1]."
    )
  }

  # Keep data easy to access.
  # `index` holds row positions in the (stratum-ordered) data; it is used to
  # look up rows of aux.IM / aux.NRM, which were reordered the same way.
  orig_data <- list(
    index = seq_len(nrow(data)),
    y = data[[study.variable]],
    resp_ind = data[[resp.ind]],
    w = 1 / data[[incl.prob]],
    pi = data[[incl.prob]],
    ksi = data$ksi,
    stratum = data[[stratum]],
    resp_prob = if (has_resp_prob) data[[resp.prob]] else NULL
  )

  ## Auxiliary/Helper Functions ##
  # Returns A^-1 %*% B with a more detailed message when an error occurs.
  solve_catch <- function(a, b, title) {
    tryCatch(
      solve(a, b),
      error = function(e) {
        stop(sprintf(
          "Could not fit the %s: the system is singular. [%s]",
          title,
          conditionMessage(e)
        ))
      }
    )
  }

  total_ht <- function(y, w, q = NULL) {
    sum(y * w)
  }

  quantile_estimator <- function(y, w, q) {
    y_w <- cbind(y, w / sum(w))[order(y), ]
    i <- 1
    cum_weight <- 0
    while (cum_weight < q) {
      cum_weight <- cum_weight + y_w[i, 2]
      i <- i + 1
    }
    if (cum_weight == q) {
      return(unname((y_w[i - 1, 1] + y_w[i, 1]) / 2))
    } else {
      return(unname(y_w[i - 1, 1]))
    }
  }
  q <- NULL
  if (parameter %in% c("total", "mean")) {
    statistic <- total_ht
  } else if (parameter == "quartile.25") {
    statistic <- quantile_estimator
    q <- 0.25
  } else if (parameter %in% c("quartile.50", "median")) {
    statistic <- quantile_estimator
    q <- 0.50
  } else if (parameter == "quartile.75") {
    statistic <- quantile_estimator
    q <- 0.75
  }
  mean_divisor <- if (parameter == "mean") sum(strata_population) else 1

  # Sampling Designs
  SRSWOR <- function(N, U_pi, n) {
    sample(1:N, size = n, replace = FALSE)
  }
  cps_poisson <- function(N, U_pi, n) {
    epsilon <- runif(N)
    s_index <- (1:N)[epsilon < U_pi]
    if (!is.null(n)) {
      while (length(s_index) != n) {
        epsilon <- runif(N)
        s_index <- (1:N)[epsilon < U_pi]
      }
    }
    return(s_index)
  }
  Design <- if (survey.design == "CPS") cps_poisson else SRSWOR

  ## Imputation Functions ##

  # Deterministic Regression Imputation MR
  impute_det_reg_mr <- function(sr_y, p_hat_MR, y_hat_MR, r_s, w_s) {
    hi <- t(cbind(1, y_hat_MR))

    T_r <- hi[, r_s == 1] %*%
      ((w_s[r_s == 1] * (1 - p_hat_MR[r_s == 1]) / p_hat_MR[r_s == 1]) *
        t(hi[, r_s == 1]))
    t_r <- hi[, r_s == 1] %*%
      ((w_s[r_s == 1] * (1 - p_hat_MR[r_s == 1]) / p_hat_MR[r_s == 1]) * sr_y)

    ksi_hat <- solve_catch(T_r, t_r, "MR deterministic regression imputation")

    imp_y <- t(hi) %*% ksi_hat

    sI_y <- imp_y
    sI_y[r_s == 1] <- sr_y

    list(sI_y = sI_y)
  }

  ## Random hot-deck imputation within C imputation classes of (approximately)
  ## equal size, formed by ordering the sample on `score` (Chen et al., 2019,
  ## 2021). Donors are chosen within a class with probability proportional to
  ## `donor_weight`. The nonrespondents of a class without any donor are set to
  ## NA, and the bootstrap sample is then ignored, as in Chen et al. (2019, 2021).
  impute_classes <- function(sI_y, r_s, score, donor_weight, n_classes) {
    n <- length(r_s)
    donor_weight <- as.vector(donor_weight)
    ordered_units <- order(score)
    ordered_class <- ceiling(seq_len(n) * n_classes / n)
    for (t in seq_len(n_classes)) {
      class_units <- ordered_units[ordered_class == t]
      donors <- class_units[r_s[class_units] == 1]
      recipients <- class_units[r_s[class_units] == 0]
      if (length(recipients) == 0) {
        next
      }
      if (length(donors) == 0 || sum(donor_weight[donors]) <= 0) {
        sI_y[recipients] <- NA
        next
      }
      sI_y[recipients] <- sI_y[donors][suppressWarnings(sample.int(
        length(donors),
        length(recipients),
        replace = TRUE,
        prob = donor_weight[donors] / sum(donor_weight[donors])
      ))]
    }
    sI_y
  }

  ## Multiply robust random hot-deck imputation: classes are formed on the
  ## compressed scores m_i (eq. 9) and donors drawn proportionally to
  ## w_i (p_i^-1 - 1), with p_i given by eq. 8 (Chen et al., 2021, Section 3).
  impute_class_hotdeck_mr <- function(
    sr_y,
    p_hat_MR,
    y_hat_MR,
    r_s,
    w_s,
    n_classes
  ) {
    sI_y <- numeric(length(r_s))
    sI_y[r_s == 1] <- sr_y
    donor_weight <- w_s * ((1 / p_hat_MR) - 1)
    list(sI_y = impute_classes(sI_y, r_s, y_hat_MR, donor_weight, n_classes))
  }

  impute_det_reg_dr <- function(sr_y, s_t, p_hat, r_s, w) {
    resp_ind <- r_s == 1
    sr_t <- s_t[resp_ind, , drop = FALSE]

    phi <- w[resp_ind] * (1 - p_hat[resp_ind]) / p_hat[resp_ind]
    beta_hat <- solve_catch(
      t(sr_t) %*% (phi * sr_t),
      t(sr_t) %*% (phi * sr_y),
      "DR deterministic-regression imputation"
    )
    sI_y <- as.numeric(s_t %*% beta_hat)
    sI_y[resp_ind] <- sr_y
    list(sI_y = sI_y, beta_hat = beta_hat)
  }

  ## Doubly robust random hot-deck imputation: classes are formed on the
  ## predicted values x_i' beta_r, with beta_r the survey-weighted solution of
  ## eq. (6) (phi_i = pi_i^-1), and donors drawn proportionally to
  ## phi_i = pi_i^-1 (p_i^-1 - 1) (Chen et al., 2019, Section 2).
  impute_class_hotdeck_dr <- function(sr_y, p_hat, s_t, r_s, w, n_classes) {
    resp_ind <- r_s == 1
    sr_t <- s_t[resp_ind, , drop = FALSE]
    beta_hat <- lm(sr_y ~ sr_t - 1, weights = w[resp_ind])[[1]]

    y_hat <- s_t %*% beta_hat
    sI_y <- numeric(length(r_s))
    sI_y[resp_ind] <- sr_y
    donor_weight <- w * (1 - p_hat) / p_hat
    list(
      sI_y = impute_classes(sI_y, r_s, y_hat, donor_weight, n_classes),
      beta_hat = beta_hat
    )
  }

  ## Process Functions ##
  # Create pseudo_pop
  generate_pseudo_pop <- function(
    index,
    inv_prob,
    ksi,
    stratum_sample_size,
    N_pop,
    same_pi
  ) {
    np <- length(index)
    if (same_pi) {
      U_boot <- c(
        rep(index, floor(N_pop / np)),
        index[sample(1:np, N_pop - np * floor(N_pop / np), replace = FALSE)]
      )
    } else {
      k <- floor(inv_prob)
      U_boot <- c(rep(index, k), index[cps_poisson(np, inv_prob - k, NULL)])
    }
    # U_boot will be a set of indices to dataset which generated 'index'
    Np <- length(U_boot)
    U_boot_ksi <- ksi[U_boot]
    U_boot_pi <- stratum_sample_size * U_boot_ksi / sum(U_boot_ksi)
    while (any(U_boot_pi > 1)) {
      for (j in seq_len(Np)) {
        if (U_boot_pi[j] >= 1) {
          U_boot_pi[j] <- 1
        } else {
          U_boot_pi[j] <- (stratum_sample_size - sum(U_boot_pi >= 1)) *
            U_boot_ksi[j] /
            sum(U_boot_ksi[(seq_len(Np)[U_boot_pi < 1])])
        }
      }
    }
    list(U = U_boot, pi = U_boot_pi)
  }

  # Calculations for Doubly Robust methods
  calculate_dr <- function(data) {
    resp_ind <- data$resp_ind == 1

    s_t <- cbind(1, aux.IM[[1]][data$index, , drop = FALSE])

    aux_nrm <- aux.NRM[[1]][data$index, , drop = FALSE]

    p_hat <- suppressWarnings(predict(
      glm(data$resp_ind ~ aux_nrm, family = binomial, weights = data$w),
      type = "response"
    ))
    p_hat[p_hat > 1] <- 1
    p_hat[p_hat < 1e-8] <- 1e-8 # guard against Inf weights when p_hat ~ 0

    phi <- data$w[resp_ind] * (1 - p_hat[resp_ind]) / p_hat[resp_ind]
    sr_t <- s_t[resp_ind, , drop = FALSE]
    beta_hat <- solve_catch(
      t(sr_t) %*% (phi * sr_t),
      t(sr_t) %*% (phi * data$y[resp_ind]),
      "DR beta_hat"
    )
    y_hat <- as.numeric(s_t %*% beta_hat)
    sr_e <- data$y[resp_ind] - y_hat[resp_ind]
    sr_e_tilde <- sr_e - mean(sr_e)
    list(
      p_hat = p_hat,
      y_pred = y_hat,
      sr_e_tilde = sr_e_tilde,
      s_t = s_t,
      beta_hat = beta_hat
    )
  }
  # Calculations for Multiply Robust Methods
  calculate_mr <- function(data) {
    resp_ind <- data$resp_ind == 1
    n_IM <- length(aux.IM)
    n_NRM <- length(aux.NRM)
    IM_est <- matrix(NA, length(data$y), n_IM)
    IM_est_resp <- matrix(NA, sum(resp_ind), n_IM)
    NRM_est <- matrix(NA, length(data$y), n_NRM)
    for (j in seq_len(n_IM)) {
      aux_j <- aux.IM[[j]][data$index, , drop = FALSE]
      fit_coef <- coefficients(lm(
        data$y[resp_ind] ~ aux_j[resp_ind, ],
        weights = data$w[resp_ind]
      ))

      IM_est[, j] <- apply(fit_coef * t(cbind(1, aux_j)), 2, sum)
      IM_est_resp[, j] <- apply(
        fit_coef * t(cbind(1, aux_j[resp_ind, ])),
        2,
        sum
      )
    }
    for (j in seq_len(n_NRM)) {
      aux_j <- aux.NRM[[j]][data$index, , drop = FALSE]
      NRM_est[, j] <- suppressWarnings(predict(
        glm(data$resp_ind ~ aux_j, family = binomial, weights = data$w),
        type = "response"
      ))
    }
    NRM_est_resp <- t(as.matrix(NRM_est[resp_ind, ]))
    eta_hat_p <- solve_catch(
      t(NRM_est) %*% (data$w * NRM_est),
      NRM_est_resp %*% (data$w[resp_ind]),
      "MR eta_hat_p"
    )
    IM_est_resp_matrix <- t(as.matrix(IM_est_resp))
    eta_hat_m <- solve_catch(
      IM_est_resp_matrix %*% (data$w[resp_ind] * t(IM_est_resp_matrix)),
      IM_est_resp_matrix %*% (data$w[resp_ind] * data$y[resp_ind]),
      "MR eta_hat_m"
    )

    p_hat <- NRM_est %*%
      (eta_hat_p^2 / as.numeric(t(eta_hat_p) %*% eta_hat_p))
    p_hat[p_hat > 1] <- 1
    p_hat[p_hat < 1e-8] <- 1e-8 # guard against Inf weights when p_hat =~ 0

    m_hat <- IM_est %*%
      (eta_hat_m^2 / as.numeric(t(eta_hat_m) %*% eta_hat_m))

    sr_e <- (data$y[resp_ind] - m_hat[resp_ind])
    sr_e_tilde <- sr_e - mean(sr_e)
    list(p_hat = p_hat, y_pred = m_hat, sr_e_tilde = sr_e_tilde)
  }
  calculate_prereq <- function(data) {
    if (imputation.method %in% c("MRDRegI", "MRCLRHDI")) {
      calculate_mr(data)
    } else if (imputation.method %in% c("DRDRegI", "DRCLRHDI")) {
      calculate_dr(data)
    }
  }

  # Shao Sitter Imputation Methods
  impute_shao_mean <- function(sr_y, s_r, s_w) {
    n <- length(s_r)
    resp_ind <- s_r == 1
    imp_y <- rep(sum(s_w[resp_ind] * sr_y) / sum(s_w[resp_ind]), n)
    imp_y[resp_ind] <- sr_y
    imp_y
  }
  impute_shao_hotdeck <- function(sr_y, s_r, s_w) {
    n <- length(s_r)
    resp_ind <- s_r == 1
    imp_y <- numeric(n)
    imp_y[resp_ind] <- sr_y
    ## Donors are selected with probability proportional to their survey
    ## weights (all equal under SRSWOR within a stratum)
    sr_w <- s_w[resp_ind]
    donor_prob <- if (length(unique(sr_w)) > 1) sr_w / sum(sr_w) else NULL
    imp_y[!resp_ind] <- sr_y[sample.int(
      length(sr_y),
      sum(!resp_ind),
      replace = TRUE,
      prob = donor_prob
    )]
    imp_y
  }
  impute_shao_regression <- function(sr_y, s_r, s_x, s_w) {
    resp_ind <- s_r == 1
    sr_x <- s_x[resp_ind, , drop = FALSE]
    sr_w <- s_w[resp_ind]
    imp_y <- predict(
      lm(sr_y ~ ., as.data.frame(sr_x), weights = sr_w),
      newdata = as.data.frame(s_x)
    )
    imp_y[resp_ind] <- sr_y
    as.vector(imp_y)
  }
  impute_shao_ratio <- function(sr_y, s_x, s_r, s_w) {
    resp_ind <- s_r == 1
    s_x <- as.vector(s_x)
    imp_y <- sum(sr_y * s_w[resp_ind]) * s_x / sum(s_x * s_r * s_w)
    imp_y[resp_ind] <- sr_y
    imp_y
  }

  # Imputation caller function
  impute <- function(data, preds) {
    resp_ind <- data$resp_ind == 1
    if (imputation.method == "MRDRegI") {
      impute_det_reg_mr(
        data$y[resp_ind],
        preds$p_hat,
        preds$y_pred,
        data$resp_ind,
        data$w
      )$sI_y
    } else if (imputation.method == "MRCLRHDI") {
      impute_class_hotdeck_mr(
        data$y[resp_ind],
        preds$p_hat,
        preds$y_pred,
        data$resp_ind,
        data$w,
        n.classes
      )$sI_y
    } else if (imputation.method == "DRDRegI") {
      impute_det_reg_dr(
        data$y[resp_ind],
        preds$s_t,
        preds$p_hat,
        data$resp_ind,
        data$w
      )$sI_y
    } else if (imputation.method == "DRCLRHDI") {
      impute_class_hotdeck_dr(
        data$y[resp_ind],
        preds$p_hat,
        preds$s_t,
        data$resp_ind,
        data$w,
        n.classes
      )$sI_y
    }
  }
  impute_shao <- function(sr_y, s_r, s_w = NULL, s_x = NULL) {
    if (imputation.method == "mean") {
      impute_shao_mean(sr_y, s_r, s_w)
    } else if (imputation.method == "hotdeck") {
      impute_shao_hotdeck(sr_y, s_r, s_w)
    } else if (imputation.method == "regression") {
      impute_shao_regression(sr_y, s_r, s_x, s_w)
    } else if (imputation.method == "ratio") {
      impute_shao_ratio(sr_y, s_x, s_r, s_w)
    }
  }

  ## Bootstrap Schemes
  # Imputed bootstrap estimate; NA (sample ignored) if a missing value could
  # not be imputed because an imputation class had no donor.
  estimate_imputed <- function(sample_b) {
    y_imp <- impute(sample_b, calculate_prereq(sample_b))
    if (anyNA(y_imp)) {
      return(NA_real_)
    }
    statistic(y_imp, sample_b$w, q) / mean_divisor
  }

  bootstrap_nrm <- function(data) {
    # If resp_probs are supplied, use those always, even when supplied full data.
    # Otherwise they are estimated
    p_hat <- data$resp_prob
    if (is.null(p_hat) && !respondents_only) {
      p_hat <- calculate_prereq(data)$p_hat
    }

    est <- matrix(0, R.pop, R.samp)
    for (population in seq_len(R.pop)) {
      U_h <- lapply(seq_len(L), function(h) {
        resp_ind_h <- which(data$stratum == strata[h] & data$resp_ind == 1)
        n_resp_h <- length(resp_ind_h)
        stratum_sample_size <- strata_sample_size[h]
        stratum_population <- strata_population[h]
        same_pi <- all(
          data$pi[data$stratum == strata[h]] ==
            stratum_sample_size / stratum_population
        )
        ki_h <- floor(1 / p_hat[resp_ind_h])
        sf_boot_h <- rep(resp_ind_h, times = ki_h)
        sp_boot <- c(
          sf_boot_h,
          resp_ind_h[cps_poisson(
            n_resp_h,
            1 / p_hat[resp_ind_h] - ki_h,
            NULL
          )]
        )
        # Create a pseudo population from the pseudo samples
        generate_pseudo_pop(
          sp_boot,
          1 / data$pi[sp_boot],
          data$ksi,
          stratum_sample_size,
          stratum_population,
          same_pi
        )
      })
      for (samp in seq_len(R.samp)) {
        sample_indices <- numeric(0)
        sample_resp_ind <- numeric(0)
        sample_weights <- numeric(0)
        for (stratum in seq_len(L)) {
          U_boot <- U_h[[stratum]]$U
          U_boot_pi <- U_h[[stratum]]$pi
          pseudo_pop_size <- length(U_boot)
          stratum_sample_size <- strata_sample_size[stratum]
          # Create a pseudo sample from respondents
          boot_index <- Design(
            pseudo_pop_size,
            U_boot_pi,
            stratum_sample_size
          )
          s_boot <- U_boot[boot_index]
          r_boot <- rbinom(stratum_sample_size, 1, p_hat[s_boot])
          while (sum(r_boot) < 2) {
            r_boot <- rbinom(stratum_sample_size, 1, p_hat[s_boot])
          }
          sample_indices <- c(sample_indices, s_boot)
          sample_resp_ind <- c(sample_resp_ind, r_boot)
          sample_weights <- c(sample_weights, 1 / U_boot_pi[boot_index])
        }
        sample_b <- list(
          index = data$index[sample_indices],
          y = data$y[sample_indices],
          resp_ind = sample_resp_ind,
          w = sample_weights,
          pi = 1 / sample_weights,
          ksi = data$ksi[sample_indices],
          stratum = data$stratum[sample_indices]
        )
        est[population, samp] <- estimate_imputed(sample_b)
      }
    }
    var_a <- apply(est, 1, function(x) var(x, na.rm = TRUE))
    list(boot_var = mean(var_a), boot_statistic = as.vector(est))
  }

  bootstrap_im <- function(data) {
    preds <- calculate_prereq(data)
    est <- matrix(0, R.pop, R.samp)
    theta <- matrix(0, R.pop, R.samp)
    var_parameter_est <- numeric(R.pop)

    for (population in seq_len(R.pop)) {
      U_h <- lapply(seq_len(L), function(stratum) {
        h_index <- which(data$stratum == strata[stratum])
        stratum_population <- strata_population[stratum]
        stratum_sample_size <- strata_sample_size[stratum]
        same_pi <- all(
          data$pi[h_index] == stratum_sample_size / stratum_population
        )
        generate_pseudo_pop(
          h_index,
          1 / data$pi[h_index],
          data$ksi,
          stratum_sample_size,
          stratum_population,
          same_pi
        )
      })
      for (samp in seq_len(R.samp)) {
        pseudo_pop_y <- numeric(0)
        pseudo_pop_weight <- numeric(0)
        sample_indices <- numeric(0)
        sample_y <- numeric(0)
        sample_resp_ind <- numeric(0)
        sample_weights <- numeric(0)

        for (stratum in seq_len(L)) {
          U_boot <- U_h[[stratum]]$U
          U_boot_pi <- U_h[[stratum]]$pi
          pseudo_pop_size <- length(U_boot)
          stratum_sample_size <- strata_sample_size[stratum]

          resid_sample <- sample(
            preds$sr_e_tilde,
            pseudo_pop_size,
            replace = TRUE
          )
          U_y_boot <- preds$y_pred[U_boot] + resid_sample
          pseudo_pop_y <- c(pseudo_pop_y, U_y_boot)
          pseudo_pop_weight <- c(pseudo_pop_weight, rep(1, pseudo_pop_size))

          boot_index <- Design(
            pseudo_pop_size,
            U_boot_pi,
            stratum_sample_size
          )
          s_boot <- U_boot[boot_index]
          sample_indices <- c(sample_indices, s_boot)
          sample_y <- c(sample_y, U_y_boot[boot_index])
          sample_resp_ind <- c(sample_resp_ind, data$resp_ind[s_boot])
          sample_weights <- c(sample_weights, 1 / U_boot_pi[boot_index])
        }
        theta[population, samp] <- statistic(
          pseudo_pop_y,
          pseudo_pop_weight,
          q
        ) /
          mean_divisor

        sample_b <- list(
          index = sample_indices,
          y = sample_y,
          resp_ind = sample_resp_ind,
          w = sample_weights,
          pi = 1 / sample_weights,
          ksi = data$ksi[sample_indices],
          stratum = data$stratum[sample_indices]
        )
        est[population, samp] <- estimate_imputed(sample_b)
      }
      var_parameter_est[population] <- mean(
        (est[population, ] - theta[population, ])^2,
        na.rm = TRUE
      )
    }
    list(
      boot_var = mean(var_parameter_est),
      boot_statistic = as.vector(est)
    )
  }

  bootstrap_mr <- function(data) {
    preds <- calculate_prereq(data)
    nrm_var <- bootstrap_nrm(data)$boot_var
    im_var <- bootstrap_im(data)$boot_var
    E <- numeric(R.MR.pop)
    for (d in seq_len(R.MR.pop)) {
      U_h <- lapply(seq_len(L), function(stratum) {
        h_index <- which(data$stratum == strata[stratum])
        stratum_population <- strata_population[stratum]
        stratum_sample_size <- strata_sample_size[stratum]
        same_pi <- all(
          data$pi[h_index] == stratum_sample_size / stratum_population
        )
        generate_pseudo_pop(
          h_index,
          1 / data$pi[h_index],
          data$ksi,
          stratum_sample_size,
          stratum_population,
          same_pi
        )
      })
      # Step 1: bootstrap values y* = m* + e* of the pseudo-population, kept
      # fixed for the C bootstrap samples of Steps 2-3 (Chen et al., 2021,
      # Section 4.3).
      U_y_h <- lapply(seq_len(L), function(stratum) {
        U_boot <- U_h[[stratum]]$U
        resid_sample <- sample(
          preds$sr_e_tilde,
          length(U_boot),
          replace = TRUE
        )
        preds$y_pred[U_boot] + resid_sample
      })
      V_c <- numeric(R.MR.samp)
      for (c in seq_len(R.MR.samp)) {
        sample_index <- numeric(0)
        sample_y <- numeric(0)
        sample_resp_ind <- numeric(0)
        sample_weight <- numeric(0)

        for (stratum in seq_len(L)) {
          U_boot <- U_h[[stratum]]$U
          U_boot_pi <- U_h[[stratum]]$pi
          pseudo_pop_size <- length(U_boot)
          stratum_sample_size <- strata_sample_size[stratum]
          U_y_boot <- U_y_h[[stratum]]

          boot_index <- Design(pseudo_pop_size, U_boot_pi, stratum_sample_size)
          s_boot <- U_boot[boot_index]
          sample_index <- c(sample_index, s_boot)
          sample_y <- c(sample_y, U_y_boot[boot_index])
          sample_resp_ind <- c(sample_resp_ind, data$resp_ind[s_boot])
          sample_weight <- c(sample_weight, 1 / U_boot_pi[boot_index])
        }
        sample_c <- list(
          index = data$index[sample_index],
          y = sample_y,
          resp_ind = sample_resp_ind,
          w = sample_weight,
          pi = 1 / sample_weight,
          ksi = data$ksi[sample_index],
          stratum = data$stratum[sample_index]
        )
        V_c[c] <- bootstrap_nrm(sample_c)$boot_var
      }
      E[d] <- mean(V_c)
    }
    list(boot_var = nrm_var - mean(E) + im_var, boot_statistic = NULL)
  }

  bootstrap_shao_sitter <- function(data) {
    aux_data <- aux.IM[[1]]

    # Perform imputation on each stratum
    # Imputation functions leave original data as is
    imp_y <- as.vector(data$y)
    for (stratum in seq_len(L)) {
      index <- which(data$stratum == strata[stratum])
      resp_ind <- data$resp_ind[index]
      sr_y <- imp_y[index][resp_ind == 1]
      s_x <- aux_data[index, , drop = FALSE]
      imp_y[index] <- impute_shao(sr_y, resp_ind, data$w[index], s_x)
    }

    est <- matrix(0, R.pop, R.samp)
    var_parameter_est <- numeric(R.pop)

    for (population in seq_len(R.pop)) {
      # Create pseudo population
      U_h <- lapply(seq_len(L), function(stratum) {
        h_index <- which(data$stratum == strata[stratum])
        stratum_population <- strata_population[stratum]
        stratum_sample_size <- strata_sample_size[stratum]
        same_pi <- all(
          data$pi[h_index] == stratum_sample_size / stratum_population
        )
        generate_pseudo_pop(
          h_index,
          1 / data$pi[h_index],
          data$ksi,
          stratum_sample_size,
          stratum_population,
          same_pi
        )
      })
      for (samp in seq_len(R.samp)) {
        sample_y <- numeric(0)
        sample_w <- numeric(0)

        for (stratum in seq_len(L)) {
          U_boot <- U_h[[stratum]]$U
          U_boot_pi <- U_h[[stratum]]$pi
          pseudo_pop_size <- length(U_boot)
          stratum_sample_size <- strata_sample_size[stratum]

          repeat {
            boot_index <- Design(
              pseudo_pop_size,
              U_boot_pi,
              stratum_sample_size
            )
            s_boot <- U_boot[boot_index]
            boot_r <- data$resp_ind[s_boot]
            if (sum(boot_r) >= 1) {
              break
            }
          }

          # Reimpute nonrespondents with the same method.
          boot_w <- 1 / U_boot_pi[boot_index]
          boot_sr_y <- imp_y[s_boot][boot_r == 1]
          boot_x <- aux_data[s_boot, , drop = FALSE]
          boot_imp_y <- impute_shao(boot_sr_y, boot_r, boot_w, boot_x)

          sample_y <- c(sample_y, boot_imp_y)
          sample_w <- c(sample_w, boot_w)
        }
        est[population, samp] <- statistic(sample_y, sample_w, q) /
          mean_divisor
      }
      # Variance of the bootstrap estimates within the pseudo-population
      # (Chen et al., 2019, Algorithm 1, Step 5)
      var_parameter_est[population] <- var(est[population, ])
    }
    list(
      boot_var = mean(var_parameter_est),
      boot_statistic = as.vector(est)
    )
  }

  ## Main Body
  bootstrap_output <- switch(
    bootstrap.method,
    NRM = bootstrap_nrm(orig_data),
    IM = bootstrap_im(orig_data),
    MR = bootstrap_mr(orig_data),
    Shao.Sitter = bootstrap_shao_sitter(orig_data)
  )

  n_ignored <- sum(is.na(bootstrap_output$boot_statistic))
  if (n_ignored > 0) {
    warning(
      n_ignored,
      " bootstrap sample(s) were ignored because an imputation class had no respondent.",
      call. = FALSE
    )
  }

  ## Output
  list(
    boot.statistic = bootstrap_output$boot_statistic,
    boot.var = bootstrap_output$boot_var,
    boot.mean = if (is.null(bootstrap_output$boot_statistic)) {
      NA_real_
    } else {
      mean(bootstrap_output$boot_statistic, na.rm = TRUE)
    },
    sample.size = strata_sample_size,
    parameter = parameter,
    imputation.method = imputation.method,
    bootstrap.method = bootstrap.method,
    number.bootstrap.pop.replicates = R.pop,
    number.bootstrap.samp.replicates = R.samp
  )
}
