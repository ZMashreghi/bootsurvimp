#' Rescaling Bootstrap for Imputed Survey Data
#'
#' @description
#' Rescaling bootstrap variance estimation for imputed survey data under
#' stratified simple random sampling without replacement. Implements the
#' independent bootstrap of Mashreghi, Leger and Haziza (2014), in both its
#' rescaled-data and bootstrap-weights forms, the
#' unequal-response-probability rescaling bootstrap of Mashreghi and Deng
#' (2023), and the Shao and Sitter (1996) method. Missing values of the study
#' variable are imputed within strata by mean, random hot-deck, regression or
#' ratio imputation, and the imputation is repeated on every bootstrap sample.
#'
#' @param data A vector, matrix or data frame containing the sample. For a
#' matrix or data frame, the study-variable column is identified through the
#' `study.variable` argument (by default, the column named `study.variable`)
#' and, under stratified sampling, the stratum column through the `stratum`
#' argument. When a vector is supplied, its values are taken as the study
#' variable. The study variable may contain missing values (`NA`) to be imputed
#' within the function, or it may be already completed with observed and
#' imputed values; in the latter case a response-indicator column (see
#' `resp.ind`) must identify the respondents.
#'
#' @param population.size A vector of strata population sizes. If it is a named
#' vector, its names must match the stratum labels; if it is unnamed, its
#' elements must follow the sorted order of the stratum labels (as in the
#' `bootsurv` package). For an unstratified design, a single population size.
#'
#' @param R The number of bootstrap replicates. The default is `200`.
#'
#' @param study.variable A character string giving the name of the column in
#' `data` that holds the study-variable values. The default is
#' `"study.variable"`.
#'
#' @param stratum A character string giving the name of the column in `data`
#' that identifies the strata under STSRSWOR. The default is `NULL`, in which
#' case a single stratum is assumed (SRSWOR). For a stratified design, set
#' `stratum` to the name of the column in `data` that identifies the strata.
#'
#' @param resp.ind A character string giving the name of the column in `data`
#' that holds the response indicators (1 = respondent, 0 = nonrespondent). The
#' default is `"resp.ind"`. If the column is absent, respondents are identified
#' by the non-missing values of the study variable. It is required when the
#' study variable is already completed, for the `"Shao.Sitter"` method, and for
#' hot-deck imputation under the `"unequal.resp.prob"` method.
#'
#' @param resp.prob A character string giving the name of the column in `data`
#' that holds the response probabilities. The default is `"resp.prob"`. It is
#' required for the `"unequal.resp.prob"` method.
#'
#' @param aux.data Auxiliary information used for regression and ratio
#' imputation, given as a vector, matrix or data frame whose rows correspond to
#' the rows of `data` in the same order. Ratio imputation requires a single
#' auxiliary variable.
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
#' @param imputation.method The imputation method used to handle the missing
#' values. One of: `"mean"` (mean imputation), `"hotdeck"` (random hot-deck
#' imputation), `"regression"` (regression imputation) or `"ratio"` (ratio
#' imputation). The `"ratio"` method cannot be used with the
#' `"unequal.resp.prob"` bootstrap method. The default is `"mean"`.
#'
#' @param bootstrap.method One of the following bootstrap methods can be
#' applied:
#' \itemize{
#'   \item `"independent"`: independent rescaling bootstrap under homogeneous response probabilities within strata, rescaled-data approach (Mashreghi, Léger and Haziza, 2014).
#'   \item `"independent.weights"`: independent rescaling bootstrap under homogeneous response probabilities within strata, bootstrap-weights approach (Mashreghi, Léger and Haziza, 2014).
#'   \item `"unequal.resp.prob"`: rescaling bootstrap for imputed data under unequal response probabilities (Mashreghi and Deng, 2023).
#'   \item `"Shao.Sitter"`: Shao and Sitter (1996) bootstrap, which keeps the original response indicators.
#' }
#' The default is `"independent"`.
#'
#' @param boot.sample.size A vector of bootstrap sample sizes within strata.
#' The length of this vector must be equal to the number of strata. If it is
#' not specified, the bootstrap sample size within each stratum will be
#' `nh-1`, where `nh` is the original sample size in stratum `h`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`boot.statistic`}{A vector of bootstrap statistics.}
#'   \item{`boot.var`}{The bootstrap variance estimator of the estimator of the parameter of interest.}
#'   \item{`boot.mean`}{The average of the bootstrap estimator of the parameter of interest.}
#'   \item{`boot.sample`}{A list of size `R`. Each element is a data frame with the bootstrap sample. For the rescaled-data methods it contains the reimputed bootstrap values (`bootstrap.sample`), the original-sample row id of each unit (`row.id`), the stratum (`stratum`), the bootstrap weights (`bootstrap.weight`), the bootstrap response indicators (`boot.resp.ind`) and, when used, the bootstrap response probabilities (`boot.resp.prob`) and rescaled auxiliary variables (`boot.aux.*`). For `"independent.weights"` it contains the original imputed values (`study.variable`), `row.id`, `stratum`, `bootstrap.weight`, the bootstrap multiplicities (`boot.multiplicity`, `boot.resp.multiplicity`), the response information (`resp.ind`, `resp.prob`) and, when used, the rescaled auxiliary variables (`rescaled.aux.*`).}
#'   \item{`sample.size`}{The number of sampled units in each stratum.}
#'   \item{`parameter`}{The parameter of interest.}
#'   \item{`imputation.method`}{The imputation method.}
#'   \item{`bootstrap.method`}{The bootstrap method.}
#'   \item{`number.bootstrap.replicates`}{The number of bootstrap replicates, `R`.}
#' }
#'
#' @seealso [pseudopop.boot.imp()] for pseudo-population bootstrap methods
#'   under doubly and multiply robust imputation.
#'
#' @references
#'
#' Mashreghi, Z. and Deng, H. (2023). A Rescaling Bootstrap Approach for Imputed Survey Data. Journal of Survey Statistics and Methodology, 11(1), 234-259.
#'
#' Mashreghi, Z., Haziza, D. and Léger, C. (2016). A Survey of Bootstrap Methods in Finite Population Sampling. Statistics Surveys, 10, 1-52.
#'
#' Mashreghi, Z., Léger, C. and Haziza, D. (2014). Bootstrap Methods for Imputed Data from Regression, Ratio and Hot-Deck Imputation. The Canadian Journal of Statistics, 42(1), 142-167.
#'
#' Rao, J. N. K., Wu, C. F. J. and Yue, K. (1992). Some recent work on resampling methods for complex surveys. Survey Methodology 18 (2), 209–217.
#'
#' Särndal, C.-E., Swensson, B. and Wretman, J. (1992). Model-Assisted Survey Sampling. NewYork: Springer.
#'
#' Shao, J. and Sitter, R.R. (1996). Bootstrap for Imputed Survey Data. Journal of the American Statistical Association, 91(435), 1278-1288.
#'
#' @examples
#' ## Stratified SRS without replacement (overall 20% sampling fraction;
#' ## strata = age groups). The package ships:
#' ##   data_samp_stsrs_unif - sample with uniform nonresponse within strata
#' ##   data_samp_stsrs_uneq - same sample, unequal response probabilities
#' ##   pop_size_st_imp      - stratum (age-group) population sizes
#' ##   data_pop_st_imp      - the full population (used below to show the truth)
#' ## R = 50 replicates are used here for speed; use more (e.g. 500+) in practice.
#'
#' ## --- Independent bootstrap, mean imputation, population total ---
#' fit <- rescale.boot.imp(data_samp_stsrs_unif, pop_size_st_imp, R = 50,
#'                         stratum = "stratum")
#' fit$boot.var              # bootstrap variance estimator
#' fit$boot.mean             # average of the bootstrap replicate estimates
#' head(fit$boot.statistic)  # the individual bootstrap replicate estimates
#'
#' ## How close is the estimate to the true population total?
#' c(estimate = fit$boot.mean, true = sum(data_pop_st_imp$study.variable))
#'
#' ## --- A strong auxiliary lowers the variance: compare mean imputation
#' ##     with regression imputation (population mean). set.seed() makes the
#' ##     two runs use the same bootstrap draws, so the variances are comparable.
#' aux <- data_samp_stsrs_unif[, c("x1", "x2", "x3")]
#' set.seed(1)
#' fit.mean <- rescale.boot.imp(data_samp_stsrs_unif, pop_size_st_imp, R = 50,
#'                              stratum = "stratum", parameter = "mean")
#' set.seed(1)
#' fit.reg <- rescale.boot.imp(data_samp_stsrs_unif, pop_size_st_imp, R = 50,
#'                             stratum = "stratum", aux.data = aux,
#'                             parameter = "mean",
#'                             imputation.method = "regression")
#' c(mean.imp = fit.mean$boot.var, regression = fit.reg$boot.var)
#'
#' ## --- Ratio imputation with the single positive auxiliary x1 ---
#' fit.ratio <- rescale.boot.imp(data_samp_stsrs_unif, pop_size_st_imp, R = 50,
#'                               stratum = "stratum",
#'                               aux.data = data_samp_stsrs_unif$x1,
#'                               imputation.method = "ratio")
#' fit.ratio$boot.var
#'
#' ## --- Bootstrap-weights version of the independent bootstrap ---
#' fit.bw <- rescale.boot.imp(data_samp_stsrs_unif, pop_size_st_imp, R = 50,
#'                            stratum = "stratum",
#'                            bootstrap.method = "independent.weights")
#' fit.bw$boot.var
#'
#' ## --- Shao-Sitter bootstrap, population median, hot-deck imputation ---
#' fit.ss <- rescale.boot.imp(data_samp_stsrs_unif, pop_size_st_imp, R = 50,
#'                            stratum = "stratum", parameter = "median",
#'                            imputation.method = "hotdeck",
#'                            bootstrap.method = "Shao.Sitter")
#' fit.ss$boot.var
#'
#' ## --- Unequal response probabilities (Mashreghi and Deng, 2023);
#' ##     uses the resp.prob column of data_samp_stsrs_uneq ---
#' fit.uneq <- rescale.boot.imp(data_samp_stsrs_uneq, pop_size_st_imp, R = 50,
#'                              stratum = "stratum",
#'                              bootstrap.method = "unequal.resp.prob")
#' fit.uneq$boot.var
#'
#' ## --- Data with other column names: give the names through the arguments.
#' ##     An unnamed population.size follows the sorted stratum labels. ---
#' my.data <- data.frame(income    = data_samp_stsrs_unif$study.variable,
#'                       age.group = data_samp_stsrs_unif$stratum,
#'                       responded = data_samp_stsrs_unif$resp.ind)
#' fit.custom <- rescale.boot.imp(my.data, c(400, 500, 600, 500), R = 50,
#'                                study.variable = "income",
#'                                stratum = "age.group",
#'                                resp.ind = "responded")
#' fit.custom$boot.var
#'
#' @export
rescale.boot.imp <- function(
  data,
  population.size,
  R = 200,
  study.variable = "study.variable",
  stratum = NULL,
  resp.ind = "resp.ind",
  resp.prob = "resp.prob",
  aux.data = NULL,
  parameter = "total",
  imputation.method = "mean",
  bootstrap.method = "independent",
  boot.sample.size = NULL
) {
  ############################## Validate the input arguments ##############################

  ## Check for missing required arguments
  if (
    is.null(data) ||
      is.null(population.size) ||
      is.null(R) ||
      is.null(parameter) ||
      is.null(imputation.method) ||
      is.null(bootstrap.method)
  ) {
    stop(
      "The following arguments must be non-null: data, population.size, R, parameter, imputation.method and bootstrap.method."
    )
  }

  ## Validate the number of bootstrap replicates
  if (!is.numeric(R) || length(R) != 1 || R %% 1 != 0 || R < 2) {
    stop(
      "R must be an integer >= 2 (a bootstrap variance needs at least two replicates)."
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
  valid_imputation_methods <- c("mean", "hotdeck", "regression", "ratio")
  if (!(imputation.method %in% valid_imputation_methods)) {
    stop(
      "The imputation method must be one of: 'mean', 'hotdeck', 'regression' or 'ratio'."
    )
  }

  ## Validate the bootstrap method
  valid_bootstrap_methods <- c(
    "independent",
    "independent.weights",
    "unequal.resp.prob",
    "Shao.Sitter"
  )
  if (!(bootstrap.method %in% valid_bootstrap_methods)) {
    stop(
      "The bootstrap method must be one of: 'independent' or 'independent.weights' (Mashreghi, Leger and Haziza, 2014), 'unequal.resp.prob' (Mashreghi and Deng, 2023) or 'Shao.Sitter' (Shao and Sitter, 1996)."
    )
  }

  ## Validate the compatibility between the bootstrap and imputation methods
  if (bootstrap.method == "unequal.resp.prob" && imputation.method == "ratio") {
    stop(
      "The 'unequal.resp.prob' bootstrap method of Mashreghi and Deng (2023) is not designed for ratio imputation."
    )
  }

  ## Convert the data and auxiliary data to a data frame
  if (is.vector(data)) {
    data <- data.frame(data)
    names(data) <- study.variable
  }
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  if (!is.null(aux.data)) {
    if (is.vector(aux.data)) {
      aux.data <- data.frame(aux.variable = aux.data)
    }
    if (is.matrix(aux.data)) aux.data <- as.data.frame(aux.data)
  }

  ## Validate the study variable
  if (!(study.variable %in% names(data))) {
    stop(
      "The study-variable column '",
      study.variable,
      "' was not found in 'data'."
    )
  }
  if (!is.numeric(data[[study.variable]])) {
    stop("The study variable must be numeric.")
  }

  ## Validate the stratification information
  if (is.null(stratum)) {
    if (length(population.size) > 1) {
      stop(
        "The 'stratum' argument must be specified when 'population.size' contains more than one value. Specify the name of the stratum column for a stratified design, or provide a single population size for an unstratified (SRSWOR) design."
      )
    }
    stratum_values <- rep(1, nrow(data))
  } else {
    if (!(stratum %in% names(data))) {
      stop(
        "The stratum column '",
        stratum,
        "' was not found in 'data'. Omit the 'stratum' argument for an unstratified (SRSWOR) design."
      )
    }
    stratum_values <- data[[stratum]]
  }

  has_resp_ind <- !is.null(resp.ind) && resp.ind %in% names(data)
  has_resp_prob <- !is.null(resp.prob) && resp.prob %in% names(data)

  ## Store the sample under the internal column names
  ## (index is the original-sample row id of each unit)
  sample_data <- data.frame(
    study_variable = data[[study.variable]],
    stratum = stratum_values,
    index = seq_len(nrow(data)),
    stringsAsFactors = FALSE
  )
  if (has_resp_ind) {
    sample_data$resp_ind <- data[[resp.ind]]
  }
  if (has_resp_prob) {
    sample_data$resp_prob <- data[[resp.prob]]
  }
  data <- sample_data

  ## Validate the response indicators
  if (has_resp_ind) {
    if (!all(data$resp_ind %in% c(0, 1))) {
      stop(
        "The response-indicator column '",
        resp.ind,
        "' must contain only 0s and 1s."
      )
    }
    if (any(is.na(data$study_variable[data$resp_ind == 1]))) {
      stop(
        "Respondents (units with ",
        resp.ind,
        " = 1) cannot have missing study-variable values."
      )
    }
  } else {
    data$resp_ind <- as.numeric(!is.na(data$study_variable))
  }

  ## Validate the response probabilities
  if (has_resp_prob) {
    if (
      any(is.na(data$resp_prob)) ||
        any(data$resp_prob <= 0 | data$resp_prob > 1)
    ) {
      stop(
        "The response-probability column '",
        resp.prob,
        "' must contain response probabilities in (0, 1]."
      )
    }
  }

  ## Validate the auxiliary information
  if (!is.null(aux.data) && nrow(aux.data) != nrow(data)) {
    stop(
      "'aux.data' must contain the same number of observations as 'data', with rows corresponding to the same sampled units in the same order."
    )
  }

  if (imputation.method == "ratio") {
    if (is.null(aux.data)) {
      stop(
        "Auxiliary information must be provided through 'aux.data' when the imputation method is 'ratio'."
      )
    }
    if (ncol(aux.data) != 1) {
      stop(
        "For ratio imputation, 'aux.data' must contain a single auxiliary variable."
      )
    }
  }

  if (imputation.method == "regression" && is.null(aux.data)) {
    stop(
      "Auxiliary information must be provided through 'aux.data' when the imputation method is 'regression'."
    )
  }

  ## Validate the population sizes
  L <- length(unique(data$stratum))

  if (
    !is.numeric(population.size) ||
      any(!is.finite(population.size)) ||
      any(population.size <= 0)
  ) {
    stop("Each population size must be a positive numeric value.")
  }

  if (length(population.size) != L) {
    stop(
      "The length of 'population.size' must be equal to the number of strata. Each element of 'population.size' is the population size of one stratum."
    )
  }

  ## Match the population sizes to the strata: by name for a named vector,
  ## otherwise in the sorted order of the stratum labels
  if (is.null(names(population.size))) {
    names(population.size) <- as.character(sort(unique(data$stratum)))
  } else if (!all(as.character(unique(data$stratum)) %in% names(population.size))) {
    stop(
      "The names of 'population.size' must match the stratum labels in the stratum column of the data set."
    )
  }

  stratum_sizes <- table(data$stratum)
  if (any(population.size[names(stratum_sizes)] < as.numeric(stratum_sizes))) {
    stop(
      "Each stratum population size must be at least as large as the corresponding sample size."
    )
  }

  ## Validate the columns required by the bootstrap methods
  if (bootstrap.method == "unequal.resp.prob") {
    if (!has_resp_prob) {
      stop(
        "The 'unequal.resp.prob' bootstrap method of Mashreghi and Deng (2023) requires a column of response probabilities in 'data' (resp.prob = '",
        resp.prob,
        "')."
      )
    }
    if (imputation.method == "hotdeck" && !has_resp_ind) {
      stop(
        "Hot-deck imputation with the 'unequal.resp.prob' bootstrap method of Mashreghi and Deng (2023) requires a column of response indicators in 'data' (resp.ind = '",
        resp.ind,
        "')."
      )
    }
  }

  ## Validate the number of respondents within strata
  resp_per_stratum <- tapply(data$resp_ind, data$stratum, sum)
  if (any(resp_per_stratum == 0)) {
    stop(
      "Stratum ",
      names(resp_per_stratum)[resp_per_stratum == 0][1],
      " contains no respondents."
    )
  }
  if (
    bootstrap.method %in% c("independent", "independent.weights") &&
      any(resp_per_stratum < 2)
  ) {
    stop(
      "The independent bootstrap requires at least two respondents in each stratum, since its rescaling factor involves nr_h - 1 (Mashreghi, Leger and Haziza, 2014, Table 1)."
    )
  }

  if (bootstrap.method == "Shao.Sitter" && !has_resp_ind) {
    stop(
      "The 'Shao.Sitter' bootstrap method of Shao and Sitter (1996) requires a column of response indicators in 'data' (resp.ind = '",
      resp.ind,
      "')."
    )
  }

  ############################## Define auxiliary functions ##############################

  ## Imputation functions
  ## Notation:
  ## sr_y : respondent values of the study variable
  ## s_x  : auxiliary values for the full sample
  ## s_r  : response indicators for the sample
  ## s_w  : survey weights for the sample
  ## imp_y: imputed values of the study variable

  # Mean imputation
  mean_impute <- function(sr_y, s_r, s_w) {
    n <- length(s_r) #sample size
    imp_y <- rep(sum(s_w[s_r == 1] * sr_y) / sum(s_w[s_r == 1]), n) #impute each y
    imp_y[s_r == 1] <- sr_y #use observed values when r=1
    return(imp_y)
  }

  # Random hot-deck imputation
  hotdeck_impute <- function(sr_y, s_r) {
    n <- length(s_r) #sample size
    imp_y <- rep(0, n) #n-vector of 0s
    imp_y[s_r == 1] <- sr_y #use observed values when r=1
    imp_y[s_r == 0] <- sr_y[sample.int(length(sr_y), sum(s_r == 0), replace = TRUE)] #imputed values
    return(imp_y)
  }

  # Regression imputation
  regression_impute <- function(sr_y, s_x, s_r, s_w) {
    sr_x <- s_x[s_r == 1, , drop = FALSE]
    sr_w <- s_w[s_r == 1]
    model <- lm(sr_y ~ ., data = as.data.frame(sr_x), weights = sr_w)
    imp_y <- predict(model, newdata = as.data.frame(s_x)) #imputed values
    imp_y[s_r == 1] <- sr_y #use observed values when r=1
    imp_y <- as.vector(imp_y)
    return(imp_y)
  }

  # Ratio imputation
  ratio_impute <- function(sr_y, s_x, s_r, s_w) {
    s_x <- as.vector(s_x)
    imp_y <- sum(sr_y * s_w[s_r == 1]) * s_x / sum(s_x * s_r * s_w) #imputed values
    imp_y[s_r == 1] <- sr_y #use observed values when r=1
    return(imp_y)
  }

  # Imputation with the selected imputation method
  impute <- function(sr_y, s_r, s_w, s_x) {
    if (imputation.method == "mean") {
      return(mean_impute(sr_y, s_r, s_w))
    }
    if (imputation.method == "hotdeck") {
      return(hotdeck_impute(sr_y, s_r))
    }
    if (imputation.method == "regression") {
      return(regression_impute(sr_y, s_x, s_r, s_w))
    }
    if (imputation.method == "ratio") {
      return(ratio_impute(sr_y, s_x, s_r, s_w))
    }
  }

  ## Parameter estimation functions
  # Horvitz-Thompson estimator of the population total
  total_ht <- function(y, w, q = NULL) sum(y * w)

  ## Quantile estimator
  ## See Särndal, Swensson, and Wretman (1992, Chapter 5)
  quantile_estimator <- function(y, w, q) {
    y_w <- cbind(y, w / sum(w))
    y_w_order <- y_w[order(y), ]
    i <- 1
    cum_weight <- 0
    while (cum_weight < q) {
      cum_weight <- cum_weight + y_w_order[i, 2]
      i <- i + 1
    }
    if (cum_weight == q) {
      return((y_w_order[i - 1, 1] + y_w_order[i, 1]) / 2)
    }

    return(y_w_order[i - 1, 1])
  }

  ## Define the estimator based on the selected parameter
  q <- NULL

  if (parameter == "total" || parameter == "mean") {
    statistic <- total_ht
  }

  if (parameter == "quartile.25") {
    statistic <- quantile_estimator
    q <- 0.25
  }

  if (parameter == "quartile.50" || parameter == "median") {
    statistic <- quantile_estimator
    q <- 0.50
  }

  if (parameter == "quartile.75") {
    statistic <- quantile_estimator
    q <- 0.75
  }

  ############################## Prepare the survey data ##############################

  ## Order the data by strata
  order_index <- order(data$stratum)

  data <- data[order_index, , drop = FALSE]

  ## Reorder the auxiliary information
  if (!is.null(aux.data)) {
    aux.data <- aux.data[order_index, , drop = FALSE]
    x <- as.matrix(aux.data)
  }

  ## Reorder the population sizes according to the strata order
  strata <- unique(data$stratum)
  population.size <- population.size[as.character(strata)]

  ## Define the total population size
  N <- sum(population.size)

  ## Define the stratum sample sizes
  sample_size <- as.numeric(table(data$stratum))

  ## Define the respondent sample sizes
  sample_resp_size <- as.numeric(tapply(data$resp_ind, data$stratum, sum))

  ## Define the survey weights and response probabilities
  data$weight <- NA
  if (!("resp_prob" %in% names(data))) {
    data$resp_prob <- NA
  }

  for (h in 1:L) {
    stratum_index <- which(data$stratum == strata[h])
    data$weight[stratum_index] <- population.size[h] / sample_size[h]

    if (all(is.na(data$resp_prob[stratum_index]))) {
      data$resp_prob[stratum_index] <- rep(
        sample_resp_size[h] / sample_size[h],
        sample_size[h]
      )
    }
  }

  ############################## Perform the imputation ##############################

  ## Initialize the completed study variable
  imp_y <- as.vector(data$study_variable)

  for (h in 1:L) {
    ## Extract the stratum indices
    stratum_index <- which(data$stratum == strata[h])

    ## Extract the study variable within stratum h
    s_y <- imp_y[stratum_index]

    ## Extract the response indicators within stratum h
    s_r <- as.numeric(data$resp_ind[stratum_index])

    ## Extract the survey weights within stratum h
    s_w <- as.numeric(data$weight[stratum_index])

    ## Extract the auxiliary information within stratum h
    if (!is.null(aux.data)) {
      s_x <- x[stratum_index, , drop = FALSE]
    }

    ## Extract the respondent values within stratum h
    if (sum(s_r == 1) == 0) {
      stop(paste("Stratum", strata[h], "contains no respondents."))
    }

    sr_y <- s_y[s_r == 1]

    ## Perform the imputation within stratum h
    imp_y[stratum_index] <- impute(sr_y, s_r, s_w, s_x)
  }

  ##################### Define the bootstrap rescaling factors and rescaled variables #####################

  ## Initialize the rescaled quantities
  data$rescaled_y <- NA
  if (!is.null(aux.data)) {
    rescaled_x <- matrix(NA, nrow = nrow(data), ncol = ncol(x))
  }

  ## Define the bootstrap sample sizes
  if (is.null(boot.sample.size)) {
    boot.sample.size <- sample_size - 1
  }

  ## Initialize the bootstrap rescaling factors
  data$rescale_factor <- NA

  ## Define the numerical stabilization constant
  epsilon <- 1e-8

  ## Flag whether the imputation method uses the auxiliary information
  uses_aux_data <- !is.null(aux.data) &&
    imputation.method %in% c("regression", "ratio")

  ############################################################
  ## Helper functions
  ############################################################

  ## Extract the survey quantities used by the rescaling within stratum h
  rescale_setup <- function(h, data) {
    ## Extract the stratum indices
    stratum_index <- which(data$stratum == strata[h])

    ## Extract the survey quantities within stratum h
    s_y <- imp_y[stratum_index]
    s_r <- as.numeric(data$resp_ind[stratum_index])
    s_p <- as.numeric(data$resp_prob[stratum_index])
    s_w <- as.numeric(data$weight[stratum_index])
    s_index <- data$index[stratum_index]

    n_h <- sample_size[h]
    nr_h <- sample_resp_size[h]
    N_h <- population.size[h]
    n_boot_h <- boot.sample.size[h]

    f_h <- n_h / N_h

    ## Define the average response probabilities
    p_bar_h <- mean(s_p)
    q_bar_h <- 1 - p_bar_h

    ## Center the study variable
    y_diff <- s_y - mean(s_y)

    ## Stabilize near-zero centered values
    y_diff[abs(y_diff) < epsilon] <- epsilon

    ## Extract and center the auxiliary information
    s_x <- NULL
    x_diff <- NULL
    if (!is.null(aux.data)) {
      s_x <- x[stratum_index, , drop = FALSE]
      x_diff <- sweep(s_x, 2, colMeans(s_x))
    }

    list(
      stratum_index = stratum_index,
      s_y = s_y,
      s_r = s_r,
      s_p = s_p,
      s_w = s_w,
      s_index = s_index,
      n_h = n_h,
      nr_h = nr_h,
      N_h = N_h,
      n_boot_h = n_boot_h,
      f_h = f_h,
      p_bar_h = p_bar_h,
      q_bar_h = q_bar_h,
      y_diff = y_diff,
      s_x = s_x,
      x_diff = x_diff
    )
  }

  ## Rescaling factor and rescaled variables of the independent bootstrap
  ## Takes in list returned by rescale_setup()
  rescale_indep <- function(rescale_h) {
    u_h <- NULL

    ## Mean imputation
    if (imputation.method == "mean") {
      C_h <- sqrt(
        (1 - rescale_h$p_bar_h * rescale_h$f_h) *
          rescale_h$n_boot_h /
          (rescale_h$nr_h - 1)
      )
      z_h <- mean(rescale_h$s_y) + C_h * rescale_h$y_diff
    }

    ## Random hot-deck imputation
    if (imputation.method == "hotdeck") {
      C_h <- sqrt(
        (rescale_h$n_boot_h * rescale_h$p_bar_h / (rescale_h$nr_h - 1)) *
          ((1 - rescale_h$p_bar_h * rescale_h$f_h) +
            rescale_h$p_bar_h * rescale_h$q_bar_h) /
          ((1 - rescale_h$q_bar_h / rescale_h$n_h) *
            (1 + rescale_h$p_bar_h * rescale_h$q_bar_h))
      )
      z_h <- mean(rescale_h$s_y) + C_h * rescale_h$y_diff
    }

    ## Regression imputation
    if (imputation.method == "regression") {
      Sxx <- t(rescale_h$x_diff) %*% rescale_h$x_diff / (rescale_h$n_h - 1)
      Sxy <- t(rescale_h$x_diff) %*% rescale_h$y_diff / (rescale_h$n_h - 1)
      CCI2_indep <- as.numeric(
        t(Sxy) %*% solve(Sxx + epsilon * diag(ncol(Sxx))) %*% Sxy
      ) /
        var(rescale_h$s_y)

      C_h <- sqrt(
        (1 -
          rescale_h$p_bar_h * rescale_h$f_h -
          (rescale_h$p_bar_h *
            rescale_h$q_bar_h *
            CCI2_indep /
            (1 + rescale_h$q_bar_h * CCI2_indep))) *
          rescale_h$n_boot_h /
          (rescale_h$nr_h - 1)
      )

      z_h <- mean(rescale_h$s_y) + C_h * rescale_h$y_diff
      u_h <- sweep(C_h * rescale_h$x_diff, 2, colMeans(rescale_h$s_x), "+")
    }

    ## Ratio imputation
    if (imputation.method == "ratio") {
      Sxx <- sum(rescale_h$x_diff * rescale_h$x_diff) / (rescale_h$n_h - 1)
      Syy <- sum(rescale_h$y_diff^2) / (rescale_h$n_h - 1)
      Sxy <- sum(rescale_h$x_diff * rescale_h$y_diff) / (rescale_h$n_h - 1)

      rho <- Sxy / sqrt(Sxx * Syy)
      R_h <- (sqrt(Sxx) / mean(rescale_h$s_x)) / (sqrt(Syy) / mean(rescale_h$s_y))

      C_h <- sqrt(
        (1 -
          (rescale_h$p_bar_h *
            rescale_h$f_h +
            (rescale_h$p_bar_h *
              (1 - rescale_h$p_bar_h) *
              (1 - rescale_h$f_h) *
              R_h^2)) /
            (1 + ((1 - rescale_h$p_bar_h) * R_h * (R_h - 2 * rho)))) *
          rescale_h$n_boot_h /
          (rescale_h$nr_h - 1)
      )

      z_h <- mean(rescale_h$s_y) + C_h * rescale_h$y_diff
      u_h <- mean(rescale_h$s_x) + C_h * rescale_h$x_diff
    }

    list(C_h = C_h, z_h = z_h, u_h = u_h)
  }

  ## Extract the survey and rescaled quantities used by a bootstrap draw within stratum h
  replicate_setup <- function(h, data, rescaled_x) {
    ## Extract the stratum indices
    stratum_index <- which(data$stratum == strata[h])

    ## Extract the survey quantities within stratum h
    s_y <- imp_y[stratum_index]
    s_r <- as.numeric(data$resp_ind[stratum_index])
    s_p <- as.numeric(data$resp_prob[stratum_index])
    s_w <- as.numeric(data$weight[stratum_index])
    s_index <- data$index[stratum_index]
    s_x <- NULL
    if (uses_aux_data) {
      s_x <- x[stratum_index, , drop = FALSE]
    }

    n_h <- sample_size[h]
    N_h <- population.size[h]
    n_boot_h <- boot.sample.size[h]

    ## Extract the rescaled study variable
    z_h <- data$rescaled_y[stratum_index]

    ## Extract the rescaled auxiliary information
    u_h <- NULL
    if (uses_aux_data) {
      u_h <- rescaled_x[stratum_index, , drop = FALSE]
    }

    list(
      stratum_index = stratum_index,
      s_y = s_y,
      s_r = s_r,
      s_p = s_p,
      s_w = s_w,
      s_index = s_index,
      s_x = s_x,
      n_h = n_h,
      N_h = N_h,
      n_boot_h = n_boot_h,
      z_h = z_h,
      u_h = u_h
    )
  }

  ## Draw a bootstrap sample and regenerate its response indicators
  # Takes in list returned by replicate_setup()
  draw_regenerated <- function(replicate_h) {
    ## Generate the bootstrap sample
    boot_index <- sample(
      1:replicate_h$n_h,
      size = replicate_h$n_boot_h,
      replace = TRUE
    )
    id_boot <- replicate_h$s_index[boot_index]

    ## Generate the bootstrap response indicators
    s_p_boot <- replicate_h$s_p[boot_index]
    r_boot <- rbinom(replicate_h$n_boot_h, 1, s_p_boot)

    while (sum(r_boot) <= 1) {
      r_boot <- rbinom(replicate_h$n_boot_h, 1, s_p_boot)
    }

    list(
      boot_index = boot_index,
      id_boot = id_boot,
      s_p_boot = s_p_boot,
      r_boot = r_boot
    )
  }

  ## Rescaled data approach: reimpute the drawn rescaled values within stratum h
  resample_rescaled <- function(h, replicate_h, draw) {
    ## Extract the bootstrap quantities
    z_boot <- replicate_h$z_h[draw$boot_index]

    u_boot <- NULL
    if (uses_aux_data) {
      u_boot <- replicate_h$u_h[draw$boot_index, , drop = FALSE]
    }

    ## Extract the respondent values
    zr_boot <- z_boot[draw$r_boot == 1]

    ## Define the bootstrap survey weights
    w_boot <- rep(replicate_h$N_h / replicate_h$n_boot_h, replicate_h$n_boot_h)

    ## Perform the bootstrap imputation
    z_imp_boot <- impute(
      zr_boot,
      draw$r_boot,
      rep(replicate_h$n_h / replicate_h$n_boot_h, replicate_h$n_boot_h),
      u_boot
    )

    ## Create the bootstrap sample
    return_new <- data.frame(
      bootstrap.sample = z_imp_boot,
      row.id = draw$id_boot,
      stratum = rep(strata[h], replicate_h$n_boot_h),
      bootstrap.weight = w_boot,
      boot.resp.ind = draw$r_boot,
      boot.resp.prob = draw$s_p_boot
    )

    ## Append the auxiliary information
    if (uses_aux_data) {
      u_boot <- as.data.frame(u_boot)
      colnames(u_boot) <- paste0("boot.aux.", colnames(replicate_h$s_x))
      return_new <- cbind(return_new, u_boot)
    }

    return_new
  }

  ############################## Bootstrap schemes ##############################

  ############################################################
  ## Shao-Sitter bootstrap
  ############################################################

  bootstrap_shao_sitter <- function(data) {
    ## Rescaling
    ## Compute the bootstrap rescaling quantities within strata
    for (h in 1:L) {
      rescale_h <- rescale_setup(h, data)

      C_h <- sqrt(
        (1 - rescale_h$f_h) * rescale_h$n_boot_h / (rescale_h$n_h - 1)
      )
      data$rescale_factor[rescale_h$stratum_index] <- C_h

      ## Rescale the study variable
      z_h <- mean(rescale_h$s_y) + C_h * rescale_h$y_diff

      ## Rescale the auxiliary information
      if (imputation.method == "regression") {
        u_h <- sweep(C_h * rescale_h$x_diff, 2, colMeans(rescale_h$s_x), "+")
      }
      if (imputation.method == "ratio") {
        u_h <- mean(rescale_h$s_x) + C_h * rescale_h$x_diff
      }

      ## Store the rescaled study variable
      data$rescaled_y[rescale_h$stratum_index] <- as.vector(z_h)

      ## Store the rescaled auxiliary information
      if (uses_aux_data) {
        rescaled_x[rescale_h$stratum_index, ] <- as.matrix(u_h)
      }
    }

    ## Bootstrap
    ## Initialize the bootstrap objects
    return_bootstrap <- vector("list", R)
    boot_statistic <- rep(0, R)

    ## Generate the bootstrap samples
    for (b in 1:R) {
      ## Initialize the bootstrap sample
      bootstrap_sample <- NULL

      ## Generate the bootstrap sample within strata
      for (h in 1:L) {
        replicate_h <- replicate_setup(h, data, rescaled_x)

        ## Generate the bootstrap sample indices, preserving the original
        ## response indicators (redrawn if the sample has no respondent)
        repeat {
          boot_index <- sample(
            1:replicate_h$n_h,
            size = replicate_h$n_boot_h,
            replace = TRUE
          )
          r_boot <- replicate_h$s_r[boot_index]
          if (sum(r_boot) >= 1) break
        }

        ## Extract the bootstrap quantities
        z_boot <- replicate_h$z_h[boot_index]
        id_boot <- replicate_h$s_index[boot_index]

        u_boot <- NULL
        if (uses_aux_data) {
          u_boot <- replicate_h$u_h[boot_index, , drop = FALSE]
        }

        ## Extract the respondent values
        zr_boot <- z_boot[r_boot == 1]

        ## Perform the bootstrap imputation
        z_imp_boot <- impute(
          zr_boot,
          r_boot,
          rep(replicate_h$n_h / replicate_h$n_boot_h, replicate_h$n_boot_h),
          u_boot
        )

        ## Define the bootstrap weights
        w_boot <- rep(
          replicate_h$N_h / replicate_h$n_boot_h,
          replicate_h$n_boot_h
        )

        ## Create the bootstrap sample
        return_new <- data.frame(
          bootstrap.sample = z_imp_boot,
          row.id = id_boot,
          stratum = rep(strata[h], replicate_h$n_boot_h),
          bootstrap.weight = w_boot,
          boot.resp.ind = r_boot
        )

        ## Append the auxiliary information
        if (uses_aux_data) {
          u_boot <- as.data.frame(u_boot)
          colnames(u_boot) <- paste0("boot.aux.", colnames(replicate_h$s_x))
          return_new <- cbind(return_new, u_boot)
        }

        ## Append the stratum contribution
        bootstrap_sample <- rbind(bootstrap_sample, return_new)
      }

      ## Store the complete bootstrap replicate
      return_bootstrap[[b]] <- bootstrap_sample

      ## Compute the bootstrap replicate statistic
      boot_statistic[b] <- statistic(
        y = bootstrap_sample$bootstrap.sample,
        w = bootstrap_sample$bootstrap.weight,
        q = q
      )
      ## Convert the bootstrap total to a bootstrap mean
      if (parameter == "mean") boot_statistic[b] <- boot_statistic[b] / N
    }

    list(boot_statistic = boot_statistic, bootstrap_samples = return_bootstrap)
  }

  ############################################################
  ## Independent bootstrap (rescaled data approach)
  ############################################################

  bootstrap_indep <- function(data) {
    ## Rescaling
    ## Compute the bootstrap rescaling quantities within strata
    for (h in 1:L) {
      rescale_h <- rescale_setup(h, data)
      rescaled_h <- rescale_indep(rescale_h)

      data$rescale_factor[rescale_h$stratum_index] <- rescaled_h$C_h

      ## Store the rescaled study variable
      data$rescaled_y[rescale_h$stratum_index] <- as.vector(rescaled_h$z_h)

      ## Store the rescaled auxiliary information
      if (uses_aux_data) {
        rescaled_x[rescale_h$stratum_index, ] <- as.matrix(rescaled_h$u_h)
      }
    }

    ## Bootstrap
    ## Initialize the bootstrap objects
    return_bootstrap <- vector("list", R)
    boot_statistic <- rep(0, R)

    ## Generate the bootstrap samples
    for (b in 1:R) {
      ## Initialize the bootstrap sample
      bootstrap_sample <- NULL

      ## Generate the bootstrap sample within strata
      for (h in 1:L) {
        replicate_h <- replicate_setup(h, data, rescaled_x)
        draw <- draw_regenerated(replicate_h)

        ## Append the stratum contribution
        bootstrap_sample <- rbind(
          bootstrap_sample,
          resample_rescaled(h, replicate_h, draw)
        )
      }

      ## Store the complete bootstrap replicate
      return_bootstrap[[b]] <- bootstrap_sample

      ## Compute the bootstrap replicate statistic
      boot_statistic[b] <- statistic(
        y = bootstrap_sample$bootstrap.sample,
        w = bootstrap_sample$bootstrap.weight,
        q = q
      )
      ## Convert the bootstrap total to a bootstrap mean
      if (parameter == "mean") boot_statistic[b] <- boot_statistic[b] / N
    }

    list(boot_statistic = boot_statistic, bootstrap_samples = return_bootstrap)
  }

  ############################################################
  ## Independent bootstrap (bootstrap weights approach)
  ############################################################

  bootstrap_indep_weight <- function(data) {
    ## Rescaling
    ## Compute the bootstrap rescaling quantities within strata
    for (h in 1:L) {
      rescale_h <- rescale_setup(h, data)
      rescaled_h <- rescale_indep(rescale_h)

      data$rescale_factor[rescale_h$stratum_index] <- rescaled_h$C_h

      ## Store the rescaled study variable
      data$rescaled_y[rescale_h$stratum_index] <- as.vector(rescaled_h$z_h)

      ## Store the rescaled auxiliary information
      if (uses_aux_data) {
        rescaled_x[rescale_h$stratum_index, ] <- as.matrix(rescaled_h$u_h)
      }
    }

    ## Bootstrap
    ## Initialize the bootstrap objects
    return_bootstrap <- vector("list", R)
    boot_statistic <- rep(0, R)

    ## Generate the bootstrap samples
    for (b in 1:R) {
      ## Initialize the bootstrap sample
      bootstrap_sample <- NULL

      ## Generate the bootstrap sample within strata
      for (h in 1:L) {
        replicate_h <- replicate_setup(h, data, rescaled_x)
        draw <- draw_regenerated(replicate_h)

        ## Independent bootstrap method (bootstrap weights approach)

        ## Extract the bootstrap multiplicities
        m_boot <- tabulate(draw$boot_index, nbins = replicate_h$n_h)

        ## Respondent multiplicities
        R_boot <- tabulate(
          draw$boot_index[draw$r_boot == 1],
          nbins = replicate_h$n_h
        )
        R_boot_total <- sum(R_boot)

        ## Extract the bootstrap rescaling factor
        C_h <- unique(data$rescale_factor[replicate_h$stratum_index])

        a_boot <- 1 +
          C_h * (replicate_h$n_h * m_boot / replicate_h$n_boot_h - 1)
        c_boot <- 1 + C_h * (replicate_h$n_h * R_boot / R_boot_total - 1)

        ## Bootstrap weights for mean imputation
        if (imputation.method == "mean") {
          w_boot <- c_boot * replicate_h$s_w
        }

        ## Bootstrap weights for ratio imputation
        if (imputation.method == "ratio") {
          ratio_adj <- sum(a_boot * replicate_h$s_w * replicate_h$s_x) /
            sum(c_boot * replicate_h$s_w * replicate_h$s_x)
          w_boot <- ratio_adj * c_boot * replicate_h$s_w
        }

        ## Bootstrap weights for regression imputation
        if (imputation.method == "regression") {
          ## Extract the rescaled auxiliary information
          u_h <- rescaled_x[replicate_h$stratum_index, , drop = FALSE]

          ## Define the bootstrap selection matrix
          B_star <- matrix(
            0,
            nrow = replicate_h$n_h,
            ncol = replicate_h$n_boot_h
          )
          for (j in 1:replicate_h$n_boot_h) {
            B_star[draw$boot_index[j], j] <- 1
          }

          ## Define the bootstrap respondents matrix
          Bhr_star <- B_star %*% diag(draw$r_boot)

          ## Define the regression adjustment vector k* (CJS 2014, Section 4.3),
          ## using the augmented x-tilde = (1, x) and u-tilde = (1, u):
          ut_h <- cbind(1, u_h)
          xt_h <- cbind(1, replicate_h$s_x)
          k_boot <- as.vector(
            ((1 / replicate_h$N_h) *
              t(t(xt_h) %*% ((a_boot - c_boot) * replicate_h$s_w))) %*%
              solve(
                t(ut_h) %*%
                  Bhr_star %*%
                  t(Bhr_star) %*%
                  ut_h +
                  epsilon * diag(ncol(ut_h))
              ) %*%
              t(ut_h) %*%
              Bhr_star %*%
              t(Bhr_star)
          )

          one_k <- sum(k_boot)

          w_boot <- (c_boot +
            one_k +
            C_h * (replicate_h$n_h * k_boot - one_k)) *
            replicate_h$s_w
        }

        ## Random hot-deck imputation
        if (imputation.method == "hotdeck") {
          donor_boot <- sample(
            which(draw$r_boot == 1),
            size = sum(draw$r_boot == 0),
            replace = TRUE
          )
          D_boot <- tabulate(
            draw$boot_index[donor_boot],
            nbins = replicate_h$n_h
          )
          w_boot <- (1 +
            C_h *
              (replicate_h$n_h *
                (R_boot + D_boot) /
                replicate_h$n_boot_h -
                1)) *
            replicate_h$s_w
        }

        ## Create the bootstrap sample
        return_new <- data.frame(
          study.variable = replicate_h$s_y,
          row.id = replicate_h$s_index,
          stratum = rep(strata[h], replicate_h$n_h),
          bootstrap.weight = w_boot,
          boot.multiplicity = m_boot,
          boot.resp.multiplicity = R_boot,
          resp.ind = replicate_h$s_r,
          resp.prob = replicate_h$s_p
        )

        ## Append the auxiliary information
        if (uses_aux_data) {
          rescaled_aux <- as.data.frame(replicate_h$u_h)
          colnames(rescaled_aux) <- paste0(
            "rescaled.aux.",
            colnames(replicate_h$s_x)
          )
          return_new <- cbind(return_new, rescaled_aux)
        }
        ## Append the stratum contribution
        bootstrap_sample <- rbind(bootstrap_sample, return_new)
      }

      ## Store the complete bootstrap replicate
      return_bootstrap[[b]] <- bootstrap_sample

      ## Compute the bootstrap replicate statistic
      boot_statistic[b] <- statistic(
        y = bootstrap_sample$study.variable,
        w = bootstrap_sample$bootstrap.weight,
        q = q
      )
      ## Convert the bootstrap total to a bootstrap mean
      if (parameter == "mean") boot_statistic[b] <- boot_statistic[b] / N
    }

    list(boot_statistic = boot_statistic, bootstrap_samples = return_bootstrap)
  }

  ############################################################
  ## Unequal response probability bootstrap
  ############################################################

  bootstrap_uneq_resp_prob <- function(data) {
    ## Rescaling
    ## Compute the bootstrap rescaling quantities within strata
    for (h in 1:L) {
      rescale_h <- rescale_setup(h, data)

      p_frac <- rescale_h$s_p / sum(rescale_h$s_p)

      A_h <- matrix(
        rep(-p_frac, rescale_h$n_h),
        nrow = rescale_h$n_h,
        byrow = TRUE
      )
      diag(A_h) <- 1 - p_frac

      ## Mean imputation
      if (imputation.method == "mean") {
        C_i <- sqrt(
          rescale_h$n_boot_h *
            (1 - rescale_h$s_p * rescale_h$f_h) /
            (rescale_h$n_h * rescale_h$s_p)
        )

        b_i <- C_i * rescale_h$y_diff
        z_h <- ginv(A_h + epsilon * diag(rescale_h$n_h)) %*% b_i
      }

      ## Random hot-deck imputation
      if (imputation.method == "hotdeck") {
        H_i <- (1 -
          rescale_h$s_p * rescale_h$f_h +
          rescale_h$p_bar_h * (1 - rescale_h$p_bar_h)) -
          (1 / rescale_h$n_h) *
            sum(
              (1 -
                rescale_h$s_p * rescale_h$f_h +
                rescale_h$p_bar_h * (1 - rescale_h$p_bar_h)) *
                (1 -
                  rescale_h$s_r * (1 - 2 / rescale_h$n_h) -
                  (1 + rescale_h$nr_h / rescale_h$n_h) / rescale_h$n_h)
            )

        C_i <- sqrt(
          (1 / (rescale_h$n_h * rescale_h$s_p)) *
            rescale_h$n_boot_h *
            (1 / (1 + rescale_h$p_bar_h * (1 - rescale_h$p_bar_h))) *
            H_i
        )

        b_i <- C_i * rescale_h$y_diff
        z_h <- ginv(A_h + epsilon * diag(rescale_h$n_h)) %*% b_i
      }

      ## Regression imputation
      if (imputation.method == "regression") {
        Sxy <- t(rescale_h$x_diff * rescale_h$s_w) %*% rescale_h$y_diff
        Sxx <- t(rescale_h$x_diff * rescale_h$s_w) %*% rescale_h$x_diff
        Syy <- sum(rescale_h$s_w * rescale_h$y_diff^2)
        CCI2_uneq <- as.numeric(
          t(Sxy) %*% solve(Sxx + epsilon * diag(ncol(Sxx))) %*% Sxy
        ) /
          Syy
        beta_hat <- solve(Sxx + epsilon * diag(ncol(Sxx))) %*% Sxy

        K_i <- (rescale_h$p_bar_h^(-2) *
          (1 - rescale_h$s_p * rescale_h$f_h) *
          (1 - rescale_h$x_diff %*% beta_hat / rescale_h$y_diff)^2) +
          ((1 - rescale_h$f_h) * CCI2_uneq)

        L_i <- (rescale_h$n_h *
          rescale_h$s_p *
          (rescale_h$p_bar_h^(-2) *
            (1 - rescale_h$x_diff %*% beta_hat / rescale_h$y_diff)^2 +
            CCI2_uneq / rescale_h$p_bar_h)) /
          rescale_h$n_boot_h

        C_i <- sqrt(K_i / L_i)
        b_i <- C_i * rescale_h$y_diff
        d_i <- rescale_h$x_diff * as.vector(C_i)

        z_h <- ginv(A_h + epsilon * diag(rescale_h$n_h)) %*% b_i
        u_h <- ginv(A_h + epsilon * diag(rescale_h$n_h)) %*% d_i
      }

      ## Store the rescaled study variable
      data$rescaled_y[rescale_h$stratum_index] <- as.vector(z_h)

      ## Store the rescaled auxiliary information
      if (uses_aux_data) {
        rescaled_x[rescale_h$stratum_index, ] <- as.matrix(u_h)
      }
    }

    ## Bootstrap

    ## Initialize the bootstrap objects
    return_bootstrap <- vector("list", R)
    boot_statistic <- rep(0, R)

    ## Generate the bootstrap samples
    for (b in 1:R) {
      ## Initialize the bootstrap sample
      bootstrap_sample <- NULL

      ## Generate the bootstrap sample within strata
      for (h in 1:L) {
        replicate_h <- replicate_setup(h, data, rescaled_x)
        draw <- draw_regenerated(replicate_h)

        ## Append the stratum contribution
        bootstrap_sample <- rbind(
          bootstrap_sample,
          resample_rescaled(h, replicate_h, draw)
        )
      }

      ## Store the complete bootstrap replicate
      return_bootstrap[[b]] <- bootstrap_sample

      ## Compute the bootstrap replicate statistic
      boot_statistic[b] <- statistic(
        y = bootstrap_sample$bootstrap.sample,
        w = bootstrap_sample$bootstrap.weight,
        q = q
      )
      ## Convert the bootstrap total to a bootstrap mean
      if (parameter == "mean") boot_statistic[b] <- boot_statistic[b] / N
    }

    list(boot_statistic = boot_statistic, bootstrap_samples = return_bootstrap)
  }

  ## Output

  bootstrap_output <- switch(
    bootstrap.method,
    Shao.Sitter = bootstrap_shao_sitter(data),
    independent = bootstrap_indep(data),
    independent.weights = bootstrap_indep_weight(data),
    unequal.resp.prob = bootstrap_uneq_resp_prob(data)
  )
  boot_statistic <- bootstrap_output$boot_statistic
  return_bootstrap <- bootstrap_output$bootstrap_samples

  ################ Compute the bootstrap mean and variance  ################
  if (any(is.na(boot_statistic))) {
    warning("Some bootstrap statistics could not be computed.")
  }

  boot_mean <- mean(boot_statistic)
  boot_var <- var(boot_statistic)

  ################ Return the bootstrap results  ################

  results <- list(
    boot.statistic = boot_statistic,
    boot.var = boot_var,
    boot.mean = boot_mean,
    boot.sample = return_bootstrap,
    sample.size = sample_size,
    parameter = parameter,
    imputation.method = imputation.method,
    bootstrap.method = bootstrap.method,
    number.bootstrap.replicates = R
  )

  return(results)
}
