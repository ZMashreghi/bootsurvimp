#' bootsurvimp: Bootstrap Variance Estimation for Imputed Survey Data
#'
#' Bootstrap methods for estimating the variance of estimators computed on
#' imputed survey data under stratified sampling designs.
#'
#' It is the companion of the `bootsurv` package, which provides bootstrap
#' methods for complete survey data, and follows the same naming conventions.
#'
#' The package provides two main functions:
#' \describe{
#'   \item{[rescale.boot.imp()]}{Rescaling bootstrap methods: the independent
#'     bootstrap of Mashreghi, Leger and Haziza (2014) (rescaled-data and
#'     bootstrap-weights forms), the unequal-response-probability bootstrap of
#'     Mashreghi and Deng (2023), and the Shao and Sitter (1996) method, under
#'     mean, random hot-deck, regression and ratio imputation.}
#'   \item{[pseudopop.boot.imp()]}{Pseudo-population bootstrap methods under
#'     doubly and multiply robust imputation (nonresponse-model,
#'     imputation-model and multiply robust schemes), plus a pseudo-population
#'     Shao-Sitter scheme.}
#' }
#'
#' Example data sets are described in [bootsurvimp.data].
#'
#' @keywords internal
#' @importFrom MASS ginv
#' @importFrom stats binomial coefficients glm lm predict rbinom runif var
"_PACKAGE"
