u <- data_samp_stsrs_unif
q <- data_samp_stsrs_uneq
aux3 <- u[, c("x1", "x2", "x3")]

fit <- function(data, ..., R = 10, seed = 1) {
  set.seed(seed)
  rescale.boot.imp(data, pop_size_st_imp, R = R, stratum = "stratum", ...)
}

test_that("every bootstrap/imputation combination runs and returns finite results", {
  cases <- list(
    list("independent", "mean", NULL, u), list("independent", "hotdeck", NULL, u),
    list("independent", "ratio", u$x1, u), list("independent", "regression", aux3, u),
    list("independent.weights", "mean", NULL, u), list("independent.weights", "hotdeck", NULL, u),
    list("independent.weights", "ratio", u$x1, u), list("independent.weights", "regression", aux3, u),
    list("Shao.Sitter", "mean", NULL, u), list("Shao.Sitter", "hotdeck", NULL, u),
    list("Shao.Sitter", "ratio", u$x1, u), list("Shao.Sitter", "regression", aux3, u),
    list("unequal.resp.prob", "mean", NULL, q), list("unequal.resp.prob", "hotdeck", NULL, q),
    list("unequal.resp.prob", "regression", q[, c("x1", "x2", "x3")], q)
  )
  for (cs in cases) {
    res <- fit(cs[[4]], aux.data = cs[[3]], imputation.method = cs[[2]],
               bootstrap.method = cs[[1]])
    label <- paste(cs[[1]], cs[[2]])
    expect_true(is.finite(res$boot.var), label = label)
    expect_length(res$boot.statistic, 10)
    expect_false(anyNA(res$boot.statistic), label = label)
  }
})

test_that("output follows the bootsurv conventions", {
  res <- fit(u, parameter = "median", R = 5)
  expect_named(res, c("boot.statistic", "boot.var", "boot.mean", "boot.sample",
                      "sample.size", "parameter", "imputation.method",
                      "bootstrap.method", "number.bootstrap.replicates"))
  expect_equal(res$parameter, "median")
  expect_equal(res$bootstrap.method, "independent")
  expect_equal(res$sample.size, c(40, 100, 180, 80))
  expect_length(res$boot.sample, 5)
  expect_true(all(c("bootstrap.sample", "row.id", "stratum", "bootstrap.weight") %in%
                    names(res$boot.sample[[1]])))
})

test_that("bootstrap weights and rescaled data give identical estimators", {
  ## Mashreghi, Leger and Haziza (2014, Appendix): for a function of means the
  ## independent.weights and independent methods give the same bootstrap
  ## estimator. Using three auxiliaries for regression also guards the
  ## column-wise centering of the auxiliary matrix.
  for (im in c("mean", "hotdeck", "ratio", "regression")) {
    ax <- switch(im, ratio = u$x1, regression = aux3, NULL)
    brs <- fit(u, aux.data = ax, imputation.method = im, parameter = "mean",
               bootstrap.method = "independent", R = 20, seed = 11)
    bw <- fit(u, aux.data = ax, imputation.method = im, parameter = "mean",
              bootstrap.method = "independent.weights", R = 20, seed = 11)
    expect_equal(bw$boot.statistic, brs$boot.statistic, tolerance = 1e-6, label = im)
  }
})

test_that("population.size is matched by name, or by sorted stratum labels if unnamed", {
  ref <- fit(u)
  set.seed(1)
  unnamed <- rescale.boot.imp(u, unname(pop_size_st_imp), R = 10, stratum = "stratum")
  set.seed(1)
  shuffled <- rescale.boot.imp(u, pop_size_st_imp[c(4, 2, 1, 3)], R = 10,
                               stratum = "stratum")
  expect_equal(unnamed$boot.statistic, ref$boot.statistic)
  expect_equal(shuffled$boot.statistic, ref$boot.statistic)
})

test_that("custom column names give the same result as the default names", {
  renamed <- data.frame(income = u$study.variable, age.group = u$stratum,
                        responded = u$resp.ind)
  set.seed(1)
  custom <- rescale.boot.imp(renamed, pop_size_st_imp, R = 10,
                             study.variable = "income", stratum = "age.group",
                             resp.ind = "responded")
  expect_equal(custom$boot.statistic, fit(u)$boot.statistic)
})

test_that("strata with very few respondents are handled", {
  one <- u
  s1 <- which(one$stratum == "18-29")
  one$resp.ind[s1] <- 0
  one$resp.ind[s1[2]] <- 1
  one$study.variable[s1] <- NA
  one$study.variable[s1[2]] <- 500
  ## the independent bootstrap needs nr_h >= 2 (Table 1 involves nr_h - 1)
  expect_error(fit(one, bootstrap.method = "independent"), "at least two respondents")
  ## Shao-Sitter redraws bootstrap samples without respondents
  res <- fit(one, imputation.method = "hotdeck", bootstrap.method = "Shao.Sitter", R = 100)
  expect_false(anyNA(res$boot.statistic))
  none <- one
  none$resp.ind[s1] <- 0
  none$study.variable[s1] <- NA
  expect_error(fit(none, bootstrap.method = "Shao.Sitter"), "no respondents")
})

test_that("an unstratified sample works with a single population size", {
  one <- u[u$stratum == "45-59", c("study.variable", "resp.ind")]
  set.seed(1)
  res <- rescale.boot.imp(one, 600, R = 5)
  expect_true(is.finite(res$boot.var))
})

test_that("invalid inputs give informative errors", {
  expect_error(rescale.boot.imp(u, pop_size_st_imp, R = 5),
               "'stratum' argument must be specified")
  expect_error(rescale.boot.imp(u, pop_size_st_imp, R = 5, stratum = "age"),
               "not found")
  expect_error(fit(u, aux.data = aux3, imputation.method = "ratio"),
               "single auxiliary")
  expect_error(fit(u, imputation.method = "regression"), "aux.data")
  expect_error(fit(u, bootstrap.method = "unequal.resp.prob"), "resp.prob")
  expect_error(fit(q, aux.data = q$x1, imputation.method = "ratio",
                   bootstrap.method = "unequal.resp.prob"),
               "not designed for ratio")
  expect_error(fit(u, parameter = "q1"), "quartile.25")
  expect_error(fit(u, R = 1), "R must")
})
