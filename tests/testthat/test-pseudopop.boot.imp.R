u <- data_samp_stsrs_unif
aux_im <- list(as.matrix(u[, c("x1", "x2", "x3")]))
aux_nrm <- list(as.matrix(u[, c("x1", "x2")]))

run <- function(data, ..., population.size = pop_size_st_imp, stratum = "stratum",
                seed = 1) {
  set.seed(seed)
  pseudopop.boot.imp(data, population.size, R.pop = 2, R.samp = 5,
                     stratum = stratum, aux.IM = aux_im, aux.NRM = aux_nrm, ...)
}

test_that("every scheme and imputation method runs and returns a finite variance", {
  cases <- list(
    list(bootstrap.method = "NRM", imputation.method = "DRDRegI"),
    list(bootstrap.method = "NRM", imputation.method = "MRDRegI"),
    list(bootstrap.method = "IM", imputation.method = "DRCLRHDI", n.classes = 4),
    list(bootstrap.method = "IM", imputation.method = "MRCLRHDI", n.classes = 4),
    list(bootstrap.method = "MR", imputation.method = "DRDRegI",
         R.MR.samp = 2, R.MR.pop = 2),
    list(bootstrap.method = "Shao.Sitter", imputation.method = "mean"),
    list(bootstrap.method = "Shao.Sitter", imputation.method = "hotdeck"),
    list(bootstrap.method = "Shao.Sitter", imputation.method = "regression")
  )
  for (cs in cases) {
    res <- do.call(run, c(list(u, parameter = "mean"), cs))
    expect_true(is.finite(res$boot.var),
                label = paste(cs$bootstrap.method, cs$imputation.method))
  }
})

test_that("output follows the bootsurv conventions", {
  res <- run(u, parameter = "quartile.25", imputation.method = "DRDRegI",
             bootstrap.method = "NRM")
  expect_named(res, c("boot.statistic", "boot.var", "boot.mean", "sample.size",
                      "parameter", "imputation.method", "bootstrap.method",
                      "number.bootstrap.pop.replicates",
                      "number.bootstrap.samp.replicates"))
  expect_equal(res$sample.size, c(40, 100, 180, 80))
  expect_length(res$boot.statistic, 2 * 5)
})

test_that("imputation classes follow Chen et al. (2019, 2021)", {
  ## classes of (approximately) equal size: n need not be a multiple of C
  res <- run(u, parameter = "median", imputation.method = "DRCLRHDI",
             bootstrap.method = "IM", n.classes = 7)
  expect_true(is.finite(res$boot.var))
  ## a bootstrap sample with a class without donors is ignored, with a warning
  expect_warning(
    res <- run(u, parameter = "median", imputation.method = "MRCLRHDI",
               bootstrap.method = "IM", n.classes = 400),
    "ignored"
  )
  expect_true(all(is.na(res$boot.statistic)))
})

test_that("the conditional Poisson design runs", {
  res <- run(u, parameter = "mean", imputation.method = "DRDRegI",
             bootstrap.method = "IM", survey.design = "CPS")
  expect_true(is.finite(res$boot.var))
})

test_that("an ID-valued unit column does not affect the results", {
  ## Regression test: an "index" column of unit IDs used to be used to look up
  ## rows of the auxiliary matrices.
  with_id <- u
  with_id$index <- with_id$unit.id
  a <- run(with_id, parameter = "mean", imputation.method = "DRDRegI",
           bootstrap.method = "NRM")
  b <- run(u, parameter = "mean", imputation.method = "DRDRegI",
           bootstrap.method = "NRM")
  expect_equal(a$boot.var, b$boot.var)
})

test_that("population.size is matched by name, or by sorted stratum labels if unnamed", {
  ref <- run(u, parameter = "mean", imputation.method = "DRDRegI",
             bootstrap.method = "NRM")
  unnamed <- run(u, population.size = unname(pop_size_st_imp), parameter = "mean",
                 imputation.method = "DRDRegI", bootstrap.method = "NRM")
  shuffled <- run(u, population.size = pop_size_st_imp[c(4, 2, 1, 3)],
                  parameter = "mean", imputation.method = "DRDRegI",
                  bootstrap.method = "NRM")
  expect_equal(unnamed$boot.var, ref$boot.var)
  expect_equal(shuffled$boot.var, ref$boot.var)
})

test_that("custom column names give the same result as the default names", {
  renamed <- data.frame(income = u$study.variable, age.group = u$stratum,
                        responded = u$resp.ind)
  default <- run(u, parameter = "mean", imputation.method = "DRDRegI",
                 bootstrap.method = "NRM")
  custom <- run(renamed, stratum = "age.group", study.variable = "income",
                resp.ind = "responded", parameter = "mean",
                imputation.method = "DRDRegI", bootstrap.method = "NRM")
  expect_equal(custom$boot.var, default$boot.var)
})

test_that("a respondents-only sample with sample.size reproduces the full sample", {
  ## With the same respondents, response probabilities and full-sample sizes,
  ## the NRM scheme must give exactly the result obtained from the full sample.
  q <- data_samp_stsrs_uneq
  resp <- q[q$resp.ind == 1, ]
  nrm <- function(d) {
    set.seed(7)
    a <- list(as.matrix(d[, c("x1", "x2")]))
    pseudopop.boot.imp(d, pop_size_st_imp, R.pop = 3, R.samp = 20, stratum = "stratum",
                       aux.IM = a, aux.NRM = a, parameter = "mean",
                       imputation.method = "DRDRegI", bootstrap.method = "NRM",
                       sample.size = if (all(d$resp.ind == 1)) c(40, 100, 180, 80))
  }
  expect_identical(nrm(resp)$boot.var, nrm(q)$boot.var)
  a_r <- list(as.matrix(resp[, c("x1", "x2")]))
  expect_error(
    pseudopop.boot.imp(resp, pop_size_st_imp, stratum = "stratum", aux.IM = a_r,
                       aux.NRM = a_r, imputation.method = "DRDRegI",
                       bootstrap.method = "NRM"),
    "sample.size"
  )
  expect_error(
    pseudopop.boot.imp(resp, pop_size_st_imp, stratum = "stratum", aux.IM = a_r,
                       aux.NRM = a_r, imputation.method = "DRDRegI",
                       bootstrap.method = "IM"),
    "full sample"
  )
})

test_that("invalid inputs give informative errors", {
  expect_error(run(u, imputation.method = "mean", bootstrap.method = "NRM"),
               "NRM, IM, and MR")
  expect_error(run(u, imputation.method = "DRDRegI", bootstrap.method = "Shao.Sitter"),
               "Shao.Sitter")
  expect_error(run(u, imputation.method = "DRCLRHDI", bootstrap.method = "IM",
                   n.classes = 401),
               "cannot exceed")
  expect_error(run(u, study.variable = "income", imputation.method = "DRDRegI",
                   bootstrap.method = "NRM"),
               "not found")
  expect_error(run(u, stratum = NULL, imputation.method = "DRDRegI",
                   bootstrap.method = "NRM"),
               "'stratum' argument must be specified")
  expect_error(run(u, survey.design = "Poisson", imputation.method = "DRDRegI",
                   bootstrap.method = "NRM"),
               "SRSWOR")
  expect_error(run(u, imputation.method = "ratio", bootstrap.method = "Shao.Sitter"),
               "single auxiliary")
  expect_error(run(u, sample.size = c(40, 100, 180, 80), imputation.method = "DRDRegI",
                   bootstrap.method = "NRM"),
               "only used for a respondents-only")
})
