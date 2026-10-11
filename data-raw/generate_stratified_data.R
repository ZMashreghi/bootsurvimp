## =====================================================================
## Data generation for the bootsurvimp package
## Stratified simple random sampling without replacement (SSRSWOR)
##
## ILLUSTRATIVE SCENARIO (so a first-time user can read the data easily):
## A synthetic survey of N = 2000 individuals. The population is stratified
## by AGE GROUP, i.e. the `stratum` column is the respondent's age band:
##   "18-29", "30-44", "45-59", "60+".
##
## Variables:
##   study.variable : the variable of interest, subject to nonresponse
##                    (read it as, e.g., annual income). It is NA for
##                    nonrespondents and is imputed by the bootstrap functions.
##   x1             : a strong, fully observed POSITIVE auxiliary
##                    (read it as, e.g., an administrative/tax income on
##                    file). Usable for ratio AND regression imputation.
##   x2, x3         : further fully observed covariates, used as extra
##                    predictors for regression imputation.
##   resp.ind       : response indicator (1 = respondent, 0 = nonrespondent).
##   resp.prob      : (unequal-response sample only) the response probability.
##
## Style follows the user's population generator
## (Population-Sigma50-Skewed-f9-f45-fr70-V5.txt): a study variable built
## from a linear model with a skewed positive auxiliary and N(0, 50) errors,
## plus a logistic response model.
##
## Produces (saved as .rda in data/):
##   data_pop_st_imp       - the full stratified population (with y and auxiliaries)
##   pop_size_st_imp       - named vector of stratum (age-group) population sizes
##   data_samp_stsrs_unif  - the drawn sample, UNIFORM nonresponse within strata
##   data_samp_stsrs_uneq  - the SAME sample, UNEQUAL response probabilities (logistic)
## =====================================================================

set.seed(2023)

## ------------------------- Design parameters -------------------------
## Strata = age groups (illustrative).
strata_labels <- c("18-29", "30-44", "45-59", "60+")
L  <- length(strata_labels)
Nh <- c("18-29" = 400, "30-44" = 500, "45-59" = 600, "60+" = 500)  # stratum sizes
N  <- sum(Nh)                                                       # N = 2000

## OVERALL sampling fraction = 20%  ->  total sample size n = 0.20 * N = 400.
## The 400 units are allocated across the age strata in a NON-proportional
## way, so the per-stratum sampling fractions f_h DIFFER (here 0.10, 0.20,
## 0.30, 0.16). Only the OVERALL fraction equals 20%. (Varying f_h is the
## regime where the independent bootstrap of Mashreghi, Leger & Haziza (2014)
## is designed to stay valid, including for large f_h.)
f_overall <- 0.20
n_total   <- round(f_overall * N)                                  # 400
nh <- c("18-29" = 40, "30-44" = 100, "45-59" = 180, "60+" = 80)    # sum = 400
stopifnot(sum(nh) == n_total)
## resulting per-stratum fractions f_h = nh / Nh = 0.10, 0.20, 0.30, 0.16

## ------------------------- Population model --------------------------
## study variable: y = alpha + b1 x1 + b2 x2 + b3 x3 + age-group effect + eps
## (skewed errors of scale 50, in the spirit of the provided population).

alpha     <- 250
b1        <- 15
b2        <- 5
b3        <- 10
strat_eff <- c("18-29" = 0, "30-44" = 60, "45-59" = 120, "60+" = 200) # age shift
x1_scale  <- c("18-29" = 3, "30-44" = 4, "45-59" = 5, "60+" = 6)      # Gamma scale

pop_list <- vector("list", L)

for (h in seq_len(L)) {
  lab  <- strata_labels[h]
  Nh_h <- Nh[lab]

  x1 <- rgamma(Nh_h, shape = 4, scale = x1_scale[lab])   # positive, skewed
  x2 <- rnorm(Nh_h, mean = 10, sd = 3)
  x3 <- rnorm(Nh_h, mean = 0,  sd = 2)
  eps <- rnorm(Nh_h, mean = 0, sd = 50)

  y <- alpha + b1 * x1 + b2 * x2 + b3 * x3 + strat_eff[lab] + eps

  pop_list[[h]] <- data.frame(
    stratum        = lab,
    x1             = x1,
    x2             = x2,
    x3             = x3,
    study.variable = y,
    stringsAsFactors = FALSE
  )
}

data_pop_st_imp <- do.call(rbind, pop_list)
data_pop_st_imp$unit.id <- seq_len(nrow(data_pop_st_imp))   # population unit id
data_pop_st_imp <- data_pop_st_imp[, c("unit.id", "stratum",
                                         "x1", "x2", "x3", "study.variable")]

## Named population-size vector (the 'population.size' argument of the bootstrap functions)
pop_size_st_imp <- Nh

## ------------------- Draw one SSRSWOR sample (overall 20%) -----------
sample_idx <- unlist(lapply(strata_labels, function(lab) {
  units <- data_pop_st_imp$unit.id[data_pop_st_imp$stratum == lab]
  sort(sample(units, size = nh[lab], replace = FALSE))
}))

strat_sample <- data_pop_st_imp[data_pop_st_imp$unit.id %in% sample_idx, ]
rownames(strat_sample) <- NULL

## ------------------------------------------------------------------
## Nonresponse mechanism 1: UNIFORM within strata
## (constant response probability per age group -> independent / independent.weights /
##  Shao.Sitter methods)
## ------------------------------------------------------------------
p_unif <- c("18-29" = 0.75, "30-44" = 0.65, "45-59" = 0.70, "60+" = 0.60)

resp_ind_unif <- rbinom(nrow(strat_sample), size = 1,
                        prob = p_unif[strat_sample$stratum])

data_samp_stsrs_unif <- data.frame(
  unit.id        = strat_sample$unit.id,
  stratum        = strat_sample$stratum,
  x1             = strat_sample$x1,
  x2             = strat_sample$x2,
  x3             = strat_sample$x3,
  study.variable = strat_sample$study.variable,
  resp.ind       = resp_ind_unif,
  stringsAsFactors = FALSE
)
## set the nonrespondents' study variable to NA (to be imputed by the bootstrap functions)
data_samp_stsrs_unif$study.variable[data_samp_stsrs_unif$resp.ind == 0] <- NA

## ------------------------------------------------------------------
## Nonresponse mechanism 2: UNEQUAL response probabilities
## logit(p_i) = g0 + g1 * z1_i + g2 * z2_i, depending on the auxiliaries
## (-> unequal.resp.prob method of Mashreghi and Deng, 2023)
## ------------------------------------------------------------------
g0 <- 0.6
g1 <- 0.05
g2 <- 0.10
z1 <- strat_sample$x1 - mean(strat_sample$x1)
z2 <- strat_sample$x2 - mean(strat_sample$x2)

resp_prob_uneq <- 1 / (1 + exp(-(g0 + g1 * z1 + g2 * z2)))
resp_ind_uneq  <- rbinom(nrow(strat_sample), size = 1, prob = resp_prob_uneq)

data_samp_stsrs_uneq <- data.frame(
  unit.id        = strat_sample$unit.id,
  stratum        = strat_sample$stratum,
  x1             = strat_sample$x1,
  x2             = strat_sample$x2,
  x3             = strat_sample$x3,
  study.variable = strat_sample$study.variable,
  resp.ind       = resp_ind_uneq,
  resp.prob      = resp_prob_uneq,
  stringsAsFactors = FALSE
)
data_samp_stsrs_uneq$study.variable[data_samp_stsrs_uneq$resp.ind == 0] <- NA

## ------------------------------- Save -------------------------------
## (Datasets to be added to the package's data/ folder.)
if (!dir.exists("data")) dir.create("data")
save(data_pop_st_imp,      file = "data/data_pop_st_imp.rda",      compress = "xz")
save(pop_size_st_imp,      file = "data/pop_size_st_imp.rda",      compress = "xz")
save(data_samp_stsrs_unif, file = "data/data_samp_stsrs_unif.rda", compress = "xz")
save(data_samp_stsrs_uneq, file = "data/data_samp_stsrs_uneq.rda", compress = "xz")
