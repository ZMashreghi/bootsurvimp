# bootsurvimp

**Bootstrap variance estimation for imputed survey data.**

Item nonresponse in surveys is usually handled by imputation. Treating the
imputed values as if they were observed leads to variance estimates that are
too small. `bootsurvimp` implements bootstrap methods that account for the
sampling design, the nonresponse mechanism and the imputation method, for
estimators of population totals, means and quartiles under stratified sampling.

`bootsurvimp` is the companion of the
[`bootsurv`](https://CRAN.R-project.org/package=bootsurv) package, which
provides bootstrap methods for complete survey data, and follows the same
conventions for function names, arguments and outputs.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("ZMashreghi/bootsurvimp")
```

## Main functions

| Function | Bootstrap methods | Imputation methods |
|---|---|---|
| `rescale.boot.imp()` | Independent bootstrap, rescaled-data (`"independent"`) and bootstrap-weights (`"independent.weights"`) forms (Mashreghi, Léger and Haziza, 2014); unequal response probability bootstrap (`"unequal.resp.prob"`, Mashreghi and Deng, 2023); Shao and Sitter (1996) (`"Shao.Sitter"`) | Mean, random hot-deck, regression, ratio |
| `pseudopop.boot.imp()` | Pseudo-population bootstrap: nonresponse-model (`"NRM"`), imputation-model (`"IM"`) and multiply robust (`"MR"`) schemes; pseudo-population Shao–Sitter (`"Shao.Sitter"`) | Doubly and multiply robust deterministic regression and class hot-deck; mean, ratio, hot-deck and regression for `"Shao.Sitter"` |

Both functions estimate the population total (`"total"`), mean (`"mean"`),
first quartile (`"quartile.25"`), median (`"quartile.50"` or `"median"`) or
third quartile (`"quartile.75"`).

## Example

The package ships a simulated stratified population and two samples drawn
from it (stratified simple random sampling without replacement, overall 20%
sampling fraction, strata = age groups).

```r
library(bootsurvimp)

## Independent bootstrap with regression imputation, population mean
aux <- data_samp_stsrs_unif[, c("x1", "x2", "x3")]
fit <- rescale.boot.imp(data_samp_stsrs_unif, pop_size_st_imp, R = 500,
                        stratum = "stratum", aux.data = aux,
                        parameter = "mean", imputation.method = "regression")
fit$boot.var    # bootstrap variance estimate
fit$boot.mean   # average of the bootstrap replicate estimates

## Pseudo-population bootstrap, doubly robust regression imputation
fit2 <- pseudopop.boot.imp(data_samp_stsrs_unif, pop_size_st_imp,
                           R.pop = 10, R.samp = 100, stratum = "stratum",
                           aux.IM  = list(as.matrix(aux)),
                           aux.NRM = list(as.matrix(aux[, c("x1", "x2")])),
                           parameter = "mean",
                           imputation.method = "DRDRegI",
                           bootstrap.method = "NRM")
fit2$boot.var
```

See `?rescale.boot.imp` and `?pseudopop.boot.imp` for all options and more
examples.

## Data sets

Described in `?bootsurvimp.data`.

| Data set | Description |
|---|---|
| `data_pop_st_imp` | Simulated population of 2,000 individuals in 4 age-group strata |
| `pop_size_st_imp` | Stratum population sizes (named vector) |
| `data_samp_stsrs_unif` | Stratified sample (n = 400) with uniform nonresponse within strata |
| `data_samp_stsrs_uneq` | Same sample with unequal response probabilities (`resp.prob` column) |

## Authors

- Zeinab Mashreghi (author, maintainer)
- Sophia Perron (author, student developer)
- Antonio Axalan (author, student developer)

## References

Chen, S., Haziza, D., Léger, C. and Mashreghi, Z. (2019). Pseudo-population
bootstrap methods for imputed survey data. *Biometrika*, 106(2), 369–384.
<https://doi.org/10.1093/biomet/asz001>

Chen, S., Haziza, D. and Mashreghi, Z. (2021). Multiply robust bootstrap
variance estimation in the presence of singly imputed survey data. *Journal of
Survey Statistics and Methodology*, 9(4), 810–832.
<https://doi.org/10.1093/jssam/smaa004>

Mashreghi, Z. and Deng, H. (2023). A rescaling bootstrap approach for imputed
survey data. *Journal of Survey Statistics and Methodology*, 11(1), 234–259.

Mashreghi, Z., Léger, C. and Haziza, D. (2014). Bootstrap methods for imputed
data from regression, ratio and hot-deck imputation. *The Canadian Journal of
Statistics*, 42(1), 142–167.

Shao, J. and Sitter, R. R. (1996). Bootstrap for imputed survey data.
*Journal of the American Statistical Association*, 91(435), 1278–1288.

## License

GPL (>= 3)
