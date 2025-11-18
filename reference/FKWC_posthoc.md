# Post-hoc hypothesis test for difference in covariance operators.

This function is post-hoc, pairwise test version of
[`fkwc_multisample()`](https://adeeb99.github.io/KWCChangepoint/reference/FKWC_multisample.md)

## Usage

``` r
fkwc_posthoc(data, derivs, g, p = 20)
```

## Arguments

- data:

  Functional data in `matrix` or `data.frame` form, where each row is an
  observation/function and the columns are the grid.

- derivs:

  First order derivative of the functional data in `matrix` or
  `data.frame` form.

- g:

  A `factor` object that indicates which sample each row of data belongs
  to.

- p:

  Number of random projections to be generated in order to compute
  random projection depths of the data.

## Value

A matrix of p-values for each pairwise comparison with a Šidák
correction applied.

## References

Ramsay, K., & Chenouri, S. (2024). Robust nonparametric hypothesis tests
for differences in the covariance structure of functional data. Canadian
Journal of Statistics, 52 (1), 43–78. https://doi.org/10.1002/cjs.11767

## See also

[fda.usc::fdata.deriv](https://moviedo5.github.io/fda.usc/reference/fdata.deriv.html):
for approximating the first order derivative if unavailable.

## Examples

``` r
set.seed(111)
t <- seq(0, 1, length.out = 200)

### Generating three sets of brownian curves with different kernels
# Brownian process 1
fd1 <- fda.usc::rproc2fdata(n = 20, t = t, sigma = "brownian",
                           par.list = list(scale = 10, theta = 1))
fd1_d <- fda.usc::fdata.deriv(fd1)

# Brownian process 2
fd2 <- fda.usc::rproc2fdata(n = 20, t = t, sigma = "brownian",
                           par.list = list(scale = 1, theta = 1))
fd2_d <- fda.usc::fdata.deriv(fd2)

# Brownian process 3
fd3 <- fda.usc::rproc2fdata(n = 20, t = t, sigma = "brownian",
                           par.list = list(scale = 1, theta = 5))
fd3_d <- fda.usc::fdata.deriv(fd3)

# Functional data in one matrix and first order derivatives in another matrix
funcdata <- rbind(fd1$data, fd2$data, fd3$data)
funcderivs <- rbind(fd1_d$data, fd2_d$data, fd3_d$data)

fkwc_posthoc(data = funcdata,
             derivs = funcderivs,
             g = factor(rep(1:3, each = 20)),
             p = 1000)
#>              1            2            3
#> 1           NA 3.439663e-07 2.198030e-07
#> 2 3.439663e-07           NA 9.999899e-01
#> 3 2.198030e-07 9.999899e-01           NA
```
