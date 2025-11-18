# Test for an epidemic period in data

Test for a temporary change in the covariance operator of functional
data using the FKWC (functional Kruskal–Wallis covariance changepoint)
procedures outlined by Ramsay and Chenouri (2025).

## Usage

``` r
epidemic_test(data, ranks = NULL, depth = c("RPD", "FM", "LTR", "FMd", "RPDd"))
```

## Arguments

- data:

  Data in `matrix` or `data.frame` form, where each row is an
  observation and each column is a dimension.

- ranks:

  Optional if data is already ranked.

- depth:

  Depth function of choice.

## Value

A list consisting of:

- `$changepoints` : Indices of the estimated start and end points for
  the epidemic period.

- `$pvalue` : The p-value based on the null distribution.

- `$method` : A `string` `"Epidemic test (KWCChangepoint)"`

## Note

The options for the `depth` argument are as follows:

- `RPD`: Random projection depth

- `FM`: Frainman-Muniz depth

- `LTR`: \\L^2\\-root depth, most suitable for detecting changes in the
  norm

- `FMd`: Frainman-Muniz depth of the data and its first order derivative

- `RPDd`: Random projection depth of the data and its first order
  derivative

  The depth arguments that incorporate the first order derivative (which
  is approximated using
  [fda.usc::fdata.deriv](https://moviedo5.github.io/fda.usc/reference/fdata.deriv.html))
  result in a more robust detection of changes in the covariance
  structure (Ramsay and Chenouri, 2025).

## References

Ramsay, K., & Chenouri, S. (2025). Robust changepoint detection in the
variability of multivariate functional data. Journal of Nonparametric
Statistics. https://doi.org/10.1080/10485252.2025.2503891

## Examples

``` r
set.seed(11)
epi_test <- rbind(replicate(3,rnorm(200)),
                  replicate(3,rnorm(200,10)),
                  replicate(3,rnorm(200,0.2)))

epidemic_test(epi_test)
#> $changepoints
#> [1] 201 402
#> 
#> $p.value
#> [1] 0
#> 
#> $method
#> [1] "Epidemic test (KWCChangepoint)"
#> 
```
