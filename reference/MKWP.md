# Find changepoints in multivariate data

The `mkwp()` function detects changepoints in multivariate data using
multivariate Kruskal-Wallis PELT (MKWP) algorithm developed by Ramsay
and Chenouri (2023).

## Usage

``` r
mkwp(data, depth = c("spat", "mahal", "mahal75", "hs"), k = 0.2)
```

## Arguments

- data:

  Data in `matrix` or `data.frame` form, where each row is an
  observation and each column is a dimension.

- depth:

  Depth function.

- k:

  Penalty constant passed to pruned exact linear time algorithm.

## Value

A list consisting of:

- `$changepoints` : Indices of the changepoints detected; will return
  `integer(0)` if no changepoints are detected.

- `$method` : A `string` `"Multivariate Kruskal-Wallis PELT (MKWP)"`

## Note

The options for the `depth` argument are as follows:

- `spat`: Spatial depth

- `hs`: Halfspace depth

- `mahal`: Mahalanobis depth

- `mahal75`: Mahalanobis depth based on re-weighted Minimum Covariance
  Determinant with 25% breakdown.

  Spatial depth is the default choice, as it computationally quicker
  than the other depths for larger data while giving similar result to
  other depths.

The penalty is of the form \$\$3.74 + k\sqrt{n}\$\$ where \\n\\ is the
number of observations. In the case that there is potentially correlated
observations, the parameter could be set to \\k=1\\. More information
could be found in the reference.

## References

Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection of
Changepoints With a Linear Computational Cost.” Journal of the American
Statistical Association 107, no. 500 (2012): 1590–98.
https://doi.org/10.1080/01621459.2012.737745.

Ramsay, K., & Chenouri, S. (2023). Robust nonparametric multiple
changepoint detection for multivariate variability. Econometrics and
Statistics. https://doi.org/10.1016/j.ecosta.2023.09.001

## Examples

``` r
set.seed(111)
multi_data <-rbind(replicate(3,rnorm(200)),
                   replicate(3,rnorm(200,10)),
                   replicate(3,rnorm(200,0.2)))
mkwp(multi_data)
#> $changepoints
#> [1] 201 398
#> 
#> $method
#> [1] "Multivariate Kruskal-Wallis PELT (MKWP)"
#> 
```
