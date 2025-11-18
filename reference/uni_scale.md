# Find scale changes in a univariate sequence

The `uni_scale()` function ranks the observations based on their
distance from the mean, then applies the pruned exact linear time
algorithm with the penalty parameter `beta` to detect changepoints.

## Usage

``` r
uni_scale(data, beta = 10)
```

## Arguments

- data:

  A vector or one-dimensional array.

- beta:

  Numeric penalty constant passed to pruned exact linear time algorithm,
  10 by default.

## Value

A list consisting of:

- `$changepoints` : Indices of the changepoints detected; will return
  `integer(0)` if no changepoints are detected.

- `$method` : A `string`
  `"Univariate Changepoint in Scale (KWCChangepoint)"`

## References

Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection of
Changepoints With a Linear Computational Cost.” Journal of the American
Statistical Association 107, no. 500 (2012): 1590–98.
https://doi.org/10.1080/01621459.2012.737745.

## Examples

``` r
set.seed(11)
scale_test <- c(rnorm(100, sd=5), # before change in sale
                rnorm(100, sd=1)) # after change in scale
uni_scale(scale_test)
#> $changepoints
#> [1] 100
#> 
#> $method
#> [1] "Univariate Changepoint in Scale (KWCChangepoint)"
#> 

```
