# Find changepoints using depth-based wild binary segmentation

Detect multiple changepoints in multivariate data using the depth-based
wild binary segmentation algorithm (Ramsay and Chenouri, 2023).

## Usage

``` r
dwbs(
  data,
  numInt = 10,
  thresh = 1.3584,
  alpha = 1,
  depth = c("spat", "hs", "mahal", "mahal75")
)
```

## Arguments

- data:

  Data in `matrix` or `data.frame` form, where each row is an
  observation and each column is a dimension.

- numInt:

  Number of intervals to be generated.

- thresh:

  Numeric scalar; detection threshold. Larger values make detection more
  conservative.

- alpha:

  Set as 1 by default, applying a standard SIC penalty. Set to a number
  larger than 1 for a strengthened SIC.

- depth:

  Depth function.

## Value

A list consisting of:

- `$changepoints` : Indicies of the change-points detected; will return
  `integer(0)` if no change-points are detected.

- `$method` : A `string` `"DWBS"`

## Note

The options for the `depth` argument are as follows:

- `spat`: Spatial depth

- `hs`: Halfspace depth

- `mahal`: Mahalanobis depth

- `mahal75`: Mahalanobis depth based on re-weighted Minimum Covariance
  Determinant with 25% breakdown.

## References

Fryzlewicz, Piotr. “Wild Binary Segmentation for Multiple Change-Point
Detection.” The Annals of Statistics 42, no. 6 (2014).
https://doi.org/10.1214/14-AOS1245.

Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection of
Changepoints With a Linear Computational Cost.” Journal of the American
Statistical Association 107, no. 500 (2012): 1590–98.
https://doi.org/10.1080/01621459.2012.737745.

Ramsay, K., & Chenouri, S. (2023). Robust nonparametric multiple
changepoint detection for multivariate variability. Econometrics and
Statistics. https://doi.org/10.1016/j.ecosta.2023.09.001

## Examples

``` r
set.seed(11)
exdata <- rbind(replicate(3,rnorm(200)),
                replicate(3,rnorm(200,10)),
                replicate(3,rnorm(200,0.2)))
dwbs(data = exdata)
#> $changepoints
#> [1] 200 403
#> 
#> $method
#> [1] "DWBS"
#> 

# Increasing `numInt` will result in more accurate detection
dwbs(data = exdata, numInt = 100)
#> $changepoints
#> [1] 200 398
#> 
#> $method
#> [1] "DWBS"
#> 
```
