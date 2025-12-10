# Detect changepoints in functional data

More specifically, `fkwc()` uses the functional Kruskal-Wallis tests for
covariance changepoint algorithm (FKWC) to detect changes in the
covariance operator.

## Usage

``` r
fkwc(data, depth = c("RPD", "FM", "LTR", "FMd", "RPDd"), k = 0.25)
```

## Arguments

- data:

  Functional data in `matrix` or `data.frame` form, where each row is an
  observation/function and the columns are the grid.

- depth:

  Depth function of choice.

- k:

  Penalty constant passed to pruned exact linear time algorithm.

## Value

A list consisting of:

- `$changepoints` : Indices of the changepoints detected; will return
  `integer(0)` if no changepoints are detected.

- `$ranks` : A `vector` of depth-based ranks for each observation.

- `$method` : A `string` `"FKWC"`

## Note

The options for the `depth` argument are as follows:

- `RPD`: Random projection depth, which generally performs best

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

The penalty is of the form \$\$3.74 + k\sqrt{n}\$\$ where \\n\\ is the
number of observations. In the case that there is potentially correlated
observations, the parameter could be set to \\k=1\\. More information
could be found in the reference.

## References

Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection of
Changepoints With a Linear Computational Cost.” Journal of the American
Statistical Association 107, no. 500 (2012): 1590–98.
https://doi.org/10.1080/01621459.2012.737745.

Ramsay, K., & Chenouri, S. (2025). Robust changepoint detection in the
variability of multivariate functional data. Journal of Nonparametric
Statistics. https://doi.org/10.1080/10485252.2025.2503891

## Examples

``` r
set.seed(2)
# Generating 80 observations, with a changepoint (in our case a change in
# kernel) at observation 40
n  <- 80
k0 <- 40
T  <- 30
t  <- seq(0, 1, length.out = T)


# Both kernels K1 and K2 are Gaussian (or squared exponential) kernels but
# with different lengthscale values, and thus we hope to detect it.
K_se <- function(s, t, ell) exp(- ( (s - t)^2 ) / (2 * ell^2))
K1   <- outer(t, t, function(a,b) K_se(a,b, ell = 0.20))
K2   <- outer(t, t, function(a,b) K_se(a,b, ell = 0.07))

L1 <- chol(K1 + 1e-8 * diag(T))
L2 <- chol(K2 + 1e-8 * diag(T))

Z1 <- matrix(rnorm(k0 * T),      k0,      T)
Z2 <- matrix(rnorm((n-k0) * T),  n - k0,  T)

# We finally have an 80 x 30 matrix where the rows are the observations and
# the columns are the grid points.
X  <- rbind(Z1 %*% t(L1), Z2 %*% t(L2))

fkwc(X)
#> $changepoints
#> [1] 19 40
#> 
#> $ranks
#>  [1] 50 19 10 14 56 21 76 47 80 35 62 17 53 18 15  6  2 43 26 67 58 64 66 61 48
#> [26] 60 72 39 75 22 65  8 68  1 74 70 77 59 78 51 13 38 24 27 63 55  4 12 42 79
#> [51]  3 29 34 52 37 11  7 31 32 23 25 57 28 71 36 40 46 54 73 44 33  9  5 16 69
#> [76] 45 30 41 49 20
#> 
#> $method
#> [1] "FKWC"
#> 
```
