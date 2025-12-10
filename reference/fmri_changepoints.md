# Detect changepoints in a resting state fMRI scan

Functional magnetic resonance imaging scans are expected to be
stationary after being pre-processed. This function attempts to find
potential changepoints using the findings of Ramsay and Chenouri (2025).

## Usage

``` r
fmri_changepoints(data, p = 100, k = 0.3)
```

## Arguments

- data:

  A four dimensional array, where the fourth dimension is time.

- p:

  Number of random vector projections, set to 100 by default.

- k:

  Penalty constant passed to pruned exact linear time algorithm.

## Value

A list consisting of:

- `$changepoints` : Indices of the change-points detected; will return
  `integer(0)` if no changepoints are detected.

- `$ranks` : A `vector` of depth-based ranks for each time stamp.

- `$method` : A `string` `"fMRI changepoints (KWCChangepoint)"`

## Note

The penalty is of the form \$\$3.74 + k\sqrt{n}\$\$ where \\n\\ is the
number of observations. In the case that there is potentially correlated
observations, the parameter could be set to \\k=1\\. More information
could be found in the reference.

The example in this document is a simple "toy example", as good fMRI
data simulation requires more dependencies. For generating fMRI data,
see `neuRosim::simVOLfmri()`, `neuRosim::simTSrestingstate()`.

## References

Ramsay, K., & Chenouri, S. (2025). Robust changepoint detection in the
variability of multivariate functional data. Journal of Nonparametric
Statistics. https://doi.org/10.1080/10485252.2025.2503891

## Examples

``` r
# In order to replicate how a changepoint would appear in a resting-state
# fMRI scan in a manner that is not computationally expensive, this example
# constructs an image of a 3D ball taken at 12 time stamps. The noise, and
# therefore the covariance function, changes at time stamp 6.
x_dim <- 24
y_dim <- 24
z_dim <- 10
time_dim <- 12
image_array <- array(0, dim = c(x_dim, y_dim, z_dim, time_dim))

center <- c(x_dim / 2, y_dim / 2, z_dim / 2)
radius <- min(x_dim, y_dim, z_dim) / 4

set.seed(42)

for (t in 1:time_dim) {
  for (x in 1:x_dim) {
    for (y in 1:y_dim) {
      for (z in 1:z_dim) {
        dist_from_center <- sqrt((x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2)
        if (dist_from_center <= radius) {
          # Adding noise with increasing variability at timestamp 6
          if (t <= 6) {
            noise <- rnorm(1, mean = 0, sd = 0.1)  # Low variability noise
          } else {
            noise <- rnorm(1, mean = 0, sd = 2)  # High variability noise
          }
          image_array[x, y, z, t] <- noise
        } else {
          # Add lower intensity noise outside the ball
          image_array[x, y, z, t] <- rnorm(1, mean = 0, sd = 0.005)
        }
      }
    }
  }
}
fmri_changepoints(image_array, k = 0.1, p = 10)
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> $changepoints
#> [1] 6
#> 
#> $ranks
#>  [1]  7 12  8  9 10 11  4  3  5  6  2  1
#> 
#> $method
#> [1] "fMRI changepoints (KWCChangepoint)"
#> 
```
