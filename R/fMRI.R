#' Detect changepoints in a resting state fMRI scan
#'
#' @description Functional magnetic resonance imaging scans are expected to be
#'   stationary after being pre-processed. This function attempts to find
#'   potential changepoints using the findings of Ramsay and Chenouri (2025).
#'
#'
#' @param data A four dimensional array, where the fourth dimension is time.
#' @param p Number of random vector projections, set to 100 by default.
#' @param k Penalty constant passed to pruned exact linear time algorithm.
#' @param gradient How the gradients are calculated; "exact" is precise but
#'   computationally expensive and will require parallelization.
#'
#' @returns A list consisting of:
#'  * `$changepoints` : Indices of the change-points detected; will return `integer(0)` if no changepoints are detected.
#'  * `$ranks` : A `vector` of depth-based ranks for each time stamp.
#'  * `$method` : A `string` `"fMRI changepoints (KWCChangepoint)"`
#' @export
#'
#' @note
#'
#' The gradient is, by default, calculated using a simple but imprecise method.
#' If accuracy is important, the argument "exact" will calculate the gradients
#' using [numDeriv::grad()], which will increase the run time significantly.
#'
#' The penalty is of the form \deqn{3.74 + k\sqrt{n}} where \eqn{n} is the
#' number of observations. In the case that there is potentially correlated
#' observations, the parameter could be set to \eqn{k=1}. More information could
#' be found in the reference.
#'
#' The example in this document is a simple "toy example", as good fMRI data
#' simulation requires more dependencies. For generating fMRI data, see
#' [neuRosim::simVOLfmri()], [neuRosim::simTSrestingstate()].
#'
#' @references Ramsay, K., & Chenouri, S. (2025). Robust changepoint detection
#'   in the variability of multivariate functional data. Journal of
#'   Nonparametric Statistics. https://doi.org/10.1080/10485252.2025.2503891
#'
#' @examples
#' # In order to replicate how a changepoint would appear in a resting-state
#' # fMRI scan in a manner that is not computationally expensive, this example
#' # constructs an image of a 3D ball taken at 12 time stamps. The noise, and
#' # therefore the covariance function, changes at time stamp 6.
#' x_dim <- 24
#' y_dim <- 24
#' z_dim <- 10
#' time_dim <- 12
#' image_array <- array(0, dim = c(x_dim, y_dim, z_dim, time_dim))
#'
#' center <- c(x_dim / 2, y_dim / 2, z_dim / 2)
#' radius <- min(x_dim, y_dim, z_dim) / 4
#'
#' set.seed(42)
#'
#' for (t in 1:time_dim) {
#'   for (x in 1:x_dim) {
#'     for (y in 1:y_dim) {
#'       for (z in 1:z_dim) {
#'         dist_from_center <- sqrt((x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2)
#'         if (dist_from_center <= radius) {
#'           # Adding noise with increasing variability at timestamp 6
#'           if (t <= 6) {
#'             noise <- rnorm(1, mean = 0, sd = 0.1)  # Low variability noise
#'           } else {
#'             noise <- rnorm(1, mean = 0, sd = 2)  # High variability noise
#'           }
#'           image_array[x, y, z, t] <- noise
#'         } else {
#'           # Add lower intensity noise outside the ball
#'           image_array[x, y, z, t] <- rnorm(1, mean = 0, sd = 0.005)
#'         }
#'       }
#'     }
#'   }
#' }
#' fmri_changepoints(image_array, k = 0.1, p = 10)
#'
fmri_changepoints <- function(data, p = 100, k = 0.3, gradient = c("estimate","exact")) {
  if (k < 0 || !is.numeric(k)){
    stop("Argument `k` must be a non-negative number")
  }
  gradient = match.arg(gradient)
  DIM1 <- dim(data)[1]
  DIM2 <- dim(data)[2]
  DIM3 <- dim(data)[3]
  DIM4 <- dim(data)[4]
  if (DIM1 < 3 || DIM2  < 3 || DIM3 < 3 || DIM4 < 3){
    stop("Argument `data` must be an array of size 3x3x3x3 or larger.")
  }
  dim_img <- dim(data)[1:3]
  random_directions <- array(0, dim = c(dim_img, p))

  d1 <- fda.usc::rproc2fdata(p, t = seq(0, 1, l = DIM1), norm = TRUE)
  d2 <- fda.usc::rproc2fdata(p, t = seq(0, 1, l = DIM2), norm = TRUE)
  d3 <- fda.usc::rproc2fdata(p, t = seq(0, 1, l = DIM3), norm = TRUE)

  for (i in 1:DIM1) {
    for (j in 1:DIM2) {
      for (kk in 1:DIM3) {
        random_directions[i, j, kk, ] <- d1$data[, i] * d2$data[, j] * d3$data[, kk]
      }
    }
  }
  # normalize
  for (i in 1:p) {
    random_directions[, , , i] <- random_directions[, , , i] / sqrt(sum(c(random_directions[, , , i])^2))
  }
  # project a 3d image onto unit function
  project <- function(uv, img3d) {
    num_bf <- dim(uv)[1]
    proj <- sum(uv * img3d)
    return(proj)
  }
  # compute projections of regular data
  projections <- sapply(1:p, function(rd) {
    apply(data, 4, function(x) {
      project(random_directions[, , , rd], x)
    })
  })
  projections

  # Calculate first order derivatives
  if (gradient == "estimate"){
    grads <- voxel_gradients(data = data)
  } else {
    grads <- numDeriv_gradients(data = data)
  }

  derivsx2 <- grads$dx
  derivsy2 <- grads$dy
  derivsz2 <- grads$dz

  # Projections of first order derivatives
  projections_derivative_x <- sapply(1:p, function(rd) {
    apply(derivsx2, 4, function(x) {
      project(random_directions[, , , rd], x)
    })
  })
  projections_derivative_y <- sapply(1:p, function(rd) {
    apply(derivsy2, 4, function(x) {
      project(random_directions[, , , rd], x)
    })
  })
  projections_derivative_z <- sapply(1:p, function(rd) {
    apply(derivsz2, 4, function(x) {
      project(random_directions[, , , rd], x)
    })
  })
  d_vals <- sapply(1:p, function(x) {
    ddalpha::depth.halfspace(
      cbind(projections[, x], projections_derivative_x[, x], projections_derivative_y[, x], projections_derivative_z[, x]),
      cbind(projections[, x], projections_derivative_x[, x], projections_derivative_y[, x], projections_derivative_z[, x])
    )
  })
  ranks <- rank(rowMeans(d_vals), ties.method = "random")
  cp <- which(PELT(ranks, length(ranks), beta = k * sqrt(length(ranks)) + 3.74) == 1) - 1
  list(changepoints = as.integer(cp[-c(1)]),
       ranks = ranks,
       method = "fMRI changepoints (KWCChangepoint)")
}
