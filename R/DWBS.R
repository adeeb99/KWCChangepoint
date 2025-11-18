#' Find changepoints using depth-based wild binary segmentation
#'
#' @description Detect multiple changepoints in multivariate data using the
#' depth-based wild binary segmentation algorithm (Ramsay and Chenouri, 2023).
#'
#'
#' @param data Data in `matrix` or `data.frame` form, where each row is an
#'   observation and each column is a dimension.
#' @param numInt Number of intervals to be generated.
#' @param thresh Numeric scalar; detection threshold. Larger values make
#'   detection more conservative.
#' @param alpha Set as 1 by default, applying a standard SIC penalty. Set to a
#'   number larger than 1 for a strengthened SIC.
#' @param depth Depth function.
#'
#' @note The options for the `depth` argument are as follows:
#'
#' * `spat`: Spatial depth
#' * `hs`: Halfspace depth
#' * `mahal`: Mahalanobis depth
#' * `mahal75`: Mahalanobis depth based on re-weighted Minimum Covariance Determinant with 25% breakdown.
#'
#' @returns A list consisting of:
#'  * `$changepoints` : Indicies of the change-points detected; will return `integer(0)` if no change-points are detected.
#'  * `$method` : A `string` `"DWBS"`
#' @export
#'
#' @examples
#' set.seed(11)
#' exdata <- rbind(replicate(3,rnorm(200)),
#'                 replicate(3,rnorm(200,10)),
#'                 replicate(3,rnorm(200,0.2)))
#' dwbs(data = exdata)
#'
#' # Increasing `numInt` will result in more accurate detection
#' dwbs(data = exdata, numInt = 100)
#'
#' @references Fryzlewicz, Piotr. “Wild Binary Segmentation for Multiple
#' Change-Point Detection.” The Annals of Statistics 42, no. 6 (2014).
#' https://doi.org/10.1214/14-AOS1245.
#'
#' Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection of
#' Changepoints With a Linear Computational Cost.” Journal of the American
#' Statistical Association 107, no. 500 (2012): 1590–98.
#' https://doi.org/10.1080/01621459.2012.737745.
#'
#' Ramsay, K., & Chenouri, S. (2023). Robust nonparametric multiple changepoint
#' detection for multivariate variability. Econometrics and Statistics.
#' https://doi.org/10.1016/j.ecosta.2023.09.001
dwbs <- function(data,
                 numInt = 10,
                 thresh = 1.3584,
                 alpha = 1,
                 depth = c("spat", "hs", "mahal", "mahal75")) {
  depth = match.arg(depth)
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop("Data must be in matrix or data frame form.")
  }
  s_test <- 1
  e_test <- nrow(data)
  lambda <- alpha * log(e_test)

  # get intervals
  intervals <- getIntervals(1:e_test, numInt)

  Xtilde <- apply(intervals, 1, testStat, data = data, depth = depth)

  # runWBS
  candidate <- WBS(intervals, 2, e_test, thresh, data, depth, Xtilde)

  depths_all <- rank(getDepths(data, depth), ties.method = "random")
  sel <- applySCH(candidate_cps = candidate, lambda = lambda, depths = depths_all)
  cp_hat <- sel$cp

  list(changepoints = as.integer(cp_hat),
       method = "DWBS")
}
