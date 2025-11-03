#' Find changepoints using data-driven depth-based wild binary segmentation
#'
#' @param data A matrix or dataframe, where each row is an observation and each
#'   column is a dimension.
#' @param numInt Number of subintervals to be generated.
#' @param thresh Numeric scalar; detection threshold. Larger values make
#'   detection more conservative.
#' @param alpha Set as 1 by default, applying a standard SIC penalty. Set to a
#'   number larger than 1 for a strengthened SIC.
#' @param depth Depth function of choice, 'hs' for halfspace depth by
#'   default. Users can also choose 'mahal' for Mahalanobis, 'mahal75' for
#'   Mahalanobis MCD, 'spat' for spatial depth.
#'
#' @returns A list of changepoints.
#' @export
#'
#' @examples
#' exdata <- rbind(replicate(3,rnorm(200)),
#'                 replicate(3,rnorm(200,10)),
#'                 replicate(3,rnorm(200,0.2)))
#' dwbs(data = exdata)
#'
#' # Increasing `numInt` will result in more accurate detection
#' dwbs(data = exdata, numInt = 100)
#'
#' @references Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection
#'   of Changepoints With a Linear Computational Cost.” Journal of the American
#'   Statistical Association 107, no. 500 (2012): 1590–98.
#'   https://doi.org/10.1080/01621459.2012.737745.
#'
#'   Ramsay, K., & Chenouri, S. (2023). Robust nonparametric multiple
#'   changepoint detection for multivariate variability. Econometrics and
#'   Statistics. https://doi.org/10.1016/j.ecosta.2023.09.001
dwbs <- function(data,
                 numInt = 10,
                 thresh = 1.3584,
                 alpha = 1,
                 depth = c("hs", "spat", "mahal", "mahal75")) {
  depth = match.arg(depth)
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop("Data must be in matrix or data frame form.")
  }
  s_test <- 1
  e_test <- nrow(data)

  # get intervals
  intervals <- getIntervals(1:e_test, numInt)

  Xtilde <- apply(intervals,
    1,
    testStat,
    data = data,
    depth = depth
  )

  # runWBS
  cp <- WBS(intervals, 2, e_test, thresh, data, depth, Xtilde)

  ## Schwartz modification
  # get the possible models, select one with minimum GSC after
  # This computes the sigma sq hats
  if (!is.null(cp)) {
    cp2 <- getScwartzCriterias(cp[order(cp[, 2], decreasing = T), 1], data, depth)
  } else {
    cp2 <- list(list("cp" = NULL, "sigSq" = 1))
  }


  cp3 <- applySCH(cp2, alpha, N = nrow(data))


  return(sort(cp3))
}
