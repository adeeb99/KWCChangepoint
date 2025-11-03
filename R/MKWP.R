#' Find changepoints in multivariate data using multivariate Kruskal-Wallis
#' PELT.
#'
#' @param data A matrix or dataframe, where each row is an observation and each
#'   column is a dimension.
#' @param depth Depth function of choice. It is 'mahal' for Mahalanobis depth by
#'   default. User can also choose 'mahal75' for Mahalanobis MCD, 'hs' for
#'   halfspace, or 'spat' for spatial depth.
#' @param k Penalty constant passed to pruned exact linear time algorithm.
#'
#' @returns A list of changepoints.
#' @export
#'
#' @references Ramsay, K., & Chenouri, S. (2023). Robust nonparametric multiple
#'   changepoint detection for multivariate variability. Econometrics and
#'   Statistics. https://doi.org/10.1016/j.ecosta.2023.09.001
#' @examples
#' set.seed(111)
#' multi_data <-rbind(replicate(3,rnorm(200)),
#'                    replicate(3,rnorm(200,10)),
#'                    replicate(3,rnorm(200,0.2)))
#' mkwp(multi_data)
#'
mkwp <- function(data,
                 depth = c("mahal", "mahal75", "spat", "hs"),
                 k = 0.2) {
  depth = match.arg(depth)
  if (!(is.matrix(data) || is.data.frame(data))) {
    stop("`data` must be a matrix or data frame.", call. = FALSE)
  }
  if (is.data.frame(data) && !all(vapply(data, is.numeric, logical(1)))) {
    stop("All columns of `data` must be numeric.", call. = FALSE)
  }
  depth.values <- getDepths(data = data, depth = depth)
  ranks <- rank(depth.values)
  beta <- 3.74 + k * sqrt(length(ranks))
  cp <- which(PELT(ranks, length(ranks), beta) == 1) - 1 # includes the changepoint 0
  if (length(cp) == 1) { # Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}
