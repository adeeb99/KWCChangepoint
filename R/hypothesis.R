#' Rank multivariate data according to depth values
#'
#' @param data A matrix or dataframe, where each row is an observation and each
#'   column is a dimension.
#' @param depth Depth function of choice. It is 'spat' for spatial depth by
#'   default. User can also choose 'mahal' for Mahalanobis, 'mahal75' for
#'   Mahalanobis MCD, or 'hs' for halfspace depth.
#'
#' @returns A list of ranks for each observation.
#' @export
#'
#' @examples
getRanks <- function(data, depth = "spat") {
  if (depth == "spat") {
    ranks <- rank(ddalpha::depth.spatial(data, data))
  } else if (depth == "hs") {
    ranks <- rank(ddalpha::depth.halfspace(data, data))
  } else if (depth == "mahal") {
    ranks <- rank(ddalpha::depth.Mahalanobis(data, data))
  } else if (depth == "mahal75") {
    ranks <- rank(ddalpha::depth.Mahalanobis(data, data), "MCD")
  } else {
    ranks <- NULL
    stop("Invalid depth function. Please choose 'mahal' for Mahalanobis, 'mahal75' for Mahalanobis MCD, 'spat' for Spatial, or 'hs' for Halfspace")
  }

  return(ranks)
}


#' Conduct an AMOC (at most one changepoint) hypothesis test
#'
#' @param data A matrix or dataframe, where each row is an observation and each
#'   column is a dimension.
#' @param ranks Optional if data is already ranked
#' @param useRank FALSE by defalut, set to TRUE to use your ranks
#' @param depth Depth function of choice. It is 'spat' for spatial depth by
#'   default. User can also choose 'mahal' for Mahalanobis, 'mahal75' for
#'   Mahalanobis MCD, or 'hs' for halfspace depth.
#' @param boundary
#'
#' @returns An estimated changepoint with a p-value.
#' @export
#'
#' @examples
AMOC_test <- function(data, ranks = NULL, useRank = FALSE, depth = "spat", boundary = 1) {
  if (useRank) {
    if (is.null(ranks)) {
      stop("Must provide ranks.")
    } else {
      n <- length(ranks)
      Znk <- function(kk) {
        abs((((n) * (n^2 - 1) / 12)^(-1 / 2)) * sum(ranks[1:kk] - (n + 1) / 2))
      }

      Zns <- sapply(boundary:(n - boundary), Znk)
      k1 <- which.max(Zns)
      k <- (boundary:(n - boundary))[k1]
      Znt <- Zns[k1]
      # return(c(Znt,k))
      ll <- -1000:1000
      cdf <- function(q) {
        sum(((-1)^ll) * exp(-2 * ll^2 * q^2))
      }
      print(paste0(
        "Estimated changepoint is ", k,
        " with a p-value: ", round(1 - cdf(Znt), 5)
      ))
      return(k)
    }
  } else {
    if (is.null(data)) {
      stop("Must provide data.")
    } else {
      ranks <- getRanks(data, depth = depth)
      n <- length(ranks)
      Znk <- function(kk) {
        abs((((n) * (n^2 - 1) / 12)^(-1 / 2)) * sum(ranks[1:kk] - (n + 1) / 2))
      }

      Zns <- sapply(boundary:(n - boundary), Znk)
      k1 <- which.max(Zns)
      k <- (boundary:(n - boundary))[k1]
      Znt <- Zns[k1]
      # return(c(Znt,k))
      ll <- -1000:1000
      cdf <- function(q) {
        sum(((-1)^ll) * exp(-2 * ll^2 * q^2))
      }
      print(paste0(
        "Estimated changepoint is ", k,
        " with a p-value: ", round(1 - cdf(Znt), 5)
      ))
      return(k)
    }
  }
}


#' Test for an epidemic period in data
#'
#' @param data A matrix or dataframe, where each row is an observation and each
#'   column is a dimension.
#' @param ranks Optional if data is already ranked.
#' @param useRank FALSE by defalut, set to TRUE to use your ranks
#' @param depth Depth function of choice. It is 'spat' for spatial depth by
#'   default. User can also choose 'mahal' for Mahalanobis, 'mahal75' for
#'   Mahalanobis MCD, or 'hs' for halfspace depth.
#'
#' @returns Estimated start and end of epidemic period with p-value.
#' @export
#'
#' @examples
Epidemic_test <- function(data, ranks = NULL, useRank = FALSE, depth = "spat") {
  if (useRank) {
    if (is.null(ranks)) {
      stop("Must provide ranks.")
    } else {
      n <- length(ranks)
      sign <- 12 / ((n) * (n + 1))
      mn <- 3 * (n + 1)

      Wk <- function(k) {
        k1 <- k[1]
        k2 <- k[2]

        -sign * (cost_cpp(0, n - k2 + k1 - 1, ranks[c(1:(k1 - 1), k2:n)]) + cost_cpp(0, k2 - k1 - 1, ranks[k1:(k2 - 1)])) - mn
      }

      pairs <- apply(combn(n, 2), 2, sort)

      Zns <- apply(pairs, 2, Wk)
      ks <- which.max(Zns)
      k <- pairs[, ks]
      Znt <- Zns[ks]

      # return(c(Znt,k)) #Test stat, then pair of changepoints - [1] 130.4072 202.0000 405.0000
      print(paste(
        "The estimated changepoint pair is ", k[1], " and ", k[2],
        " with a p-value: ", mean(maxes >= Znt)
      ))
      return(k)
    }
  } else {
    if (is.null(data)) {
      stop("Must provide data.")
    } else {
      ranks <- getRanks(data, depth = depth)
      n <- length(ranks)
      sign <- 12 / ((n) * (n + 1))
      mn <- 3 * (n + 1)

      Wk <- function(k) {
        k1 <- k[1]
        k2 <- k[2]

        -sign * (cost_cpp(0, n - k2 + k1 - 1, ranks[c(1:(k1 - 1), k2:n)]) + cost_cpp(0, k2 - k1 - 1, ranks[k1:(k2 - 1)])) - mn
      }

      pairs <- apply(combn(n, 2), 2, sort)

      Zns <- apply(pairs, 2, Wk)
      ks <- which.max(Zns)
      k <- pairs[, ks]
      Znt <- Zns[ks]

      # return(c(Znt,k)) #Test stat, then pair of changepoints - [1] 130.4072 202.0000 405.0000
      print(paste(
        "Estimated changepoint pair is ", k[1], " and ", k[2],
        " with a p-value: ", round(mean(maxes >= Znt), 5)
      ))
      return(k)
    }
  }
}

#' Find changepoint in univariate data based on mean
#'
#' @param data Univariate data.
#' @param beta Numeric penalty constant passed to pruned exact linear time
#'   algorithm.
#'
#' @returns List of changepoints.
#' @export
#'
#' @references Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection
#'   of Changepoints With a Linear Computational Cost.” Journal of the American
#'   Statistical Association 107, no. 500 (2012): 1590–98.
#'   https://doi.org/10.1080/01621459.2012.737745.
#'
#' @examples
uniMean <- function(data, beta = 10) {
  ranks <- rank(data)
  cp <- which(PELT(ranks, length(ranks), beta) == 1) - 1 # includes the changepoint 0
  if (length(cp) == 1) { # Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}


#' Find changepoints in univariate data based on scale
#'
#' @param data Univariate data.
#' @param beta Numeric penalty constant passed to pruned exact linear time
#'   algorithm, 10 by default.
#'
#' @returns List of changepoints.
#' @export
#'
#' @examples
uniScale <- function(data, beta = 10) {
  ranks <- rank(abs(data - mean(data)))
  cp <- which(PELT(ranks, length(ranks), beta) == 1) - 1 # includes the changepoint 0
  if (length(cp) == 1) { # Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}
