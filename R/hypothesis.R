#' Conduct an AMOC (at most one changepoint) hypothesis test
#'
#' @param data A matrix or dataframe, where each row is an observation and each
#'   column is a dimension.
#' @param ranks Optional if data is already ranked.
#' @param depth Depth function of choice.
#' @note
#' The `depth` arguments are as follows:
#'
#' * `FM`: Frainman-Muniz depth
#' * `RPD`: Random projection depth
#' * `LTR`: \eqn{L^2} norm depth, most suitable for detecting changes in the norm
#' * `FMd`: Frainman-Muniz depth of the data and its first order derivative
#' * `RPDd`: Random projection depth of the data and its first order derivative
#'
#' @returns An estimated changepoint with a p-value.
#' @export
#'
#' @examples
#' test_data <- rbind(replicate(3,rnorm(200,1,1)), #before changepoint
#'                    replicate(3,rnorm(200,1,5))) #after changepoint
#'
#' amoc_test(test_data)
amoc_test <- function(data,
                      ranks = NULL,
                      depth = c("FM", "RPD", "LTR", "FMd", "RPDd")) {
  depth = match.arg(depth)
  if (!is.null(ranks)) {
    ranks <- as.numeric(ranks)
  } else {
    if (missing(data) || is.null(data)) {
      stop("Provide `data` or `ranks`.", call. = FALSE)
    }
    if (depth == "FM") {
      depths <- fda.usc::depth.FM(data)$dep
    } else if (depth == "RPD") {
      depths <- RPD(data)
    } else if (depth == "LTR") {
      depths <- fda.usc::norm.fdata(fda.usc::fdata(data))
    } else if (depth == "FMd") {
      derivs <- fda.usc::fdata.deriv(fda.usc::fdata(t(data))$data)$data
      depths <- FMp(data, t(derivs))
    } else if (depth == "RPDd") {
      derivs <- fda.usc::fdata.deriv(fda.usc::fdata(t(data))$data)$data
      depths <- RPDd(data, t(derivs))
    }
    ranks <- rank(depths)
  }
  n <- length(ranks)
  Znk <- function(kk) {
    abs((((n) * (n^2 - 1) / 12)^(-1 / 2)) * sum(ranks[1:kk] - (n + 1) / 2))
  }
  Zns <- sapply(1:(n - 1), Znk)
  k1 <- which.max(Zns)
  k <- (1:(n - 1))[k1]
  Znt <- Zns[k1]
  ll <- -1000:1000
  cdf <- function(q) {
    sum(((-1)^ll) * exp(-2 * ll^2 * q^2))
  }
  list(changepoint = as.integer(k),
       p.value = as.numeric(round(1 - cdf(Znt), 5)),
       method = "AMOC test (KWCChangepoint)")
}


#' Test for an epidemic period in data
#'
#' @param data A matrix or dataframe, where each row is an observation and each
#'   column is a dimension.
#' @param ranks Optional if data is already ranked.
#' @param depth Depth function of choice.
#' @note
#' The `depth` arguments are as follows:
#'
#' * `FM`: Frainman-Muniz depth
#' * `RPD`: Random projection depth
#' * `LTR`: \eqn{L^2} norm depth, most suitable for detecting changes in the norm
#' * `FMd`: Frainman-Muniz depth of the data and its first order derivative
#' * `RPDd`: Random projection depth of the data and its first order derivative
#'
#' @returns Estimated start and end of epidemic period with p-value.
#' @export
#'
#' @examples
#' epi_test <- rbind(replicate(3,rnorm(200)),
#'                   replicate(3,rnorm(200,10)),
#'                   replicate(3,rnorm(200,0.2)))
#'
#' epidemic_test(epi_test)
#'
epidemic_test <- function(data,
                          ranks = NULL,
                          depth = c("FM", "RPD", "LTR", "FMd", "RPDd")) {
  depth = match.arg(depth)
  if (!is.null(ranks)) {
    ranks <- as.numeric(ranks)
  } else {
    if (missing(data) || is.null(data)) {
      stop("Provide `data` or `ranks`.", call. = FALSE)
    }
    if (depth == "FM") {
      depths <- fda.usc::depth.FM(data)$dep
    } else if (depth == "RPD") {
      depths <- RPD(data)
    } else if (depth == "LTR") {
      depths <- fda.usc::norm.fdata(fda.usc::fdata(data))
    } else if (depth == "FMd") {
      derivs <- fda.usc::fdata.deriv(fda.usc::fdata(t(data))$data)$data
      depths <- FMp(data, t(derivs))
    } else if (depth == "RPDd") {
      derivs <- fda.usc::fdata.deriv(fda.usc::fdata(t(data))$data)$data
      depths <- RPDd(data, t(derivs))
    }
    ranks <- rank(depths)
  }
  n <- length(ranks)
  sign <- 12 / ((n) * (n + 1))
  mn <- 3 * (n + 1)

  rank_cumsum <- c(0, cumsum(ranks))
  total <- rank_cumsum[n + 1]

  best <- -Inf
  best_k1 <- 1L
  best_k2 <- 2L

  for (k1 in 1:(n - 1L)) {
    k2_cadidates <- (k1 + 1L):n
    len_mid  <- k2_cadidates - k1
    sum_mid  <- rank_cumsum[k2_cadidates] - rank_cumsum[k1]
    len_out  <- n - len_mid
    sum_out  <- total - sum_mid

    val_vec <-  sign * ((sum_out * sum_out) / len_out + (sum_mid * sum_mid) / len_mid) - mn

    i_local <- which.max(val_vec)
    v_local <- val_vec[i_local]
    if (v_local > best) {
      best   <- v_local
      best_k1 <- k1
      best_k2 <- k2_cadidates[i_local]
    }
  }
  list(changepoints = as.integer(c(best_k1, best_k2)),
       p.value = as.numeric(round(mean(maxes >= best), 5)),
       method = "Epidemic test (KWCChangepoint)")
}

#' Find changepoint in univariate data based on mean
#'
#' @param data A vector or one-dimensional array.
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
#' mean_test <- c(rnorm(100, mean = 0), # before change in mean
#'                rnorm(100, mean = 5)) # after change in mean
#' uni_mean(mean_test)
#'
#'
uni_mean <- function(data, beta = 10) {
  ranks <- rank(data)
  cp <- which(PELT(ranks, length(ranks), beta) == 1) - 1
  if (length(cp) == 1) {
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}


#' Find changepoints in univariate data based on scale
#'
#' @param data A vector or one-dimensional array.
#' @param beta Numeric penalty constant passed to pruned exact linear time
#'   algorithm, 10 by default.
#'
#' @returns List of changepoints.
#' @export
#'
#' @examples
#' scale_test <- c(rnorm(100, sd=5), # before change in sale
#'                 rnorm(100, sd=1)) # after change in scale
#' uni_scale(scale_test)
#'
#'
uni_scale <- function(data, beta = 10) {
  ranks <- rank(abs(data - mean(data)))
  cp <- which(PELT(ranks, length(ranks), beta) == 1) - 1
  if (length(cp) == 1) {
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}
