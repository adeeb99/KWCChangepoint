#' Conduct an AMOC hypothesis test
#'
#' @description Conduct an at-most one changepoint hypothesis test for changes
#'   in the covariance operator of functional data based on the FKWC (functional
#'   Kruskal–Wallis covariance changepoint) procedures outlined by Ramsay
#'   and Chenouri (2025).
#'
#'
#' @param data Data in `matrix` or `data.frame` form, where each row is an
#'   observation and each column is a dimension.
#' @param ranks Optional if data is already ranked.
#' @param depth Depth function of choice.
#' @note The options for the `depth` argument are as follows:
#'
#' * `RPD`: Random projection depth
#' * `FM`: Frainman-Muniz depth
#' * `LTR`: \eqn{L^2}-root depth, most suitable for detecting changes in the norm
#' * `FMd`: Frainman-Muniz depth of the data and its first order derivative
#' * `RPDd`: Random projection depth of the data and its first order derivative
#'
#'   The depth arguments that incorporate the first order derivative (which is
#'   approximated using [fda.usc::fdata.deriv]) result in a more robust
#'   detection of changes in the covariance structure (Ramsay and Chenouri,
#'   2025).
#'
#' @returns A list consisting of:
#'  * `$changepoint` : Index of the estimated changepoint.
#'  * `$pvalue` : The p-value based on the null distribution.
#'  * `$ranks` : A `vector` of depth-based ranks for each observation.
#'  * `$method` : A `string` `"AMOC test (KWCChangepoint)"`
#' @export
#' @references Ramsay, K., & Chenouri, S. (2025). Robust changepoint detection
#'   in the variability of multivariate functional data. Journal of
#'   Nonparametric Statistics. https://doi.org/10.1080/10485252.2025.2503891
#' @examples
#' set.seed(11)
#' test_data <- rbind(replicate(3,rnorm(200,1,1)), #before changepoint
#'                    replicate(3,rnorm(200,1,5))) #after changepoint
#'
#' amoc_test(test_data)
amoc_test <- function(data,
                      ranks = NULL,
                      depth = c("RPD", "FM", "LTR", "FMd", "RPDd")) {
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
       pvalue = as.numeric(round(1 - cdf(Znt), 5)),
       ranks = ranks,
       method = "AMOC test (KWCChangepoint)")
}


#' Test for an epidemic period in data
#'
#' @description Test for a temporary change in the covariance operator of
#'   functional data using the FKWC (functional Kruskal–Wallis covariance
#'   changepoint) procedures outlined by Ramsay and Chenouri (2025).
#'
#'
#' @param data Data in `matrix` or `data.frame` form, where each row is an
#'   observation and each column is a dimension.
#' @param ranks Optional if data is already ranked.
#' @param depth Depth function of choice.
#' @note The options for the `depth` argument are as follows:
#'
#' * `RPD`: Random projection depth
#' * `FM`: Frainman-Muniz depth
#' * `LTR`: \eqn{L^2}-root depth, most suitable for detecting changes in the norm
#' * `FMd`: Frainman-Muniz depth of the data and its first order derivative
#' * `RPDd`: Random projection depth of the data and its first order derivative
#'
#'   The depth arguments that incorporate the first order derivative (which is
#'   approximated using [fda.usc::fdata.deriv]) result in a more robust
#'   detection of changes in the covariance structure (Ramsay and Chenouri,
#'   2025).
#'
#' @returns A list consisting of:
#'  * `$changepoints` : Indices of the estimated start and end points for the epidemic period.
#'  * `$pvalue` : The p-value based on the null distribution.
#'  * `$ranks` : A `vector` of depth-based ranks for each observation.
#'  * `$method` : A `string` `"Epidemic test (KWCChangepoint)"`
#' @export
#' @references Ramsay, K., & Chenouri, S. (2025). Robust changepoint detection
#'   in the variability of multivariate functional data. Journal of
#'   Nonparametric Statistics. https://doi.org/10.1080/10485252.2025.2503891
#' @examples
#' set.seed(11)
#' epi_test <- rbind(replicate(3,rnorm(200)),
#'                   replicate(3,rnorm(200,10)),
#'                   replicate(3,rnorm(200,0.2)))
#'
#' epidemic_test(epi_test)
#'
epidemic_test <- function(data,
                          ranks = NULL,
                          depth = c("RPD", "FM", "LTR", "FMd", "RPDd")) {
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
  min_sep <- max(1L, ceiling(0.10 * n))   # at least 10% apart
  edge    <- max(0L, ceiling(0.05 * n))   # at least 5% from start/end

  k1_min <- 1L + edge
  k2_max <- n  - edge

  if (k2_max - k1_min < min_sep) {
    stop("Not enough observations to enforce the 5% edge and 10% separation rules.", call. = FALSE)
  }
  sign <- 12 / ((n) * (n + 1))
  mn <- 3 * (n + 1)

  rank_cumsum <- c(0, cumsum(ranks))
  total <- rank_cumsum[n + 1]

  best <- -Inf
  best_k1 <- 1L
  best_k2 <- 2L

  for (k1 in k1_min:(k2_max - min_sep)) {
    k2_candidates <- (k1 + min_sep):k2_max
    len_mid  <- k2_candidates - k1
    sum_mid  <- rank_cumsum[k2_candidates] - rank_cumsum[k1]
    len_out  <- n - len_mid
    sum_out  <- total - sum_mid

    val_vec <- sign * ((sum_out * sum_out) / len_out + (sum_mid * sum_mid) / len_mid) - mn

    i_local <- which.max(val_vec)
    v_local <- val_vec[i_local]
    if (v_local > best) {
      best    <- v_local
      best_k1 <- k1
      best_k2 <- k2_candidates[i_local]
    }
  }
  list(changepoints = as.integer(c(best_k1, best_k2)),
       p.value = as.numeric(round(mean(maxes >= best), 5)),
       ranks = ranks,
       method = "Epidemic test (KWCChangepoint)")
}

#' Find mean changes in a univariate sequence
#'
#' The `uni_mean()` function ranks the observations from smallest to largest,
#' then applies the pruned exact linear time algorithm with the penalty
#' parameter `beta` to detect changepoints.
#'
#' @param data A vector or one-dimensional array.
#' @param beta Numeric penalty constant passed to pruned exact linear time
#'   algorithm.
#'
#' @returns A list consisting of:
#'  * `$changepoints` : Indices of the changepoints detected; will return `integer(0)` if no changepoints are detected.
#'  * `$method` : A `string` `"Univariate Changepoint in Mean (FKWC)"`
#' @export
#'
#' @references Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection
#'   of Changepoints With a Linear Computational Cost.” Journal of the American
#'   Statistical Association 107, no. 500 (2012): 1590–98.
#'   https://doi.org/10.1080/01621459.2012.737745.
#'
#' @examples
#' set.seed(11)
#' mean_test <- c(rnorm(100, mean = 0), # before change in mean
#'                rnorm(100, mean = 5)) # after change in mean
#' uni_mean(mean_test)
#'
#'
uni_mean <- function(data, beta = 10) {
  ranks <- rank(data)
  cp <- which(PELT(ranks, length(ranks), beta) == 1) - 1
  list(changepoints = as.integer(cp[-c(1)]),
       method = "Univariate Changepoint in Mean (KWCChangepoint)")
}


#' Find scale changes in a univariate sequence
#'
#' The `uni_scale()` function ranks the observations based on their distance
#' from the mean, then applies the pruned exact linear time algorithm with the
#' penalty parameter `beta` to detect changepoints.
#'
#' @param data A vector or one-dimensional array.
#' @param beta Numeric penalty constant passed to pruned exact linear time
#'   algorithm, 10 by default.
#'
#' @returns A list consisting of:
#'  * `$changepoints` : Indices of the changepoints detected; will return `integer(0)` if no changepoints are detected.
#'  * `$method` : A `string` `"Univariate Changepoint in Scale (KWCChangepoint)"`
#' @export
#'
#' @references Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection
#'   of Changepoints With a Linear Computational Cost.” Journal of the American
#'   Statistical Association 107, no. 500 (2012): 1590–98.
#'   https://doi.org/10.1080/01621459.2012.737745.
#'
#' @examples
#' set.seed(11)
#' scale_test <- c(rnorm(100, sd=5), # before change in sale
#'                 rnorm(100, sd=1)) # after change in scale
#' uni_scale(scale_test)
#'
#'
uni_scale <- function(data, beta = 10) {
  ranks <- rank(abs(data - mean(data)))
  cp <- which(PELT(ranks, length(ranks), beta) == 1) - 1
  list(changepoints = as.integer(cp[-c(1)]),
       method = "Univariate Changepoint in Scale (KWCChangepoint)")
}
