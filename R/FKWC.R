#' Detect changepoints in functional data
#'
#' More specifically, `fkwc()` uses the functional Kruskal-Wallis tests for
#' covariance changepoint algorithm (FKWC) to detect changes in the covariance
#' operator.
#'
#' @param data Functional data in `matrix` or `data.frame` form, where each row
#'   is an observation/function and the columns are the grid.
#' @param depth Depth function of choice.
#' @param k Penalty constant passed to pruned exact linear time algorithm.
#'
#' @returns A list consisting of:
#'  * `$changepoints` : Indices of the changepoints detected; will return `integer(0)` if no changepoints are detected.
#'  * `$method` : A `string` `"FKWC"`
#' @export
#'
#' @note The options for the `depth` argument are as follows:
#'
#' * `RPD`: Random projection depth, which generally performs best
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
#'  The penalty is of the form \deqn{3.74 + k\sqrt{n}} where \eqn{n} is the
#'   number of observations. In the case that there is potentially correlated
#'   observations, the parameter could be set to \eqn{k=1}. More information
#'   could be found in the reference.
#'
#'
#' @references Killick, R., P. Fearnhead, and I. A. Eckley. “Optimal Detection
#'   of Changepoints With a Linear Computational Cost.” Journal of the American
#'   Statistical Association 107, no. 500 (2012): 1590–98.
#'   https://doi.org/10.1080/01621459.2012.737745.
#'
#'   Ramsay, K., & Chenouri, S. (2025). Robust changepoint detection in the
#'   variability of multivariate functional data. Journal of Nonparametric
#'   Statistics. https://doi.org/10.1080/10485252.2025.2503891
#' @examples
#'
#' set.seed(2)
#' # Generating 80 observations, with a changepoint (in our case a change in
#' # kernel) at observation 40
#' n  <- 80
#' k0 <- 40
#' T  <- 30
#' t  <- seq(0, 1, length.out = T)
#'
#'
#' # Both kernels K1 and K2 are Gaussian (or squared exponential) kernels but
#' # with different lengthscale values, and thus we hope to detect it.
#' K_se <- function(s, t, ell) exp(- ( (s - t)^2 ) / (2 * ell^2))
#' K1   <- outer(t, t, function(a,b) K_se(a,b, ell = 0.20))
#' K2   <- outer(t, t, function(a,b) K_se(a,b, ell = 0.07))
#'
#' L1 <- chol(K1 + 1e-8 * diag(T))
#' L2 <- chol(K2 + 1e-8 * diag(T))
#'
#' Z1 <- matrix(rnorm(k0 * T),      k0,      T)
#' Z2 <- matrix(rnorm((n-k0) * T),  n - k0,  T)
#'
#' # We finally have an 80 x 30 matrix where the rows are the observations and
#' # the columns are the grid points.
#' X  <- rbind(Z1 %*% t(L1), Z2 %*% t(L2))
#'
#' fkwc(X)
#'
fkwc <- function(data,
                 depth = c("RPD", "FM", "LTR", "FMd", "RPDd"),
                 k = 0.25) {
  depth = match.arg(depth)
  if (!(k >= 0L && is.numeric(k))){
    stop("`k` should be a nonnegative number. See function documentation for more information.", call. = FALSE)
  }
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop("Data must be in matrix or data frame form.")
  }
  if(tibble::is_tibble(data)) {
    data <- as.data.frame(data)
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
  beta <- 3.74 + k * sqrt(length(ranks))
  cp <- which(PELT(ranks, length(ranks), beta) == 1) - 1
  list(changepoints = as.integer(cp[-c(1)]),
       method = "FKWC")
}


#' Multisample hypothesis test for difference in covariance operators
#'
#' Executes a multisample hypothesis test for differences in covariance
#' operators using functional Kruskal–Wallis tests for covariance (FKWC) as
#' outlined by Ramsay and Chenouri (2024). The function requires the first order
#' derivative of the functional data in order to better detect changes.
#'
#' @param data Functional data in `matrix` or `data.frame` form, where each row
#'   is an observation/function and the columns are the grid.
#' @param derivs First order derivative of the functional data in `matrix` or
#'   `data.frame` form.
#' @param g A `factor` object that indicates which sample each row of data
#'   belongs to.
#' @param p Number of random projections to be generated in order to compute
#'   random projection depths of the data.
#'
#' @returns A list consisting of:
#'  * `$statistic` : The observed test statistic.
#'  * `$pvalue` : The p-value based on the null distribution.
#'  * `$method` : A `string` `"FKWC"`
#' @export
#'
#' @references Ramsay, K., & Chenouri, S. (2024). Robust nonparametric
#'   hypothesis tests for differences in the covariance structure of functional
#'   data. Canadian Journal of Statistics, 52 (1), 43–78.
#'   https://doi.org/10.1002/cjs.11767
#'
#' @examples
#' set.seed(111)
#' t <- seq(0, 1, length.out = 200)
#'
#' ### Generating three sets of Brownian curves with different kernels, each
#' ### kernel generating 20 observations
#' # Brownian process 1
#' fd1 <- fda.usc::rproc2fdata(n = 20, t = t, sigma = "brownian",
#'                            par.list = list(scale = 10, theta = 1))
#' fd1_d <- fda.usc::fdata.deriv(fd1)
#'
#' # Brownian process 2
#' fd2 <- fda.usc::rproc2fdata(n = 20, t = t, sigma = "brownian",
#'                            par.list = list(scale = 1, theta = 1))
#' fd2_d <- fda.usc::fdata.deriv(fd2)
#'
#' # Brownian process 3
#' fd3 <- fda.usc::rproc2fdata(n = 20, t = t, sigma = "brownian",
#'                            par.list = list(scale = 1, theta = 5))
#' fd3_d <- fda.usc::fdata.deriv(fd3)
#'
#' # Functional data in one matrix and first order derivatives in another matrix
#' funcdata <- rbind(fd1$data, fd2$data, fd3$data)
#' funcderivs <- rbind(fd1_d$data, fd2_d$data, fd3_d$data)
#'
#' fkwc_multisample(data = funcdata,
#'                  derivs = funcderivs,
#'                  g = factor(rep(1:3, each = 20)),
#'                  p = 1000)
#' @seealso [fda.usc::fdata.deriv()]: for approximating the first order
#'   derivative if unavailable.
#'
#'   [fkwc_posthoc()]: for a post-hoc version of this test
fkwc_multisample <- function(data, derivs, g, p = 20) {
  if (!(is.matrix(data) || is.data.frame(data))) {
    stop("Argument `data` must be a matrix or data frame.", call. = FALSE)
  }
  if (!(is.matrix(derivs) || is.data.frame(derivs))) {
    stop("Argument `derivs` must be a matrix or data frame.", call. = FALSE)
  }
  if (!identical(dim(data), dim(derivs))) {
    stop("Arguments `data` and `derivs` must have the same dimensions.", call. = FALSE)
  }
  if (nrow(data) != length(g)) {
    stop("Argument `g` must have length equal to number of observations (rows) in `data`.", call. = FALSE)
  }
  depths <- RPDd(data = data, derivs = derivs, p = p)
  ranks <- rank(depths)
  kw <- stats::kruskal.test(ranks, g = g)
  list(statistic = as.numeric(kw$statistic),
       pvalue = kw$p.value,
       method = "FKWC multi-sample hypothesis test")
}


#' Post-hoc hypothesis test for difference in covariance operators.
#'
#' This function is post-hoc, pairwise test version of [fkwc_multisample()]
#'
#' @param data Functional data in `matrix` or `data.frame` form, where each row
#'   is an observation/function and the columns are the grid.
#' @param derivs First order derivative of the functional data in `matrix` or
#'   `data.frame` form.
#' @param g A `factor` object that indicates which sample each row of data
#'   belongs to.
#' @param p Number of random projections to be generated in order to compute
#'   random projection depths of the data.
#'
#' @returns A matrix of p-values for each pairwise comparison with a Šidák
#'   correction applied.
#' @export
#'
#' @references Ramsay, K., & Chenouri, S. (2024). Robust nonparametric
#'   hypothesis tests for differences in the covariance structure of functional
#'   data. Canadian Journal of Statistics, 52 (1), 43–78.
#'   https://doi.org/10.1002/cjs.11767
#'
#' @examples
#' set.seed(111)
#' t <- seq(0, 1, length.out = 200)
#'
#' ### Generating three sets of brownian curves with different kernels
#' # Brownian process 1
#' fd1 <- fda.usc::rproc2fdata(n = 20, t = t, sigma = "brownian",
#'                            par.list = list(scale = 10, theta = 1))
#' fd1_d <- fda.usc::fdata.deriv(fd1)
#'
#' # Brownian process 2
#' fd2 <- fda.usc::rproc2fdata(n = 20, t = t, sigma = "brownian",
#'                            par.list = list(scale = 1, theta = 1))
#' fd2_d <- fda.usc::fdata.deriv(fd2)
#'
#' # Brownian process 3
#' fd3 <- fda.usc::rproc2fdata(n = 20, t = t, sigma = "brownian",
#'                            par.list = list(scale = 1, theta = 5))
#' fd3_d <- fda.usc::fdata.deriv(fd3)
#'
#' # Functional data in one matrix and first order derivatives in another matrix
#' funcdata <- rbind(fd1$data, fd2$data, fd3$data)
#' funcderivs <- rbind(fd1_d$data, fd2_d$data, fd3_d$data)
#'
#' fkwc_posthoc(data = funcdata,
#'              derivs = funcderivs,
#'              g = factor(rep(1:3, each = 20)),
#'              p = 1000)
#' @seealso
#' [fda.usc::fdata.deriv]: for approximating the first order derivative if unavailable.
#'
#'
fkwc_posthoc <- function(data, derivs, g, p = 20) {
  if (!(is.matrix(data) || is.data.frame(data))) {
    stop("Argument `data` must be a matrix or data frame.", call. = FALSE)
  }
  if (!(is.matrix(derivs) || is.data.frame(derivs))) {
    stop("Argument `derivs` must be a matrix or data frame.", call. = FALSE)
  }
  if (!identical(dim(data), dim(derivs))) {
    stop("Arguments `data` and `derivs` must have the same dimensions.", call. = FALSE)
  }
  if (nrow(data) != length(g)) {
    stop("Argument `g` must have length equal to number of observations (rows) in `data`.", call. = FALSE)
  }
  if (!is.factor(g)) {
    stop("Argument 'g' must be a factor.", call. = FALSE)
  }
  all_pairs <- utils::combn(levels(g), m = 2)
  result_matrix <- diag(NA_real_, nrow = nlevels(g))
  rownames(result_matrix) <- levels(g)
  colnames(result_matrix) <- levels(g)
  apply(all_pairs,
    MARGIN = 2,
    FUN = function(pair) {
      index <- g %in% pair
      new_data <- data[index, ]
      new_derivs <- derivs[index, ]
      new_g <- droplevels(g[index])
      pair_results <- fkwc_multisample(
        data = new_data,
        derivs = new_derivs,
        g = new_g,
        p = p
      )
      result_matrix[pair[1], pair[2]] <<- pair_results$pvalue
      result_matrix[pair[2], pair[1]] <<- pair_results$pvalue
      invisible(NULL)
    }
  )
  m <- ncol(all_pairs)
  upper <- upper.tri(result_matrix)
  result_matrix[upper] <- pmin(1, 1 - (1 - result_matrix[upper])^m) #Sidak correction
  result_matrix[lower.tri(result_matrix)] <- t(result_matrix)[lower.tri(result_matrix)]
  return(result_matrix)
}
