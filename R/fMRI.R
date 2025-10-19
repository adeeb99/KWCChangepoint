# Takes a 4 dimensional input image array - 4th dimension is time
# p is the number of unit vectors - see the paper.

#' Detect changepoints in a resting state fMRI scan
#'
#' @param data A four dimensional array, where the fourth dimension is each time
#'   stamp.
#' @param p Number of random vector projections, set to 100 by default
#' @param k Part of penalty constant passed to pruned exact linear time
#'   algorithm.
#'
#' @returns A list of changepoints
#' @export
#'
#' @note
#'
#' The penalty is of the form \deqn{3.74 + k*\sqrt{n}} where \eqn{n} is the
#' number of observations. In the case that there is potentially correlated
#' observations, the parameter could be set to \eqn{k=1}. More information could
#' be found in the reference.
#'
#' @references Ramsay, K., & Chenouri, S. (2025). Robust changepoint detection
#'   in the variability of multivariate functional data. Journal of
#'   Nonparametric Statistics. https://doi.org/10.1080/10485252.2025.2503891
get_4D_changepoints <- function(data, p = 100, k = 0.3) {
  DIM1 <- dim(data)[1]
  DIM2 <- dim(data)[2]
  DIM3 <- dim(data)[3]


  dim_img <- dim(data)[1:3]
  # prod(dim_img)
  random_directions <- array(0, dim = c(dim_img, p))

  # lets try multiplying together

  d1 <- fda.usc::rproc2fdata(p, t = seq(0, 1, l = DIM1), norm = TRUE)
  d2 <- fda.usc::rproc2fdata(p, t = seq(0, 1, l = DIM2), norm = TRUE)
  d3 <- fda.usc::rproc2fdata(p, t = seq(0, 1, l = DIM3), norm = TRUE)

  # Now combine into one function
  for (i in 1:DIM1) {
    for (j in 1:DIM2) {
      for (k in 1:DIM3) {
        random_directions[i, j, k, ] <- d1$data[, i] * d2$data[, j] * d3$data[, k]
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

  # print("Done computing projections of regular data")

  # compute derivatives and their projections
  pic <- function(x, pic_num = 1) {
    # print(x)
    xx <- floor(x[1] * DIM1)
    yy <- floor(x[2] * DIM2)
    zz <- floor(x[3] * DIM3)
    if (xx < 1) {
      xx <- 1
    }
    if (zz < 1) {
      zz <- 1
    }
    if (yy < 1) {
      yy <- 1
    }

    if (xx > DIM1) {
      xx <- DIM1
    }
    if (zz > DIM3) {
      zz <- DIM3
    }
    if (yy > DIM2) {
      yy <- DIM2
    }


    return(data[xx, yy, zz, pic_num])
  }

  all_pnts <- expand.grid((1:DIM1) / DIM1, (1:DIM2) / DIM2, (1:DIM3) / DIM3)
  # print("Done all_pnts")

  all_pnts2 <- expand.grid((1:DIM1), (1:DIM2), (1:DIM3))
  # print("Done all_pnts2")

  get_derivatives <- function(pic_num) {
    deriv1 <- apply(all_pnts, 1, function(z) {
      numDeriv::grad(pic, unlist(z, use.names = FALSE),
        method = "Richardson",
        method.args = list(eps = 0.01, d = .1),
        pic_num = pic_num
      )
    })

    dxx <- dyy <- dzz <- data[, , , pic_num]

    for (x in 1:nrow(all_pnts2)) {
      y <- unlist(all_pnts2[x, ], use.names = FALSE)
      dxx[y[1], y[2], y[3]] <- deriv1[1, x]
      dyy[y[1], y[2], y[3]] <- deriv1[2, x]
      dzz[y[1], y[2], y[3]] <- deriv1[3, x]
    }

    return(list(dxx, dyy, dzz))
  }
  # print("Defined get_derviatives")

  num_scans <- dim(data)[4]
  derivs <- sapply(1:num_scans, get_derivatives)

  # print("Applied get_derivatives")


  # Compute derivatives
  derivsx2 <- data
  derivsy2 <- data
  derivsz2 <- data

  for (i in 1:(dim(data)[4])) {
    derivsx2[, , , i] <- derivs[1, i][[1]]
    derivsy2[, , , i] <- derivs[2, i][[1]]
    derivsz2[, , , i] <- derivs[3, i][[1]]
  }
  # print("Done computing derivs_2")

  projections_derivative_x <- sapply(1:p, function(rd) {
    apply(derivsx2, 4, function(x) {
      project(random_directions[, , , rd], x)
    })
  }) # projections_derivative_x
  projections_derivative_y <- sapply(1:p, function(rd) {
    apply(derivsy2, 4, function(x) {
      project(random_directions[, , , rd], x)
    })
  }) # projections_derivative_y
  projections_derivative_z <- sapply(1:p, function(rd) {
    apply(derivsz2, 4, function(x) {
      project(random_directions[, , , rd], x)
    })
  }) # projections_derivative_z


  # print("Done computing projects_derivative_")


  # for with derivatives, d_vals is depth values
  d_vals <- sapply(1:p, function(x) {
    ddalpha::depth.halfspace(
      cbind(projections[, x], projections_derivative_x[, x], projections_derivative_y[, x], projections_derivative_z[, x]),
      cbind(projections[, x], projections_derivative_x[, x], projections_derivative_y[, x], projections_derivative_z[, x])
    )
  })

  # print("Done depth values")

  # plot(d_vals,type='l')
  ranks <- rank(rowMeans(d_vals))

  # print("Done ordering ranks")

  cp <- which(PELT(ranks, length(ranks), beta = k * length(ranks) + 3.74) == 1) - 1

  if (length(cp) == 1) { # Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}
