## HELPERS FOR DWBS ##

# schwartz criteria for choosing the threshold
# returns all the different models which are the selected change-points and the sigmahatsquared value
# (see paper, section on choosing threshold)
#' @noRd
getScwartzCriterias <- function(cp, data, depth) {
  # get indidvidual criteria for a set of cp
  sqhat <- getsqhat(NULL, data, depth)

  abc <- list("cp" = NULL, "sigSq" = sqhat)

  models <- list(abc)

  for (i in 1:length(cp)) {
    sqhat <- getsqhat(cp[1:i], data, depth)

    abc$cp <- cp[1:i]
    abc$sigSq <- sqhat

    models <- append(models, list(abc))
  }

  return(models)
}

#' @keywords internal
#' @noRd
getsqhat <- function(cp, data, depth) {
  depths <- rank(getDepths(data, depth), ties.method = "random")
  # depths<-getDepths(data,depth)

  # if at least 1 cp
  if (length(cp) >= 1) {
    N <- length(depths)
    indicies <- cbind(c(1, sort(cp)), c(sort(cp), N + 1))

    getGroups <- function(vec, dat) {
      return(dat[vec[1]:(vec[2] - 1)])
    }

    breaks <- apply(indicies, 1, getGroups, dat = depths)

    absSum <- lapply(breaks, function(x) {
      sum((x - mean(x))^2)
    })
    sighatSq <- sum(unlist(absSum)) / N


    return(sighatSq)
  } else {
    N <- length(depths)
    sighatSq <- mean((depths - mean(depths))^2)
    return(sighatSq)
  }
}

# get indidvidual criteria for a set of cp
#' @keywords internal
#' @noRd
getsSic2 <- function(cp, sighatSq, alpha, N) {
  # if at least 1 cp
  if (length(cp) >= 1) {
    sSic <- (N / 2) * log(sighatSq) + length(cp) * (log(N))^alpha
  } else {
    sSic <- (N / 2) * log(sighatSq)
  }

  return(sSic)
}


#' @keywords internal
#' @noRd
applySCH <- function(models, alpha, N) {
  # get SIC for all amounts of cp
  sSic <- getsSic2(models[[1]]$cp, models[[1]]$sigSq, alpha, N)


  if (length(models) > 1) {
    for (j in 2:length(models)) {
      SIC_j <- getsSic2(models[[j]]$cp, models[[j]]$sigSq, alpha, N)
      sSic <- c(sSic, SIC_j)
    }
  }
  minVal <- which.min(sSic)

  return(models[[minVal]]$cp)
}
#' @keywords internal
#' @noRd
WBS <- function(intervals, s, e, threshold, data, depth, Xtilde) {
  if ((e - s) < 1) {
    return(NULL)
  } else {
    # intervals contained in s,e
    Mes <- which(apply(intervals, 1, checkIfSubInterval, super = c(s, e)))
    ## Get subintervals in (s,e); returns index of intervals matrix that is in (s,e)

    Xtilde.abs <- Xtilde[Mes]

    if (length(Mes) > 1) {
      ## get testStat(range = intervals[Mes, ], data, depth)
      ## which in turn calculates:
      ## The data depths of then data in the interval intervals[Mes, ]

      if (!is.null(dim(Xtilde.abs))) {
        Xtilde.absT <- Xtilde.abs
        Xtilde.abs <- list()
        for (i in 1:dim(Xtilde.absT)[2]) {
          Xtilde.abs <- append(Xtilde.abs, list(Xtilde.absT[, i]))
        }
      }


      bs <- lapply(Xtilde.abs, which.max)
      m0 <- which.max(lapply(Xtilde.abs, max))
      b0 <- bs[[m0]] + intervals[Mes[m0], 1] - 1

      maxX <- Xtilde.abs[[m0]][bs[[m0]]]
    } else if (length(Mes) == 1) {
      Xtilde.abs <- unlist(Xtilde.abs)
      bs <- which.max(Xtilde.abs)
      m0 <- 1
      b0 <- bs[[m0]] + intervals[Mes[m0], 1] - 1
      maxX <- max(Xtilde.abs)
    } else {
      return(NULL)
    }
  }

  if (maxX > threshold) { ## if point b0 has CUMSUM value maxX greater than the threshold
    return(rbind(
      c(b0, maxX),
      WBS(intervals, s, b0, threshold, data, depth, Xtilde),
      WBS(intervals, b0 + 1, e, threshold, data, depth, Xtilde)
    ))
  } else {
    return(NULL)
  }
}
#' @keywords internal
#' @noRd
testStat <- function(range, data, depth) {
  if (depth == "spat") {
    ts <- testStatSpat(range, data)
  } else if (depth == "hs") {
    ts <- testStatHs(range, data)
  } else if (depth == "mahal") {
    ts <- testStatMahal(range, data)
  } else if (depth == "mahal75") {
    ts <- testStatMahal75(range, data)
  } else {
    ts <- NULL
    stop("Invalid depth function. Please choose 'mahal' for Mahalanobis, 'mahal75' for Mahalanobis MCD, 'spat' for Spatial, or 'hs' for Halfspace")
  }
  return(ts)
}

# test cusum from depth values
#' @noRd
getStatFromDepths <- function(depths, N) {
  ranks <- rank(depths, ties.method = "random")
  expected.val <- (N - 1) / 2
  std.dev <- sqrt((N^2 - 1) / 12)
  cusum <- cumsum(N^(-0.5) * (ranks - expected.val) / std.dev)
  return(abs(cusum)[1:(length(cusum) - 1)])
}

#' @keywords internal
#' @noRd
testStatHs <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) + 1)) {
    range <- range[1]:range[2]
    N <- nrow(data[range, ])
    depths <- ddalpha::depth.halfspace(data[range, ], data[range, ])
    return(getStatFromDepths(depths, N))
  } else {
    return(0)
  }
}

#' @keywords internal
#' @noRd
testStatSpat <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) + 1)) {
    range <- range[1]:range[2]
    N <- nrow(data[range, ])
    depths <- ddalpha::depth.spatial(data[range, ], data[range, ])

    return(getStatFromDepths(depths, N))
  } else {
    return(0)
  }
}

#' @keywords internal
#' @noRd
testStatMahal75 <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) * 2)) {
    range <- range[1]:range[2]

    N <- nrow(data[range, ])
    depths <- ddalpha::depth.Mahalanobis(data[range, ], data[range, ], "MCD")

    return(getStatFromDepths(depths, N))
  } else {
    return(0)
  }
}

#' @keywords internal
#' @noRd
testStatMahal <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) * 2)) {
    range <- range[1]:range[2]

    N <- nrow(data[range, ])
    depths <- ddalpha::depth.Mahalanobis(data[range, ], data[range, ])
    return(getStatFromDepths(depths, N))
  } else {
    return(0)
  }
}


# returns indices of the intervals selected, M is the number of intervals
#' @keywords internal
#' @noRd
getIntervals <- function(indices, M) {
  ints <- t(replicate(M, sort(sample(indices, 2))))
  diffs <- (ints[, 2] - ints[, 1]) == 1
  if (any(diffs)) {
    ints[diffs, ] <- getIntervals(indices, sum(diffs))
    return(ints)
  } else {
    return(ints)
  }
}

# checks if an interval is a sub of another
#' @keywords internal
#' @noRd
checkIfSubInterval <- function(sub, super) {
  return(sub[1] >= super[1] && sub[2] <= super[2])
}



#' @noRd
getDepths <- function(data, depth) {
  if (depth == "spat") {
    ts <- ddalpha::depth.spatial(data, data)
  } else if (depth == "hs") {
    ts <- ddalpha::depth.halfspace(data, data)
  } else if (depth == "mahal") {
    ts <- ddalpha::depth.Mahalanobis(data, data)
  } else if (depth == "mahal75") {
    ts <- ddalpha::depth.Mahalanobis(data, data, "MCD")
  } else {
    ts <- NULL
    stop("Invalid depth function. Please choose 'mahal' for Mahalanobis, 'mahal75' for Mahalanobis MCD, 'spat' for Spatial, or 'hs' for Halfspace")
  }

  return(ts)
}

#' @keywords internal
#' @noRd
RPD <- function(data, p = 20,
                depth_fun = ddalpha::depth.simplicial) {
  X <- as.matrix(data)
  n <- nrow(X)
  m <- ncol(X)
  # Generate p random directions
  U <- replicate(p, stats::rnorm(m))
  U <- apply(U, 2, function(v) v / sqrt(sum(v^2)))
  bd <- apply(U, 2, function(u) {
    bdd <- X %*% u
    rnk <- rank(bdd)
    rnk * (length(bdd) - rnk)
  })
  return(rowMeans(bd))
}

#' @keywords internal
#' @noRd
FMp <- function(data, derivs) {
  dp <- sapply(1:100, function(x) {
    ddalpha::depth.halfspace(cbind(data[, x], derivs[, x]), cbind(data[, x], derivs[, x]), num.directions = 100)
  })
  return(rowMeans(dp))
}


#' @keywords internal
#' @noRd
RPDd <- function(data, derivs, p = 20,
                 depth_fun = ddalpha::depth.simplicial) {
  X <- as.matrix(data)
  D <- as.matrix(derivs)
  n <- nrow(X)
  m <- ncol(X)
  # Generate p random directions
  U <- replicate(p, stats::rnorm(m))
  U <- apply(U, 2, function(v) v / sqrt(sum(v^2)))
  # For each direction u_k, compute 2D projections (x·u_k, x'·u_k),
  # then take depth of the n points w.r.t. themselves
  bd <- apply(U, 2, function(u) {
    bdd <- cbind(as.vector(X %*% u), as.vector(D %*% u))
    depth_fun(x = bdd, data = bdd)
  })

  drop(rowMeans(bd))
}
#' Rank multivariate data according to depth values
#'
#' @param data A matrix or dataframe, where each row is an observation and each
#'   column is a dimension.
#' @param depth Depth function of choice. It is 'spat' for spatial depth by
#'   default. User can also choose 'mahal' for Mahalanobis, 'mahal75' for
#'   Mahalanobis MCD, or 'hs' for halfspace depth.
#'
#' @returns A list of ranks for each observation.
#' @noRd
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
