DWBS = function(data,
                N,
                d,
                numInt = 10,
                thresh = 1.3584,
                depth = "hs") {
  s_test <- 1
  e_test <- nrow(data)

  #get intervals
  intervals <- getIntervals(1:e_test, numInt) ## AR: Generate numInt subintervals in (1:e_test)

  #get test stats
  Xtilde = apply(intervals,
                 1,
                 testStat,
                 data = data,
                 depth = depth)

  #runWBS
  cp <- WBS(intervals , 2, e_test, thresh, data, depth, Xtilde)

  return(cp)
}



WBS <- function(intervals, s, e, threshold, data, depth, Xtilde) {
  if ((e - s) < 1)
    return(NULL)

  else{
    #intervals contained in s,e
    Mes <- which(apply(intervals, 1, checkIfSubInterval, super = c(s, e)))
    ## AR: Get subintervals in (s,e); returns index of intervals matrix that is in (s,e)

    Xtilde.abs <- Xtilde[Mes]

    if (length(Mes) > 1) {

      ## AR: get testStat(range = intervals[Mes, ], data, depth)
      ## which in turn calculates:
      ## The data depths of then data in the interval intervals[Mes, ]

      #weird bug where x tilde comes back as matrix
      if (!is.null(dim(Xtilde.abs))) {
        if (dim(Xtilde.abs)[2] > 2) {
          print("dim ")
          print(dim(Xtilde.abs))
        }
        Xtilde.absT <- Xtilde.abs
        Xtilde.abs <- list()
        for (i in 1:dim(Xtilde.absT)[2])
          Xtilde.abs <- append(Xtilde.abs, list(Xtilde.absT[, i]))
      }



      bs <- lapply(Xtilde.abs, which.max)
      m0 <- which.max(lapply(Xtilde.abs, max))
      b0 <- bs[[m0]] + intervals[Mes[m0], 1] - 1

      maxX <- Xtilde.abs[[m0]][bs[[m0]]]

    }


    else if (length(Mes) == 1) {
      Xtilde.abs = unlist(Xtilde.abs)
      bs <- which.max(Xtilde.abs)
      m0 <- 1
      b0 <- bs[[m0]] + intervals[Mes[m0], 1] - 1
      maxX <- max(Xtilde.abs)

    }

    else{
      return(NULL)
    }

  }

  if (maxX > threshold) { ##AR: if point b0 has CUMSUM value maxX greater than the threshold
    return(rbind(
      c(b0, maxX),
      WBS(intervals, s, b0, threshold, data,depth,Xtilde),
      WBS(intervals, b0 + 1, e, threshold, data,depth,Xtilde)
    ))
  }

  else
    return(NULL)
}






#------------------------------------------
#DWBS(rbind(data1,data2))


#packages
#library(fMultivar)
library(MASS)
#library(depth)
#library(rrcov)
#library(mrfDepth)
library(ddalpha)


##wild binary segmentation
#calculate the test statistic
#spat, hs, mahal,mahal75 are the depth parameters
testStat <- function(range, data, depth) {
  if (depth == "spat") {
    ts = testStatSpat(range, data)
  }
  else if (depth == "hs") {
    ts = testStatHs(range, data)
  }
  else if (depth == "mahal") {
    ts = testStatMahal(range, data)
  }
  else if (depth == "mahal75") {
    ts = testStatMahal75(range, data)
  }
  else{
    ts = NULL
    print("bad depth specification")
  }

  return(ts)
}

#test cusum from depth values
getStatFromDepths <- function(depths, N) { ### AR: takes depth values, returns CUMSUM
  ranks <- rank(depths, ties.method = "random")
  expected.val <- (N - 1) / 2
  std.dev <- sqrt((N ^ 2 - 1) / 12)
  cusum <- cumsum(N ^ (-0.5) * (ranks - expected.val) / std.dev)
  return(abs(cusum)[1:(length(cusum) - 1)])
}

testStatHs <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) + 1)) {
    range <- range[1]:range[2]
    N <- nrow(data[range, ])
    depths <- ddalpha::depth.halfspace(data[range, ], data[range, ])
    return(getStatFromDepths(depths, N))
  }
  else
    return(-10)
}

testStatSpat <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) + 1)) {
    range <- range[1]:range[2]
    N <- nrow(data[range, ])
    depths <- ddalpha::depth.spatial(data[range, ], data[range, ])

    return(getStatFromDepths(depths, N))
  }
  else
    return(-10)
}

testStatMahal75 <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) * 2)) {
    range <- range[1]:range[2]

    N <- nrow(data[range, ])
    depths <- ddalpha::depth.Mahalanobis(data[range, ], data[range, ], "MCD")

    return(getStatFromDepths(depths, N))
  }
  else
    return(-10)
}

testStatMahal <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) * 2)) {
    range <- range[1]:range[2]

    N <- nrow(data[range, ])
    depths <- ddalpha::depth.Mahalanobis(data[range, ], data[range, ])
    return(getStatFromDepths(depths, N))
  }
  else
    return(-10)
}


#returns indices of the intervals selected, M is the number of intervals
getIntervals <- function(indices, M) {
  ints <- t(replicate(M, sort(sample(indices, 2))))
  diffs <- (ints[, 2] - ints[, 1]) == 1
  if (any(diffs)) {
    ints[diffs, ] = getIntervals(indices, sum(diffs)) ### AR: Question for Prof. Is this so that the minimum interval length is 2?
    return(ints)
  }
  else{
    return(ints)
  }
}

#checks if an interval is a sub of another
checkIfSubInterval <- function(sub, super) { ### AR: Returns true or false
  return(sub[1] >= super[1] && sub[2] <= super[2])
}










#set.seed(440)
#test_data=rbind(replicate(2,rnorm(200)),replicate(2,rnorm(200,5)),replicate(2,rnorm(200,0.2)))



#thresh=2
#DWBS(test_data, N=nrow(test_data), d=2, numInt=100, thresh = thresh, depth = "spat")



