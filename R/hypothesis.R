
#' Rank Multivariate Data on Depth Values
#'
#' @param data Matrix of data
#' @param depth Depth function of choice
#'
#' @returns A list of ranks for each observation
#' @export
#'
#' @examples
getRanks = function(data, depth = 'spat'){
  if (depth == "spat") {
    ranks = rank(ddalpha::depth.spatial(data, data))
  }
  else if (depth == "hs") {
    ranks = rank(ddalpha::depth.halfspace(data, data))
  }
  else if (depth == "mahal") {
    ranks = rank(ddalpha::depth.Mahalanobis(data, data))
  }
  else if (depth == "mahal75") {
    ranks = rank(ddalpha::depth.Mahalanobis(data, data), "MCD")
  }
  else{
    ranks = NULL
    stop("Invalid depth function. Please choose 'mahal' for Mahalanobis, 'mahal75' for Mahalanobis MCD, 'spat' for Spatial, or 'hs' for Halfspace")
  }

  return(ranks)
}




#' Conduct an AMOC (At Most 0ne Changepoint) Hypothesis Test
#'
#' @param data Data in matrix form
#' @param ranks Optional if data is already ranked
#' @param useRank FALSE by defalut, set to TRUE to use your ranks
#' @param setDepth Depth function of choice
#' @param boundary
#'
#' @returns An estimated changepoint with a p-value
#' @export
#'
#' @examples
AMOC_test = function(data, ranks = NULL, useRank = FALSE, setDepth = 'spat', boundary = 1){
  if (useRank){
    if (is.null(ranks)){
      stop("Must provide ranks.")
    }
    else {
      n = length(ranks)
      Znk=function(kk){

        abs((((n)*(n^2-1)/12)^(-1/2))*sum(ranks[1:kk]-(n+1)/2))

      }

      Zns=sapply(boundary:(n-boundary),Znk)
      k1=which.max(Zns)
      k=(boundary:(n-boundary))[k1]
      Znt=Zns[k1]
      #return(c(Znt,k))
      ll <- -1000:1000
      cdf <- function(q){sum(((-1)^ll)*exp(-2*ll^2*q^2))}
      print(paste0("The estimated changepoint is ", k,
                   " with a p-value: ", 1 - cdf(Znt)))
      return(k)
    }
  }
  else{
    if (is.null(data)){
      stop("Must provide data.")
    }
    else {
      ranks = getRanks(data, depth = setDepth)
      n = length(ranks)
      Znk=function(kk){

        abs((((n)*(n^2-1)/12)^(-1/2))*sum(ranks[1:kk]-(n+1)/2))

      }

      Zns=sapply(boundary:(n-boundary),Znk)
      k1=which.max(Zns)
      k=(boundary:(n-boundary))[k1]
      Znt=Zns[k1]
      #return(c(Znt,k))
      ll <- -1000:1000
      cdf <- function(q){sum(((-1)^ll)*exp(-2*ll^2*q^2))}
      print(paste0("The estimated changepoint is ", k,
                   " with a p-value: ", 1 - cdf(Znt)))
      return(k)
    }
  }
}


#' Test For an Epidemic Period in Data
#'
#' @param data Matrix of data
#' @param ranks Optional if data is already ranked
#' @param useRank FALSE by defalut, set to TRUE to use your ranks
#' @param setDepth Depth function of choice
#'
#' @returns Estimated start and end of epidemic period with p-value
#' @export
#'
#' @examples
Epidemic_test = function(data, ranks = NULL, useRank = FALSE, setDepth = 'spat'){
  if (useRank){
    if (is.null(ranks)){
      stop("Must provide ranks.")
    }
    else {
      n=length(ranks)
      sign=12/((n)*(n+1))
      mn=3*(n+1)

      Wk=function(k){
        k1=k[1]
        k2=k[2]

        -sign*(cost_cpp(0,n-k2+k1-1, ranks[c(1:(k1-1),k2:n)])+cost_cpp(0,k2-k1-1, ranks[k1:(k2-1)]))-mn


      }

      pairs=apply(combn(n,2),2,sort)

      Zns=apply(pairs,2,Wk)
      ks=which.max(Zns)
      k=pairs[,ks]
      Znt=Zns[ks]

      #return(c(Znt,k)) #Test stat, then pair of changepoints - [1] 130.4072 202.0000 405.0000
      print(paste("The estimated changepoint pair is ", k[1], " and ", k[2],
                  " with a p-value: ", mean(maxes>=Znt)))
      return(k)
    }
  }
  else{
    if (is.null(data)){
      stop("Must provide data.")
    }
    else {
      ranks = getRanks(data, depth = setDepth)
      n=length(ranks)
      sign=12/((n)*(n+1))
      mn=3*(n+1)

      Wk=function(k){
        k1=k[1]
        k2=k[2]

        -sign*(cost_cpp(0,n-k2+k1-1, ranks[c(1:(k1-1),k2:n)])+cost_cpp(0,k2-k1-1, ranks[k1:(k2-1)]))-mn


      }

      pairs=apply(combn(n,2),2,sort)

      Zns=apply(pairs,2,Wk)
      ks=which.max(Zns)
      k=pairs[,ks]
      Znt=Zns[ks]

      #return(c(Znt,k)) #Test stat, then pair of changepoints - [1] 130.4072 202.0000 405.0000
      print(paste("The estimated changepoint pair is ", k[1], " and ", k[2],
                  " with a p-value: ", mean(maxes>=Znt)))
      return(k)
    }
  }
}

#' Find Changepoint in Univariate Data Based on Mean
#'
#' @param data Univariate data
#' @param beta Numeric penalty constant passed to PELT
#'
#' @returns List of changepoints
#' @export
#'
#' @examples
uniMean = function(data, beta = 10){
  ranks = rank(data)
  cp = which(PELT(ranks,length(ranks),beta)==1)-1 #includes the changepoint 0
  if (length(cp) == 1){ #Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}


#' Find Changepoints in Univariate Data Based on Scale
#'
#' @param data Univariate data
#' @param beta Numeric penalty constant passed to PELT, 10 by default
#'
#' @returns List of changepoints
#' @export
#'
#' @examples
uniScale = function(data, beta = 10){
  ranks = rank(abs(data - mean(data)))
  cp = which(PELT(ranks,length(ranks),beta)==1)-1 #includes the changepoint 0
  if (length(cp) == 1){ #Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}

