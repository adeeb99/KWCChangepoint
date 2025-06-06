####################### AMOC #######################

#estimate cdf of BB cdf

#ll=-1000:1000
#cdf<-function(q){sum(((-1)^ll)*exp(-2*ll^2*q^2))}
#qBB<-function(p){return(uniroot(function(q){cdf(q)-p},c(1,1.5))$root)}
#thresh=qBB(0.95)
#thresh



####################### Epidemic #######################

#estimate cdf

# gridd = 100
# bridges=replicate(10000,sde::BBridge(x=0, y=0, t0=0, T=1, N=gridd))
# dim(bridges)



maximize_bb=function(bi){
  print(bi)
  bridge_fn=Vectorize(function(k1,k2,bi=1){

    ((1-k2/gridd+k1/gridd)^(-1)+(k2/gridd-k1/gridd)^(-1))*(bridges[k2,bi]-bridges[k1,bi])^2

  },vectorize.args = c('k1','k2'))

  tst=outer(2:gridd,2:gridd,bridge_fn,bi=bi)
  tst[lower.tri(tst,T)]=-gridd
  tst2=unlist(tst)
  tst2[is.nan(tst2)]=-gridd
  # max(tst2)
  val=max(tst2)
  print(val)
  return(val)
}

# maxes=sapply(1:10000,maximize_bb)
# quantile(maxes,.95)
# 20.4

# ecdf_maxes = ecdf(maxes)




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
    stop("Invalid depth function. Please choose 'mahal' for Mahalanobis, 'mahal175' for Mahalanobis175, 'spat' for Spatial, or 'hs' for Halfspace")
  }

  return(ranks)
}

#testranks = getRanks(Multivariate_test, 'mahal')


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
      print(paste0("The estimated changepoint is ", k,
                   " with a p-value: ", 1 - cdf(Znt)))
      return(k)
    }
  }
}

# AMOC_test(Multivariate_test, setDepth = 'mahal')
# AMOC_test(ranks = testranks, useRank = T)



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
                  " with a p-value: ", 1 - ecdf_maxes(Znt)))
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
                  " with a p-value: ", 1 - ecdf_maxes(Znt)))
      return(k)
    }
  }
}

# Epidemic_test(Multivariate_test, setDepth = 'mahal')
# Epidemic_test(ranks = testranks, useRank = T)

uniMean = function(data, beta = 10){
  ranks = rank(data)
  cp = which(PELT(ranks,length(ranks),beta)==1)-1 #includes the changepoint 0
  if (length(cp) == 1){ #Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}

#uniMean_test = c(rnorm(n = 50), rnorm(n = 50, mean = 3), rnorm(n = 50, mean = 5))

#uniMean(uniMean_test)

uniScale = function(data, beta = 10){
  ranks = rank(abs(data - mean(data)))
  cp = which(PELT(ranks,length(ranks),beta)==1)-1 #includes the changepoint 0
  if (length(cp) == 1){ #Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}

#uniScale_test = c(rnorm(n = 50), rnorm(n = 50, sd = 5), rnorm(n = 50, sd = 0.1))
#uniScale(uniScale_test)
