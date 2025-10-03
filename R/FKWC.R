#' @keywords internal
RPD=function(data){
  gen_direction=replicate(20,rnorm(100))
  gen_direction=apply(gen_direction,2,function(x){x/sqrt(sum(x^2))})

  #data is n by m
  # dim(data)
  bd=apply(gen_direction,2,function(u){
    bdd=data%*%u
    rnk=rank(bdd)
    rnk*(length(bdd)-rnk)})
  return(rowMeans(bd))
}

#' @keywords internal
FMp=function(data,derivs){
  dp=sapply(1:100,function(x){ddalpha::depth.halfspace(cbind( data[,x],derivs[,x]),cbind( data[,x],derivs[,x]),num.directions =100)   })
  return(rowMeans(dp))
}



#' Find Changepoints Using Functional Kruskall-Wallis Tests for Covariance
#' Algorithm
#'
#' @param funcdata Functional data in fdata form, where each row is an
#'   observation and each column is a dimension.
#' @param depth Depth function of choice
#' @param k Part of penalty constant passed to PELT
#'
#' @returns A list of changepoints
#' @export
#'
#' @examples
FKWC <- function(funcdata, depth = "FM_depth", k = 0.25){

  if (!fda.usc::is.fdata(funcdata)){
    stop("Data must be in fdata form.")
  }
  if (!depth %in% c("FM_depth","RPD_depth","LTR_depth","FM_depth_d")) {
    stop("Invalid depth function. Please choose 'FM_depth', 'RPD_depth', 'LTR_depth', or
         'FM_depth_d'")
  }

  rownames(funcdata$data) = as.character(1:(nrow(funcdata$data)))

  if (depth == "FM_depth"){
    depths = fda.usc::depth.FM(funcdata)$dep
  } else if (depth == "RPD_depth"){
    depths = RPD(funcdata$data)
  } else if (depth == "LTR_depth"){
    depths = c(fda.usc::norm.fdata(funcdata))
  } else if (depth == "FM_depth_d"){
    derivs = MFHD::derivcurves(funcdata$data)
    depths = FMp(funcdata$data,derivs)
  } else if (depth == "RPD_depth_d"){
    derivs = MFHD::derivcurves(funcdata$data)
    depths = RPDd(funcdata$data,derivs)
  }
  ranks = rank(depths)
  beta = 3.74 + k*sqrt(length(ranks))
  cp = which(PELT(ranks,length(ranks),beta)==1)-1 #includes the changepoint 0

  if (length(cp) == 1){ #Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}

#' @keywords internal
RPDd <- function(data, derivs, p = 20,
                  depth_fun = ddalpha::depth.simplicial) {
  X <- as.matrix(data)
  D <- as.matrix(derivs)
  n <- nrow(X)
  m <- ncol(X)
  # Generate p random directions
  U <- replicate(p, rnorm(m))
  U <- apply(U, 2, function(v) v / sqrt(sum(v^2)))
  # For each direction u_k, compute 2D projections (x·u_k, x'·u_k),
  # then take depth of the n points w.r.t. themselves
  bd <- apply(U, 2, function(u) {
    bdd <- cbind(as.vector(X %*% u), as.vector(D %*% u))
    depth_fun(x = bdd, data = bdd)
  })

  drop(rowMeans(bd))
}


#' Multisample hypothesis test for differences in covariance operators using
#' Functional Kruskal–Wallis Tests for Covariance.
#'
#' @param data Functional data in fdata form, where each row is an observation
#'   and each column is a dimension.
#' @param derivs First order derivative of the functional data in fdata form.
#' @param g Factor object that indicates which sample each row of data belongs
#'   to.
#' @param p Number of random projections to be generated in order to compute
#'   random projection depths of the data.
#'
#' @returns A list containing the test statistic and the p-value.
#' @export
#'
#' @examples
FKWC_multisample <- function(data, derivs, g, p = 20){
  depths <- RPDd2(data = data, derivs = derivs, p = p)
  ranks <- rank(depths)
  kw <- kruskal.test(ranks, g = g)
  list(statistic = as.numeric(kw$statistic), p.value = kw$p.value)
}


#' Pairwise comparison post-hoc hypothesis test for differences in covariance
#' operators using Functional Kruskal–Wallis Tests for Covariance.
#'
#' @param data Functional data in fdata form, where each row is an observation
#'   and each column is a dimension.
#' @param derivs First order derivative of the functional data in fdata form.
#' @param g Factor object that indicates which sample each row of data belongs
#'   to.
#' @param p Number of random projections to be generated in order to compute
#'   random projection depths of the data.
#'
#' @returns A matrix of p-values for each pairwise comparison with a Šidák
#'   correction applied.
#' @export
#'
#' @examples
FKWC_posthoc <- function(data, derivs, g, p = 20){
  all_pairs <- combn(levels(g), m = 2) # Get all possible pairs
  result_matrix <- diag(NA_real_, nrow = nlevels(g))
  rownames(result_matrix) <- levels(g)
  colnames(result_matrix) <- levels(g)
  apply(all_pairs, MARGIN = 2,
        FUN = function(pair){
          index <- g %in% pair
          new_data <- data[index,]
          new_derivs <- derivs[index,]
          new_g <- droplevels(g[index])
          pair_results <- FKWC_multisample(data = new_data,
                                           derivs = new_derivs,
                                           g = new_g,
                                           p = p)
          result_matrix[pair[1], pair[2]] <<- pair_results$p.value
          result_matrix[pair[2], pair[1]] <<- pair_results$p.value
          invisible(NULL)
        })
  m  <- ncol(all_pairs)
  upper <- upper.tri(result_matrix)
  result_matrix[upper] <- pmin(1, 1 - (1 - result_matrix[upper])^m)
  result_matrix[lower.tri(result_matrix)] <- t(result_matrix)[lower.tri(result_matrix)]
  print(result_matrix)
}
