#Rcpp::sourceCpp('src/PELT_CPP.cpp') #Line 75 in FKWC_Changepoint
                                #Takes ~10 seconds to load
#library(mvtnorm)
library(ddalpha)




#--------------------- CONSTRUCTING FUNCTION FOR MKWP (v3) -------------------

# Additions to previous version:
# 1) Allow user to specify custom depth function

MKWP.Changepoints = function(data, depth = "Mahalanobis", beta = 10, custom_depth_function = NULL){
  if (!is.matrix(data)){
    stop("Data must be in matrix form.")
  }
  if (!is.numeric(beta)){
    stop("Parameter 'beta' must be numeric.")
  }

  if (!is.null(custom_depth_function)) {
    if (!is.function(custom_depth_function)) {
      stop("custom_depth_function must be a valid function.")
    }

    depth.values <- custom_depth_function(data)

    if (!is.numeric(depth.values) || length(depth.values) != length(data)) {
      stop("The custom depth function must return a numeric vector of the same length as the number of observations in data.")
    }
  } else {
    if (!depth %in% c("Mahalanobis", "Spatial", "Halfspace")) {
      stop("Invalid depth function. Please choose 'Mahalanobis', 'Spatial', or 'Halfspace'")
    }


    if (depth == "Mahalanobis") {
      depth.values = ddalpha::depth.Mahalanobis(x = data, data = data)
    } else if (depth == "Spatial") {
      depth.values = ddalpha::epth.spatial(x = data, data = data)
    } else if (depth == "Halfspace") {
      depth.values = ddalpha::depth.halfspace(x = data, data = data)
    }
  }
  ranks = rank(depth.values)
  cp = which(PELT(ranks,length(ranks),beta)==1)-1 #includes the changepoint 0
  if (length(cp) == 1){ #Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}

#MKWP.Changepoints(data = rbind(data1,data2))

