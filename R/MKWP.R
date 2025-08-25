
#' Find Changepoints Using Multivariate Kruskal-Wallis PELT
#'
#' @param data Matrix of data
#' @param depth Depth function of choice
#' @param beta Numeric penalty constant passed to PELT
#' @param custom_depth_function Use your own custom depth function (NULL by default)
#'
#' @returns A list of changepoints
#' @export
#'
#' @examples
MKWP = function(data, depth = "mahal", beta = 10, custom_depth_function = NULL){
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
    if (!depth %in% c("mahal", "mahal75", "spat", "hs")) {
      stop("Invalid depth function. Please choose 'mahal' for Mahalanobis, 'mahal75' for Mahalanobis MCD, 'spat' for Spatial, or 'hs' for Halfspace")
    }


    if (depth == "mahal") {
      depth.values = ddalpha::depth.Mahalanobis(x = data, data = data)
    } else if (depth == "mahal75") {
      depth.values = ddalpha::depth.spatial(x = data, data = data, "MCD")
    } else if (depth == "spat") {
      depth.values = ddalpha::depth.spatial(x = data, data = data)
    } else if (depth == "hs") {
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



