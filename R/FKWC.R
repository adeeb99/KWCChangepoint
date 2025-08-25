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

#' @keywords internal
RPDd=function(data,derivs){
  gen_direction=replicate(20,rnorm(100))
  gen_direction=apply(gen_direction,2,function(x){x/sqrt(sum(x^2))})

  #data is n by m
  # dim(data)
  bd=apply(gen_direction,2,function(u){
    bdd=cbind(data%*%u,derivs%*%u)
    ddalpha::depth.simplicial(bdd,bdd)})
  return(rowMeans(bd))
}


#' Find Changepoints Using Functional Kruskall-Wallis for Covariance Algorithm
#'
#' @param funcdata Data in fdata form
#' @param depth Depth function of choice
#' @param beta Numeric penalty constant passed to PELT
#'
#' @returns A list of changepoints
#' @export
#'
#' @examples
FKWC <- function(funcdata, depth = "FM_depth", beta = 10){

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
  } #else if (depth == "RPD_depth_d"){
    #derivs = MFHD::derivcurves(funcdata$data)
    #depths = RPDd(funcdata$data,derivs)
  #}
  ranks = rank(depths)
  cp = which(PELT(ranks,length(ranks),beta)==1)-1 #includes the changepoint 0

  if (length(cp) == 1){ #Since first value in cp is 0
    return("No changepoint detected")
  } else {
    return(cp[-c(1)])
  }
}

