#' Normal confidence interval
#' 
#' Returns confidence interval assuming estimator has a normal distribution
#' 
#' @param x The estimate
#' @param se The standard error of the estimate
#' @param alpha The significance level (defaults to \code{alpha}=0.05 for 95\% level)
#' 
#' @return  
#' stuff
#' 
#' @export
norm.ci = function(x,se,alpha=0.05) {
  x + qnorm(1-alpha/2)*c(-1,1)*se
}

#' Lognormal confidence interval
#' 
#' Returns confidence interval assuming estimator has a lognormal distribution
#' 
#' @param x The estimate
#' @param se The standard error of the estimate
#' @param alpha The significance level (defaults to \code{alpha}=0.05 for 95\% level)
#' 
#' @return
#' The function returns a list with the following elements:
#' 
#' @export
lognorm.ci = function(x,se,alpha=0.05) {
  x*exp(qnorm(1-alpha/2)*c(-1,1)*sqrt(log(1+(se/x)^2)))
}

