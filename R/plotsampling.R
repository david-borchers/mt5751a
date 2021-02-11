#' Plot sampling estimation
#' 
#' Estimates abundance as n/p, where n is total count and p inclusion probability.
#' Also returns a confidence interval, estimated standard error and CV, using one of
#' four possible methods (see argument \code{method}.)
#' 
#' @param n Vector of counts in plots
#' @param p Proportion of survey region covered by plots
#' @param alpha The significance level (defaults to \code{alpha}=0.05 for 95\% level)
#' @param method If 'normal', the estimator is assumed to be normally distributed with 
#' variance calculated assuming the distribution given by \code{dbn}, if 
#' 'lognormal', the estimator is assumed to be lognormally distributed with 
#' variance calculated assuming the distribution given by \code{dbn}, if 'exact', the 
#' estimator is assumed to have the parameteric distribution specified by \code{dbn}. If
#' 'percentile' and \code{dbn} is 'bootstrap', the percentile method is used.
#' @param dbn If 'binomial', the count is assumed to be binomially distributed (for 
#' variance estimation, and if \code{method} is 'exact', then for confidence interval 
#' estimation too.  If 'binomial', the count is assumed to be binomially distributed (for 
#' variance estimation, and if \code{method} is 'exact', then for confidence interval 
#' estimation too.  If 'bootstrap' the variance is obtained by nonparametric bootstrap 
#' of the plots, and the confidence interval is obtained using this and assuming that the 
#' estimator is normally distributed (if \code{method} is 'bionomial'), or  assuming that the 
#' estimator is Poisson distributed (if \code{method} is 'poisson'),  or making no 
#' distributional assumption and using the pecentile method (if \code{method} is 'percentile')
#' @param B Number of bootstrap replicates.
#' @param Nmult multiple of estimate beyond which probability of getting observed data
#' is assumed to be zero (just for computational covenience when calculating exact CI).
#' 
#' @return
#' The function returns a list with the following elements:
#' \itemize{
#'  \item{$Nhat :}{ Estimated abundance.}
#'  \item{$se.Nhat :}{ Standard error of the estimator.}
#'  \item{$cv.Nhat :}{ Coefficient of variance of the estimator.}
#'  \item{$ci.Nhat:}{ Confidence interval.}
#' }
#' 
#' @export
plotsample_est = function(n,p,alpha=0.05,method="lognormal",dbn="binomial",B=999,Nmult=3) {
  ntot = sum(n)
  Nhat = ntot/p
  if(dbn=="bootstrap") {
    b.Nhat = rep(NA,B)
    for(i in 1:B) b.Nhat[i] = sum(sample(n,size=length(n),replace=TRUE))/p
    se.Nhat = sqrt(var(b.Nhat))
    Nhat = sum(n)/p
    cv.Nhat = se.Nhat/mean(b.Nhat)
    if(method=="normal") {
      ci.Nhat = round(norm.ci(Nhat,se.Nhat,alpha))
    }else if(method=="lognormal") {
      ci.Nhat = round(lognorm.ci(Nhat,se.Nhat,alpha))
    }else if(method=="percentile") {
      ci.Nhat = round(quantile(b.Nhat,probs=c(alpha/2,1-alpha/2)))
    }else {
      stop("Invalid method. For boostrap, method must be 'normal', 'lognormal' or 'percentile'.")
    }
  }else {
    if(dbn=="binomial") {
      var.Nhat = Nhat*(1-p)/p
    }else if (dbn=="poisson") {
      var.Nhat = Nhat/p
    }else {
      stop("Invalid dbn. It must be 'binomial', 'poisson' or 'bootstrap.")
    }
    se.Nhat = sqrt(var.Nhat)
    cv.Nhat = se.Nhat/Nhat
    if(method=="normal") {
      ci.Nhat = round(norm.ci(Nhat,se.Nhat,alpha))
    }else if(method=="lognormal") {
      ci.Nhat = round(lognorm.ci(Nhat,se.Nhat,alpha))
    }else if(method=="exact") {
      Ns = c(ntot:(Nmult*Nhat))
      cdf = rep(NA,length(Ns))
      if(dbn=="binomial") for(i in 1:length(Ns)) cdf[i] = pbinom(ntot,Ns[i],p)
      else if(dbn=="poisson") for(i in 1:length(Ns)) cdf[i] = ppois(ntot,Ns[i]*p)
      else stop("Invalid dbn. With method 'exact', dbn must be 'binomial' or 'poisson'.")
      lo = Ns[max(which(cdf>=(1-alpha/2)))]
      up = Ns[min(which(cdf<=(alpha/2)))]
      ci.Nhat = round(c(lo,up))
    }
  }
  return(list(Nhat=Nhat,se.Nhat=se.Nhat,cv.Nhat=cv.Nhat,ci.Nhat=ci.Nhat,method=method))
}


