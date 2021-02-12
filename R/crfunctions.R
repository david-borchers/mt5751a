#' Capture-recapture estimation
#' 
#' Estimates abundance using capture-recapture methods. Five methods are implemented: 
#' M0, Mt, Mb, Mtb and Mh (see argument \code{model}.)
#' 
#' @param ch Binary capture history matrix, with rows being individuals and columns
#' being occasions.
#' @param model One of 
#' \itemize{
#'  \item{'M0'}{ Model with constant capture probability.}
#'  \item{'Mt'}{ Model with differewnt capture probability on each occasion.}
#'  \item{'Mb'}{ Model with behavioural response to capture that increases or 
#'  decreases capture probability after first capture.}
#'  \item{'Mtb'}{ Model with differewnt capture probability on each occasion and
#'  behavioural response to capture that increases or decreases capture probability 
#'  after first capture.}
#'  \item{'Mh2'}{ Model a finite 2-part mixture capture probability, i.e. 2 levels of 
#'  latent capture probability.}
#' }
#' @param alpha The significance level for confidence intervals 
#' (defaults to \code{alpha}=0.05 for 95\% level)
#' @param start.n0 Starting value for number of undetected animals for numerical optimisation of likelihood.
#' @param start.p Starting value for detection probability for numerical optimisation of likelihood. For 
#' model Mt all probabilities are assumed to be the same to start with.
#' @param start.b Starting value for the behavioural effect of capture. Positive means capture increases
#' detection probability, negative decreases. (\code{start.b} is the multiplicative effect of previous
#'  capture on the ODDS of capture, i.e. on p/(1-p).)
#' @param start.mu1 Starting value for the first latent capture probability in model Mh2.
#' @param start.mu2 Starting value for the second latent capture probability in model Mh2.
#' @param start.phi Starting value for the proportion of the population that has the first level of the 
#' latent capture probability.
#' @param start.mu1 Starting value for number of undetected animals for numerical optimisation of likelihood.
#' 
#' @return
#' The function returns a list with the following elements:
#' \itemize{
#'  \item{$model :}{ Just reflects the \code{model} argument passed to the function.}
#'  \item{$Nhat :}{ Estimated abundance and 95\% confidence interval.}
#'  \item{$phat :}{ Estimated detection probabilities and 95\% confidence interval(s). NOTE: in the
#'  case of models Mb and Mtb, 'b(odds)' is given - this is the multiplicative effect of previous
#'  capture on the ODDS of capture, i.e. on p/(1-p).}
#'  \item{$beta.point.ests:}{ Point estimates of the model parameters (on the log scale in the 
#'  case of abundance, and on the logit scale in the case of probabilities.).}
#'  \item{$beta.varcovar.ests:}{ The estimated variance-covariance matrix of the model parameters 
#'  (on the log scale in the case of abundance, and on the logit scale in the case of probabilities.)}
#'  \item{$$beta.corrmatrix:}{ The estimated correlation matrix of the model parameters 
#'  (on the log scale in the case of abundance, and on the logit scale in the case of probabilities.).}
#'  \item{$nobs:}{ Number of unique captures.}
#'  \item{$loglik:}{ The log-likelihood at the MLE.}
#'  \item{$AIC:}{ Akaike's information criterion for the model.}
#' }
#' 
#' @export
fit.cr = function(ch, model,alpha=0.05,
                  start.n0=NULL,start.p=NULL,start.b=NULL,
                  start.mu1=NULL,start.mu2=NULL,start.phi=NULL) {
    chmat = as.matrix(ch)
    n = nrow(chmat)
    K = ncol(chmat)
    if(is.null(start.p)) start.p = 0.5
    if(is.null(start.b)) start.b = 0.5
    if(is.null(start.n0)) start.n0 = n/start.p - n
    if(is.null(start.mu1)) start.mu1=start.p
    if(is.null(start.mu2)) start.mu2=start.p
    if(is.null(start.phi)) start.phi=0.5
    prevcap = timecov = NULL
    
    # construct pars and timecov, prevcap
    if(model=="Mh2") { # here for model Mh2
        pars = c(log(start.n0),qlogis(c(start.mu1,start.mu2,start.phi)))
    } else if(model=="M0") {
        pars = c(log(start.n0),qlogis(start.p))
    } else if(model=="Mt") {
        #        if(time.factor) timecov = model.matrix(~occ,data=data.frame(occ=as.factor(1:K)))
        #        else timecov = model.matrix(~occ,data=data.frame(occ=1:K))
        timecov = model.matrix(~occ,data=data.frame(occ=as.factor(1:K)))
        pars = c(log(start.n0),rep(qlogis(start.p),K))
    } else if(model=="Mb") {
        prevcap = make.prevcap(chmat)
        pars = c(log(start.n0),qlogis(start.p),exp(start.b))
    } else if(model=="Mtb") {
        #        if(time.factor) timecov = model.matrix(~occ,data=data.frame(occ=as.factor(1:K)))
        #        else timecov = model.matrix(~occ,data=data.frame(occ=1:K))
        timecov = model.matrix(~occ,data=data.frame(occ=as.factor(1:K)))
        prevcap = make.prevcap(chmat)
        pars = c(log(start.n0),rep(qlogis(start.p),K),exp(start.b))
    } else if(model=="Mh2") {
        pars = c(log(start.n0),qlogis(start.mu1),qlogis(start.mu2),qlogis(start.phi))
    } else stop("Invalid model; must be 'M0' or 'Mt' or 'Mb' or 'Mtb' or 'Mh2'.")
    
    if(model=="Mh2") {
        start.nll = loglik.Mh2(pars, chmat)
        if(is.na(start.nll) | is.infinite(start.nll)) stop("Infeasible starting value(s).")
        mle = nlm(loglik.Mh2 , pars, chmat, hessian=TRUE) 
    } else {
        start.nll = loglik.Mtb(pars, chmat, prevcap, timecov)
        if(is.na(start.nll) | is.infinite(start.nll)) stop("Infeasible starting value(s).")
        mle = nlm(loglik.Mtb , pars, chmat, hessian=TRUE, prevcap=prevcap, timecov=timecov) 
    }
    
    # Unpack N estimate
    Nhat<- nrow(ch) + exp(mle$estimate[1])
    vcv = solve(mle$hessian)
    se.Nhat = sqrt(delta.var(mle$estimate[1],vcv[1,1],"log"))
    ci.Nhat = nrow(ch) + exp(norm.ci(mle$estimate[1],sqrt(vcv[1,1]),alpha))
    
    # name parameters and get CIs for det pars
    if(model=="M0") {
        parnames = c("N","p")
        pests = make.pests(parnames,mle$estimate,vcv,alpha=alpha)
    } else if(model=="Mt") {
        parnames = c("N",paste("p",1:K,sep=""))
        pests = make.pests(parnames,mle$estimate,vcv,alpha=alpha)
    } else if(model=="Mb") {
        parnames = c("N","p","b")
        pests = make.pests(parnames,mle$estimate,vcv,b=TRUE,alpha=alpha)
    } else if(model=="Mtb") {
        parnames = c("N",paste("p",1:K,sep=""),"b")
        pests = make.pests(parnames,mle$estimate,vcv,b=TRUE,alpha=alpha)
    } else if(model=="Mh2") {
        parnames = c("N","p1","p2","phi")
        pests = make.pests(parnames,mle$estimate,vcv,alpha=alpha)
    }
    
    point.ests = mle$estimate
    names(point.ests) = parnames
    row.names(vcv) = colnames(vcv) = parnames
    
    outp =list(
        model = model,
        Nhat = c(Nhat=Nhat,
                 se.Nhat=se.Nhat,
                 lcl.Nhat=ci.Nhat[1],
                 ucl.Nhat=ci.Nhat[2]),
        phat = pests,
        beta.point.ests = point.ests,
        beta.varcovar.ests = vcv,
        beta.corrmatrix = cov2cor(vcv),
        nobs = n,
        loglik = -1*mle$minimum,
        AIC = 2*mle$minimum + 2*length(mle$estimate))
    class(outp) = c("crfit",class(outp))
    return(outp)
}


#' Negative log-likelihood for models M0, Mt, Mb, Mtb
#' 
#' Calculates the negative log-likelihood for models M0, Mt, Mb, Mtb.
#' 
#' @param pars Parameter vector with the correct number of parameters in the correct order 
#' (see Details).
#' @param chmat Binary capture history matrix, with one row per individual detected and 
#' one column per capture occasion.
#' @param prevcap Binary matrix of same size as \code{chmat}, with 0 indicating that the
#' individual has not been captured by that occasion, and 1 indicating that it has. Used
#' by the function to decide what the model is, so if it is not null, the model is taken
#' to be Mb (if \code{timecov} is NULL) or Mtb (if \code{timecov} is not NULL).
#' @param timecov Matrix of same size as \code{chmat} indicating the time associated with the
#' capture occasion. At present it is treated as a factor.
#' 
#' @details 
#' The function decides which model to use on the basis of which of \code{prevcap} and 
#' \code{timecov} are null. It assumes that \code{pars} has the appropriate number of parameters
#' for the model, and in the appropriate order, as follows:
#' \itemize{
#'  \item{'M0'}{ \code{prevcap} and \code{timecov} both NULL: \code{pars} must be 
#'  \code{(log(N),qlogis(p))}}.
#'  \item{'Mt'}{ \code{prevcap} NULL, \code{timecov} not: \code{pars} must be 
#'  \code{(log(N),qlogis(p1)...qlogis(pK))}, where there are \code{K} occasions.}
#'  \item{'Mb'}{ \code{prevcap} not NUL, \code{timecov} is NULL: \code{pars} must be 
#'  \code{(log(N),qlogis(p),log(b))}, where \code{log(b)} is the change in log-odds due to 
#'  previous capture.}
#'  \item{'Mtb'}{ Neither \code{prevcap} nor \code{timecov} are NULL: \code{pars} must be 
#'  \code{(log(N),qlogis(p1)...qlogis(pK),log(b))}, where there are \code{K} occasions and 
#'  \code{log(b)} is the change in log-odds due to previous capture.}
#' }
#' 
#' @return
#' The negative log-likelihood for the appropriate model.
#' 
#' @export
loglik.Mtb <- function(pars, chmat, prevcap=NULL, timecov=NULL ) {

    K = ncol(chmat)
    n = nrow(chmat)
    npar = length(pars)
    
    # deal with abundance parameter:
    n0 <- exp(pars[1])
    N <- n + n0
    
    chmat <- rbind(chmat, rep(0, K)) # add row of zeros for undetected
    nplus = nrow(chmat)
    
    # M0
    if(is.null(prevcap) & is.null(timecov)) {
        lp = matrix(rep(pars[2],nplus*K),nrow=nplus)
    }
    # Mt
    if(is.null(prevcap) & !is.null(timecov)) {
        lp = matrix(rep(timecov%*%pars[2:(K+1)],nplus),ncol=K,byrow=TRUE)
    }
    # Mb
    if(!is.null(prevcap) & is.null(timecov)) {
        prevcap = rbind(prevcap,rep(0,K)) # add row of zeros for undetected
        lp = rep(pars[2],nplus) + prevcap*pars[npar]
    }
    # Mtb
    if(!is.null(prevcap) & !is.null(timecov)) {
        prevcap = rbind(prevcap,rep(0,K)) # add row of zeros for undetected
        lp = matrix(rep(timecov%*%pars[2:(K+1)],nplus),ncol=K,byrow=TRUE) + prevcap*pars[npar]
    }
    
    # Calculate p with linear predictor
    p <- plogis(lp)
    
    # calculate negatiev log-likelihood
    negllik <- rep(NA, nplus)
    for (i in 1:nplus) {
        negllik[i] <- sum(dbinom(chmat[i, ], 1, p[i, ], log = TRUE))
    }
    nll = -1 * (lgamma(N + 1) - lgamma(n0 + 1) + sum(c(rep(1, n),  n0) * negllik))
    
    return(nll)
}

    

#' Negative log-likelihood for model Mh2
#' 
#' Calculates the negative log-likelihood for a 2-part finite mixture model.
#' 
#' @param pars Parameter vector comprising (qlogis(p1), qlogis(p2), qlogis(phi)), where
#' p1 and p2 are the two latent capture probabilities, and phi is the proportion of the population
#' with capture probability p1.
#' @param chmat Binary capture history matrix, with one row per individual detected and 
#' one column per capture occasion.
#' 
#' @return
#' The negative log-likelihood.
#' 
#' @export
loglik.Mh2 <- function(pars, chmat) {
    mu1 <- pars[2]
    mu2 <- pars[3]
    p1 <- exp(mu1)/(1 + exp(mu1))
    p2 <- exp(mu2)/(1 + exp(mu2))
    psi <- exp(pars[4])/(1 + exp(pars[4]))
    n0 <- exp(pars[1])
    
    nind <- nrow(chmat)
    K <- ncol(chmat)
    chmat <- rbind(chmat, rep(0, K))
    ml <- rep(NA, nrow(chmat))
    for (i in 1:nrow(chmat)) {
        ml[i] <- exp(sum(dbinom(chmat[i, ], 1, p1, log = TRUE))) * 
            psi + exp(sum(dbinom(chmat[i, ], 1, p2, log = TRUE))) * (1 - psi)
    }
    nll = -1 * (lgamma(n0 + nind + 1) - lgamma(n0 + 1) + sum(c(rep(1,  nind), n0) * log(ml)))
    return(nll)
}



#' Make binary covariate for previous capture
#' 
#' Makes an indicator variable for whether or not an animal has been captured before
#' each occasion.
#' 
#' @param chmat Capture history matrix.
#' 
#' @return
#' The function returns a matrix with one row per individual, one column per occasion, with 
#' a 0 indicating the animal had not been captured before an a 1 indicating that it had.
#' 
#' @export
make.prevcap = function(chmat) {
    prevcap = chmat*0
    K = ncol(chmat)
    for(k in 2:K) prevcap[,k] = as.numeric(rowSums(as.matrix(chmat[,1:(k-1)]))>0)
    return(prevcap)
}


#' Make matrix of probability estimates
#' 
#' Make matrix of probability estimates on the natural scale, from estimates on the linke
#' function scale, and their variance-covariance matrix.
#' 
#' @param parnames Names of the parameters in \code{ests}.
#' @param ests Vector of parameter estimates on the link function scale.
#' @param vcv The variance-covariance matrix of \code{ests}.
#' @param b Indicator variable saying whether or not the last parameter in \code{ests} is
#' the behavioural effect parameter, in which case it is interpreted as change in log-odds.
#' @param alpha Confidence interval level.
#' 
#' @return
#' The function returns a matrix with one row per parameter and columns for the point estimate
#' on the natural scale and associated confidence intervals.
#' 
#' @export
make.pests = function(parnames,ests,vcv,b=FALSE,alpha=0.05) {
    np = length(parnames)-1
    phat = plogis(ests[2:(1+np)])
    pests = matrix(c(phat,rep(NA,length(phat)*3)),ncol=4)
    row.names(pests) = parnames[2:(np+1)]
    colnames(pests) = c("est","se","lcl","ucl")
    for(i in 1:np) {
        if(i==np & b) { # report b as ODDS
            row.names(pests)[np] = "b(odds)"
            pests[np,1] = exp(ests[np+1])
            pests[i,2] = sqrt(delta.var(ests[np+1],vcv[np+1,np+1],"log"))
            pests[np,3:4] = exp(norm.ci(ests[np+1],sqrt(vcv[np+1,np+1]),alpha))
        } else {
            pests[i,2] = sqrt(delta.var(ests[i+1],vcv[i+1,i+1],"logit"))
            pests[i,3:4] = plogis(norm.ci(ests[i+1],sqrt(vcv[i+1,i+1]),alpha))
        }
    }
    return(pests)
}

#' Delta method variance
#' 
#' Calculated Delta Method variance for inverse log and inverse logistic link functions.
#' 
#' @param x Random variable that inverse transformation is applied to
#' @param var.x Variance of \code{x}.
#' @param link Link function ('log' or 'logit').
#' 
#' @return
#' The Delta Method (first order Taylor series) approximation.
#' 
#' @export
delta.var = function(x,var.x,link) {
    if(link=="logit") {
        ex = exp(x)
        deltavar = var.x*(ex/(1+ex)^2)^2
    }else if(link=="log") {
        deltavar = exp(x)^2*var.x
    }else stop("Invalid link; only 'logit' and 'log' implemented.")
    return(deltavar)
}
