% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crfunctions.R
\name{fit.cr}
\alias{fit.cr}
\title{Capture-recapture estimation}
\usage{
fit.cr(
  ch,
  model,
  alpha = 0.05,
  start.n0 = NULL,
  start.p = NULL,
  start.b = NULL,
  start.mu1 = NULL,
  start.mu2 = NULL,
  start.phi = NULL
)
}
\arguments{
\item{ch}{Binary capture history matrix, with rows being individuals and columns
being occasions.}

\item{model}{One of 
\itemize{
 \item{'M0'}{ Model with constant capture probability.}
 \item{'Mt'}{ Model with differewnt capture probability on each occasion.}
 \item{'Mb'}{ Model with behavioural response to capture that increases or 
 decreases capture probability after first capture.}
 \item{'Mtb'}{ Model with differewnt capture probability on each occasion and
 behavioural response to capture that increases or decreases capture probability 
 after first capture.}
 \item{'Mh2'}{ Model a finite 2-part mixture capture probability, i.e. 2 levels of 
 latent capture probability.}
}}

\item{alpha}{The significance level for confidence intervals 
(defaults to \code{alpha}=0.05 for 95\% level)}

\item{start.n0}{Starting value for number of undetected animals for numerical optimisation of likelihood.}

\item{start.p}{Starting value for detection probability for numerical optimisation of likelihood. For 
model Mt all probabilities are assumed to be the same to start with.}

\item{start.b}{Starting value for the behavioural effect of capture. Positive means capture increases
detection probability, negative decreases. (\code{start.b} is the multiplicative effect of previous
 capture on the ODDS of capture, i.e. on p/(1-p).)}

\item{start.mu1}{Starting value for number of undetected animals for numerical optimisation of likelihood.}

\item{start.mu2}{Starting value for the second latent capture probability in model Mh2.}

\item{start.phi}{Starting value for the proportion of the population that has the first level of the 
latent capture probability.}
}
\value{
The function returns a list with the following elements:
\itemize{
 \item{$model :}{ Just reflects the \code{model} argument passed to the function.}
 \item{$Nhat :}{ Estimated abundance and 95\% confidence interval.}
 \item{$phat :}{ Estimated detection probabilities and 95\% confidence interval(s). NOTE: in the
 case of models Mb and Mtb, 'b(odds)' is given - this is the multiplicative effect of previous
 capture on the ODDS of capture, i.e. on p/(1-p).}
 \item{$beta.point.ests:}{ Point estimates of the model parameters (on the log scale in the 
 case of abundance, and on the logit scale in the case of probabilities.).}
 \item{$beta.varcovar.ests:}{ The estimated variance-covariance matrix of the model parameters 
 (on the log scale in the case of abundance, and on the logit scale in the case of probabilities.)}
 \item{$$beta.corrmatrix:}{ The estimated correlation matrix of the model parameters 
 (on the log scale in the case of abundance, and on the logit scale in the case of probabilities.).}
 \item{$nobs:}{ Number of unique captures.}
 \item{$loglik:}{ The log-likelihood at the MLE.}
 \item{$AIC:}{ Akaike's information criterion for the model.}
}
}
\description{
Estimates abundance using capture-recapture methods. Five methods are implemented: 
M0, Mt, Mb, Mtb and Mh (see argument \code{model}.)
}
