\name{sampling.model}
\alias{sampling.model}
\title{Specifying a Sampling Model}
\description{
These functions set the mass, discount, and temperature parameters and, in the case of them being random, specify the parameters of their distribution.}
\usage{
sampling.model(sample.parameter, log.density)
}
\arguments{
\item{sample.parameter}{A function taking two arguments with names \code{indices} and \code{parameter} having default values \code{c()} and \code{NULL}, respectively.  With those default values, the function should return a sample from the centering distribution as an \R object of any type.  Otherwise, \code{indices} is a vector of integers indicating the cluster elements and \code{parameter} is the current value of the parameter for the cluster, and the function should return an updated
value for the cluster parameter based on a valid MCMC update.}
\item{log.density}{A function taking three arguments with names \code{i}, \code{indices}, and \code{parameter} and returning a double giving the log of the likelihood contribution of item \code{i} to a clustering with members \code{indices} and parameter \code{parameter}.}
}
\value{
An object of class \code{shallot.distribution.data}.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\examples{
## Model inputs.
data <- c(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
sigma  <- 0.1
mu0    <- 0.0
sigma0 <- 1.0

## Derived values.
s2 <- sigma * sigma
s02 <- sigma0 * sigma0
s02Inv <- 1.0 / s02
c <- -1.0 / (2.0 * s2)

## Sampling model of Neal (JCGS, 2009)
## Function to perform an MCMC update of the parameter.
sample.parameter <- function(indices=c(), parameter=NULL) {
  sum <- sum(data[indices])
  variance <- 1 / (s02Inv + length(indices) / s2)
  mean <- variance * (mu0 / s02 + sum / s2)
  rnorm(1, mean=mean, sd=sqrt(variance))
}

## Function to evaluate the likelihood contribution for an observation.
log.density <- function(i, indices, parameter) {
  resid <- data[i] - parameter
  c * resid * resid
}

sm <- sampling.model(sample.parameter, log.density)
sm
}
