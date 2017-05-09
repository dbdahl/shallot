\name{pairwise.probabilities}
\alias{pairwise.probabilities}
\alias{print.shallot.pairwiseProbability}
\alias{as.matrix.shallot.pairwiseProbability}
\title{Pairwise Probabilities}
\description{
These functions relate to the pairwise probabilities that two items are clustered together, e.g., belong to the same subset in a partition.
}
\usage{
pairwise.probabilities(x, parallel = TRUE)
\method{print}{shallot.pairwiseProbability}(x, ...)
\method{as.matrix}{shallot.pairwiseProbability}(x, ...)
}
\arguments{
\item{x}{An object of class \code{shallot.samples.raw} when supplied to \code{\link{pairwise.probabilities}} or an object of class \code{shallot.pairwiseProbability} when supplied to \code{\link{print}} and \code{\link{as.matrix}}.}
\item{parallel}{Should sampling be done in parallel by simultaneously using all CPU cores?}
\item{...}{Currently ignored.}
}
\details{
\code{\link{pairwise.probabilities}} calculates the pairwise probabilities that two items are clustered together, e.g., belong to the same subset in a partition.

\code{\link[=as.matrix.shallot.pairwiseProbability]{as.matrix}} converts the results of \code{\link{pairwise.probabilities}} to an \R matrix.
}
\value{
\code{\link{pairwise.probabilities}} function returns an object of class \code{shallot.pairwiseProbability}.

\code{\link[=as.matrix.shallot.pairwiseProbability]{as.matrix}} return a square matrix.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{sample.partition}}
\code{\link{process.partition}}
\code{\link{estimate.partition}}
}
% DBD: 2017-05-09
