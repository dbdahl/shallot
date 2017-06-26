\name{pairwise.probabilities}
\alias{pairwise.probabilities}
\alias{print.shallot.pairwiseProbability}
\alias{as.matrix.shallot.pairwiseProbability}
\title{Pairwise Probabilities}
\description{
These functions relate to the pairwise probabilities that two items are clustered together, i.e., belong to the same subset in a partition.
}
\usage{
pairwise.probabilities(x, parallel = TRUE)
\method{print}{shallot.pairwiseProbability}(x, ...)
\method{as.matrix}{shallot.pairwiseProbability}(x, ...)
}
\arguments{
\item{x}{An object of class \code{shallot.samples.raw} when supplied to \code{\link{pairwise.probabilities}} or an object of class \code{shallot.pairwiseProbability} when supplied to \code{\link{print}} and \code{\link{as.matrix}}.}
\item{parallel}{Should all of the CPU cores should be used?  Defaults to \code{TRUE}.}
\item{...}{Currently ignored.}
}
\details{
\code{\link{pairwise.probabilities}} calculates the pairwise probabilities that two items are clustered together, i.e., belong to the same subset in a partition.

\code{\link[=as.matrix.shallot.pairwiseProbability]{as.matrix}} converts the results of \code{\link{pairwise.probabilities}} to an \R matrix.
}
\value{
The \code{\link{pairwise.probabilities}} function returns an object of class \code{shallot.pairwiseProbability}.

The \code{\link[=as.matrix.shallot.pairwiseProbability]{as.matrix}} function returns a square matrix.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{sample.partitions}},
\code{\link{process.partitions}},
\code{\link{estimate.partition}}
}
\examples{
\dontrun{
example(shallot)
}
}
