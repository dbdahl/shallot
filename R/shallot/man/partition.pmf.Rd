\name{partition.pmf}
\alias{partition.pmf}
\title{Obtain the Probability Mass Function of a Partition Distribution}
\description{
This function returns the probability mass function (pmf) of a partition distribution.
}
\usage{
partition.pmf(x)
}
\arguments{
\item{x}{An object of class \code{shallot.distribution} obtained, for example, from the \code{\link{ewens.pitman.attraction}} function.}
}
\value{
A function that takes a partition (as a vector in cluster label notation) and returns the probability --- or, if \code{log=TRUE}, the log of the probability --- of the supplied partition.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\examples{
\dontrun{
example(shallot)
}
}
