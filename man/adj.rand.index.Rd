\name{adj.rand.index}
\alias{adj.rand.index}
\title{Adjusted Rand Index}
\description{
This function calculates the adjusted Rand index between two clusterings/partitions.}
\usage{
adj.rand.index(c1, c2)
}
\arguments{
\item{c1}{A vector containing cluster labels for a clustering/partition.}
\item{c2}{A vector containing cluster labels for a clustering/partition.}
}
\details{
The \code{\link{adj.rand.index}} function takes as its input two clusterings/partitions in cluster label notation and computes the adjusted Rand index.
The adjusted Rand index is at most \code{1.0} and larger number indicate high similarity.
}
\value{
A numeric vector of length one representing the adjusted Rand index between the two clusterings/partition.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{estimate.partition}}
}
\examples{
## Load example dataset.
truth    <- c(1,1,2,2,2,1,3,3,3)
estimate <- c(1,2,2,2,2,1,2,3,3)
adj.rand.index(truth,estimate)
}

