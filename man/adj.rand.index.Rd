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
The adjusted Rand index is at most \code{1.0} and large numbers indicate high similarity.
}
\value{
A numeric vector of length one representing the adjusted Rand index between the two clusterings/partitions.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{estimate.partition}}
}
\references{
Hubert, L. and Arabie, P. (1985), Comparing partitions, \emph{Journal of Classification}, \bold{2}, 193-218.
}
\examples{
truth    <- c(1,1,2,2,2,1,3,3,3)
estimate <- c(1,2,2,2,2,1,2,3,3)
adj.rand.index(truth,estimate)
}
% DBD: 2017-05-02
