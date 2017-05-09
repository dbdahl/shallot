\name{process.partition}
\alias{process.partition}
\title{Process Sampled Partitions}
\description{
  This function extracts the partitions from the results of the \code{\link{sample.partition} function.}
}
\usage{
process.partition(x, sample.parameter = NULL, as.matrix = FALSE)
}
\arguments{
\item{x}{An object from the \code{\link{sample.partition}} function.}
\item{sample.parameter}{If not \code{NULL}, a function taking no arguments and returning a random sample from the centering distribution.}
\item{as.matrix}{Should the results be returned as a matrix of cluster labels or lists of partitions?}
}
\details{
This function extracts the sampled partitions from the results of the \code{\link{sample.partition} function.}
}
\value{
A list containing either: i. a matrix of cluster labels in which each row
represents a clusterings, or ii. a list containing partitions.  The list also
contains sampled model parameters if \code{sample.parameter} is not
\code{NULL}.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{sample.partition}}
}
\examples{
\dontrun{
## Demonstrate process.
data <- iris[,-ncol(iris)]
truth <- as.integer(iris[,ncol(iris)])
distance <- as.dist(as.matrix(dist(scale(data))+0.001))

decay <- decay.exponential(9.0, fixed = TRUE)
attract <- attraction(distance, permutation(n.items = nrow(data), fixed = FALSE), decay)
mass <- mass(1.0, fixed = TRUE)
discount <- discount(0.2, fixed = TRUE)
distribution <- ewens.pitman.attraction(mass, discount, attract)

samples.raw <- sample.partition(distribution, 500)
samples <- process.partition(samples.raw,as.matrix=TRUE)

est <- estimate.partition(samples.raw,max.subsets=3,max.scans=0)
pp <- pairwise.probabilities(samples.raw)

conf <- confidence(pp,est)
plot(conf)
plot(conf,data=data)
adj.rand.index(est,truth)
}
}

