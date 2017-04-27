\name{process.partition}
\alias{process.partition}
\alias{print.shallot.samples}
\title{Process}
\description{
  A function that unpacks the number of partitions and the number of subsets per sample from the collected results.}
\usage{
process(x, ...)
\method{print}{shallot.samples}(x, ...)
}
\arguments{
\item{x}{An object of class \code{shallot.samples.raw}.}
\item{...}{An unused argument.}
}
\details{
This function, though not really necessary in the shallot program, allows the user to see the distribution of the number of subsets across the samples prior to factoring in maxSubsets and maxScans in estimating the partition.  The function also processes the data further if an internal variable from the \code{\link{sample.partition}} function called augmentedSamples is not \code{NULL}.

The \code{\link{print.shallot.samples}} function is an internal function that prints a message listing the attributes of the results of the \code{\link{process.partition}} function.
}
\value{
An object of class \code{shallot.samples} storing the number of partitions and the number of subsets per sample.
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

