\name{shallot-package}
\alias{shallot-package}
\alias{shallot}
\docType{package}
\title{
Random Partition Distribution Indexed by Pairwise Information
}
\description{
This package implements models described in the paper 
\href{https://doi.org/10.1080/01621459.2016.1165103}{Dahl, D. B., Day, R., and Tsai, J. (2017), Random Partition Distribution Indexed by Pairwise Information, \emph{Journal of the American Statistical Association}, accepted.}
The Ewens, Ewens-Pitman, Ewens attraction, Ewens-Pitman attraction, and ddCRP distributions are available for prior simulation.  We hope in the future to add posterior simulation with a user-supplied likelihood.  Supporting functions for partition estimation and plotting are also planned.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\references{
\href{https://doi.org/10.1080/01621459.2016.1165103}{Dahl, D. B., Day, R., and Tsai, J. (2017), Random Partition Distribution Indexed by Pairwise Information, \emph{Journal of the American Statistical Association}, accepted. <DOI:10.1080/01621459.2016.1165103>}
}
\seealso{
  \code{\link{ewens.pitman.attraction}},
  \code{\link{sample.partitions}},
  \code{\link{estimate.partition}}
}
\examples{
\donttest{
data <- iris[,-ncol(iris)]
truth <- as.integer(iris[,ncol(iris)])
distance <- as.dist(as.matrix(dist(scale(data))+0.001))

decay <- decay.exponential(temperature(9.0, fixed=TRUE), distance)
permutation <- permutation(n.items=nrow(data), fixed = FALSE)
attraction <- attraction(permutation, decay)
mass <- mass(1.0, fixed = TRUE)
discount <- discount(0.2, fixed = TRUE)
distribution <- ewens.pitman.attraction(mass, discount, attraction)

raw <- sample.partitions(distribution, 500, parallel=FALSE)
samples <- process.samples(raw)

library(sdols)
pp <- expectedPairwiseAllocationMatrix(samples$labels)
est <- salso(pp)
conf <- confidence(est,pp)
plot(conf)
plot(conf,data=data)

\dontshow{
rscala::scalaDisconnect(shallot:::s)
}
}
}
\keyword{package}
