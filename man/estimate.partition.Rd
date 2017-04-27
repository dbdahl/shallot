\name{estimate.partition}
\alias{estimate.partition}
\title{Estimate Partition}
\description{
    A function that calculates the partition by factoring in maxSubsets and maxScans.}
\usage{
estimate.partition(x, pairwise.probabilities = NULL,
    max.subsets = 0, max.scans = 0, parallel = TRUE)
}
\arguments{
\item{x}{An object of class \code{shallot.samples.raw}.}
\item{pairwise.probabilities}{An object of class \code{shallot.pairwiseProbability}.  Defaults to \code{NULL}.}
\item{max.subsets}{An integer listing the maximum number of subsets.  Defaults to \code{0}.}
\item{max.scans}{An integer listing the number of times to scan.  Defaults to \code{0}, can also be set to \code{10}.}
\item{parallel}{A boolean object indicating whether all CPU cores should be used.  Defaults to \code{TRUE}.}
}
\details{
This function takes the given partition and factors in maxSubsets and maxScans to determine what the partition looks like given these two additional variables.  The variable known as \var{max.subsets} represents the number of subsets established by the truth column (although it is recommended to use the default value of \code{0}), while \var{max.scans} represents the number of times to scan.  Output does not appear to be altered by setting \var{parallel} to \code{FALSE}.  Output is altered slightly by including \var{pairwise.probabilities} as an input value.
}
\value{
An array listing the partition - that is, indicating which numbered subset each element is assigned to.  This array is then one of the two input values for \code{\link{adj.rand.index}} (the other being the truth).
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{sample.partition}}
\code{\link{plotting}}
\code{\link{adj.rand.index}}
}
\examples{
\donttest{
## Demonstrate estimate.partition.
data <- Theoph[, c("Wt", "Dose", "Time", "conc", "Subject")]
truth <- as.numeric(levels(data[, ncol(data)]))[data[, ncol(data)]]
distance <- dist(scale(data[, -ncol(data)]))
decay <- decay.exponential(9.0, fixed = TRUE)
mass <- mass(1.0, fixed = TRUE)
discount <- discount(0.1, fixed = TRUE)
distribution <- ewens.pitman(mass, discount, n.items=nrow(data))
samples <- sample.partition(distribution, 500)
est <- estimate.partition(samples, max.subsets = 12, max.scans = 10)
print(est)
}
}

