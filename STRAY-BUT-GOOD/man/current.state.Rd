\name{current.state}
\alias{current.state}
\alias{print.shallot.state}
\title{Current State}
\description{
  A function that displays the partition in its current state.}
\usage{
current.state(x)
\method{print}{shallot.state}(x, ...)
}
\arguments{
\item{x}{An object of class \code{shallot.samples.raw}.}
\item{...}{An additional argument, not used.}
}
\details{
This function allows the user to see the partition in its current state prior to factoring in maxSubsets and maxScans in estimating the partition.  The function also processes the data further based on whether or not \code{\link{mcmc.parameters}} and an internal variable from the \code{\link{collect}} function called augmentedSamples are \code{NULL}.

The \code{\link{print.shallot.state}} function is an internal function that prints the partition - that is, it prints which numbered subset each observation is assigned to.
}
\value{
An object of class \code{shallot.state} storing the partition.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{collect}}
}
\examples{
\dontrun{
## Demonstrate current.state.
data <- Theoph[, c("Wt", "Dose", "Time", "conc", "Subject")]
truth <- as.numeric(levels(data[, ncol(data)]))[data[, ncol(data)]]
distance <- dist(scale(data[, -ncol(data)]))
decay <- decay.exponential(9.0, fixed = TRUE)
attract <- attraction(distance(as.matrix(distance)),
    permutation(n.items = nrow(data), fixed = FALSE), decay)
mass <- mass(1.0, fixed = TRUE)
distribution <- ewens(mass, n.items = attract$n.items)
samples <- collect(distribution, 500)
print(current.state(samples))
}
}

