\name{plot.partition}
\alias{plot.partition}
\alias{confidence}
\alias{plot.shallot.confidence}
\title{Confidence and Exemplar Plotting}
\description{
Functions to produce confidence plots (e.g., heatmaps of pairwise probabilities of clustering) and exemplar plots.
}
\usage{
confidence(pairwise.probabilities, partition)
\method{plot}{shallot.confidence}(x, partition = NULL, data = NULL,
  show.labels = length(x$partition) <= 50, ...)
}
\arguments{
\item{partition}{A vector of cluster labels, or \code{NULL}.}
\item{x}{An object of class \code{shallot.confidence}.}
\item{pairwise.probabilities}{An object of class \code{shallot.pairwiseProbability}, or \code{NULL}.}
\item{data}{The data from which the distances were computed.}
\item{show.labels}{Show the items names be shown in the plot?}
\item{...}{Currently ignored.}
}
\value{
The \code{\link{confidence}} function outputs an object of class \code{shallot.confidence} being a list of confidence information.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{sample.partition}}
\code{\link{process.partition}}
\code{\link{estimate.partition}}
}
\examples{
\dontrun{
## Demonstrate pairwise.probabilities, as.matrix, confidence, and plot.
data <- Theoph[, c("Wt", "Dose", "Time", "conc", "Subject")]
truth <- as.numeric(levels(data[, ncol(data)]))[data[, ncol(data)]]
distance <- dist(as.matrix(dist(scale(data[, -ncol(data)]))+0.001))
decay <- decay.exponential(temperature(9.0, fixed = TRUE), distance)
attract <- attraction(permutation(n.items = nrow(data), fixed = FALSE), decay)
mass <- mass(1.0, fixed = TRUE)
discount <- discount(0.1, fixed = TRUE)
distribution <- ewens.pitman.attraction(mass, discount, attract)
samples <- sample.partition(distribution, 500)
est <- estimate.partition(samples, max.subsets = 12, max.scans = 10)
pp <- pairwise.probabilities(samples)
print(pp)
print(as.matrix(pp))
conf <- confidence(pp, est)
print(conf)
plot(conf)
plot(conf, data = data)
}
}
% DBD: 2017-05-09
