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
\code{\link{sample.partitions}},
\code{\link{process.partitions}},
\code{\link{estimate.partition}}
}
\examples{
\dontrun{
example(shallot)
}
}
