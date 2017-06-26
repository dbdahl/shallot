\name{process.partitions}
\alias{process.partitions}
\title{Process Sampled Partitions}
\description{
This function extracts the partitions from the results of the \code{\link{sample.partitions} function.}
}
\usage{
process.partitions(x, sample.parameter = NULL, as.matrix = TRUE)
}
\arguments{
\item{x}{An object from the \code{\link{sample.partitions}} function.}
\item{sample.parameter}{If not \code{NULL}, a function taking no arguments and returning a random sample from the centering distribution.}
\item{as.matrix}{Should the results be returned as a matrix of cluster labels or list of partitions?}
}
\details{
This function extracts the sampled partitions from the results of the \code{\link{sample.partitions} function.}
}
\value{
A list containing either: i. a matrix of cluster labels in which each row
represents a clusterings, or ii. a list containing partitions.  The list also
contains sampled model parameters if \code{sample.parameter} is not
\code{NULL}.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{sample.partitions}}
}
\examples{
\dontrun{
example(shallot)
}
}
