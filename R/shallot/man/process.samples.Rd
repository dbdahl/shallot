\name{process.samples}
\alias{process.samples}
\title{Process Sampled Partitions}
\description{
This function extracts the partitions from the results of the \code{\link{sample.partitions} function.}
}
\usage{
process.samples(x)
}
\arguments{
\item{x}{An object from the \code{\link{sample.partitions}} function.}
}
\details{
This function extracts the sampled partitions from the results of the \code{\link{sample.partitions} function.}
}
\value{
A list containing a matrix of cluster labels in which each row
represents a clusterings.  The list also
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
