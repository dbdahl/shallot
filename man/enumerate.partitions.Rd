\name{enumerate.partitions}
\alias{enumerate.partitions}
\title{Enumerate Partitions}
\description{
This function extracts the partitions from the results of the \code{\link{sample.partitions} function.}
}
\usage{
enumerate.partitions(n.items, as.matrix = TRUE)
}
\arguments{
\item{n.items}{An integer given then number of items to partition.}
\item{as.matrix}{Should the results be returned as a matrix of cluster labels or list of partitions?}
}
\details{
This function returns an enumeration of the partition of \code{n.items} items.
}
\value{
A list containing either: i. a matrix of cluster labels in which each row
represents a clusterings, or ii. a list containing partitions.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{process.samples}}
}
\examples{
\dontrun{
example(shallot)
}
}
