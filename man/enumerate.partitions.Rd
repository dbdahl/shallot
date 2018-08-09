\name{enumerate.partitions}
\alias{enumerate.partitions}
\title{Enumerate Partitions}
\description{
This function enumerates all possible partitions for a given number of items.}
\usage{
enumerate.partitions(n.items)
}
\arguments{
\item{n.items}{An integer given then number of items to partition.}
}
\details{
This function returns an enumeration of the partition of \code{n.items} items.
}
\value{
A matrix of cluster labels in which each row represents a clusterings.
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
