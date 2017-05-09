\name{attraction}
\alias{attraction}
\alias{print.shallot.attraction}
\alias{as.matrix.shallot.attraction}
\title{Attraction}
\description{
This function creates an attraction from a distance, permutation, and decay in preparation for use in the \code{\link{ewens.pitman.attraction}} function.
For details on each of these arguments, please see the links below.
}
\usage{
attraction(distance, permutation, decay)
\method{print}{shallot.attraction}(x, ...)
\method{as.matrix}{shallot.attraction}(x, ...)
}
\arguments{
\item{distance}{An object of class \code{dist} encoding the distances between the items.}
\item{permutation}{An object of class \code{shallot.permutation} encoding the permutation of the items.}
\item{decay}{An object of class \code{shallot.decay} detailing the transformation from distances to attractions.}
\item{x}{An object of class \code{shallot.attraction}.}
\item{...}{Additional arguments (that are currently ignored).}
}
\value{
An object of class \code{shallot.attraction}.  
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{stats::dist}}
\code{\link{permutation}}
\code{\link{decay}}
}
\examples{
\dontrun{
## Load example dataset.
data <- Theoph[, c("Wt", "Dose", "Time", "conc", "Subject")]
truth <- as.numeric(levels(data[, ncol(data)]))[data[, ncol(data)]]
distance <- dist(scale(data[, -ncol(data)]))

## Demonstrate attraction.
attraction(distance, permutation(n.items = nrow(data) ,fixed = FALSE),
    decay.exponential(9.0, fixed = TRUE))
attraction(distance, permutation(n.items = nrow(data), fixed = FALSE),
    decay.subtraction(9.0, distance=distance, fixed = TRUE))
attraction(distance, permutation(n.items = nrow(data), fixed = FALSE),
    decay.reciprocal(9.0, fixed = TRUE))
}
}

