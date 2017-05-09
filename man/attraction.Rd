\name{attraction}
\alias{attraction}
\alias{print.shallot.attraction}
\alias{as.matrix.shallot.attraction}
\title{Attraction}
\description{
This function creates an attraction from a permutation and a decay in preparation for use in the \code{\link{ewens.attraction}}, \code{\link{ewens.pitman.attraction}}, and \code{\link{ddcrp}} functions.
For details on each of these arguments, please see the links below.
}
\usage{
attraction(permutation, decay)
\method{print}{shallot.attraction}(x, ...)
\method{as.matrix}{shallot.attraction}(x, ...)
}
\arguments{
\item{permutation}{An object of class \code{shallot.permutation} encoding the permutation of the items.}
\item{decay}{An object of class \code{shallot.decay} detailing the transformation from distances to attractions.}
\item{x}{An object of class \code{shallot.attraction}.}
\item{...}{Currently ignored.}
}
\value{
An object of class \code{shallot.attraction}.  
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{ddcrp}},
\code{\link{decay}},
\code{\link{ewens.attraction}},
\code{\link{ewens.pitman.attraction}},
\code{\link{permutation}}
}
\references{
\href{https://doi.org/10.1080/01621459.2016.1165103}{Dahl, D. B., Day, R., and Tsai, J. (2017), Random Partition Distribution Indexed by Pairwise Information, \emph{Journal of the American Statistical Association}, accepted.}
}
\examples{
permutation <- permutation(n.items=50,fixed=FALSE)
decay <- decay.exponential(temperature(1.0),dist(scale(USArrests)))
attraction(permutation, decay)
}
% DBD: 2017-05-09
