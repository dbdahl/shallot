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
\item{...}{Additional arguments (that are currently ignored).}
}
\value{
An object of class \code{shallot.attraction}.  
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{ddcrp}}}
\code{\link{decay}}
\code{\link{ewens.attraction}}
\code{\link{ewens.pitman.attraction}}
\code{\link{permutation}}
}
% DBD: 2017-05-09
