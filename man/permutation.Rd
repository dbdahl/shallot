\name{permutation}
\alias{permutation}
\alias{print.shallot.permutation}
\title{Permutation}
\description{
Function defining a permutation for subsequent use.}
\usage{
permutation(..., n.items = NULL, fixed = TRUE)
\method{print}{shallot.permutation}(x, ...)
}
\arguments{
\item{...}{For the fuction \code{\link{permutation}}, a permutation of the integers \code{1}, \code{2},... \code{n}, where \code{n} is the length of the vector.  For the function \code{\link{print.shallot.permutation}}, this is ignored.}
\item{n.items}{An optional argument provided instead of \code{...} to request a random partition.  The argument \code{fixed} must be \code{FALSE}.}
\item{x}{An object of class \code{shallot.permutation}.}
\item{fixed}{Should the permutation be fixed?}
}
\details{
A valid permutation of length \code{n} is an integer vector of length \code{n} containing each integer \code{1}, \code{2},... \code{n} only once.
}
\value{
An object of class \code{shallot.permutation}.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{attraction}}
}
\examples{
\dontrun{
## Demonstrate permutation.
permutation(c(3, 1, 2, 5, 4))
permutation(c(3, 1, 2, 5, 4), fixed=FALSE)
permutation(n.items=5, fixed=FALSE)
}
}
% DBD: 2017-05-02
