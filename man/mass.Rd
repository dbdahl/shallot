\name{mass}
\alias{mass}
\alias{discount}
\alias{print.shallot.mass}
\alias{print.shallot.discount}
\title{Mass and Discount Parameters}
\description{
Functions setting mass and discount parameters and, in the case of them being random, specifying the parameters of their distribution.}
\usage{
mass(..., fixed = TRUE)
discount(..., fixed = TRUE)
\method{print}{shallot.mass}(x, ...)
\method{discount}{shallot.mass}(x, ...)
}
\arguments{
\item{...}{A number greater than \code{0.0} representing the value of the mass or discount parameters.  Or, in the case of them being random, a vector of two numbers representing either: i. the shape and rate parameters of the gamma distribution for the mass, or ii. the shape parameters of the beta distribution for the discount.  This argument is currently ignored for the \code{\link{print.shallot.mass}} and \code{\link{print.shallot.discount}}, an unused argument.}
\item{x}{An object from the \code{\link{mass}} or \code{\link{discount}} functions.}
\item{fixed}{If \code{TRUE}, the parameter is fixed.  If \code{FALSE}, the parameter value is samples from either: i. a gamma distribution for the mass, or ii. a beta distribution for the discount.}
}
\details{
If no parameters are specified, the mass parameter defaults to \code{1.2} and the discount parameter defaults to \code{0.05}.
If the mass parameter is random, the default parameters shape and rate of the gamma distribution are \code{2.5} and \code{2}, respectively.
If the discount parameter is random, the default parameters shape parameters of the beta distribution are \code{1.0} and \code{1.0}.

  The mass variable itself defaults to \code{1.2}.  If no input is specified, then the function will fix the mass at the default value of \code{1.2}.  If an input for the mass value is specified, then the mass function sets that number as the mass value instead of \code{1.2}.  If the input value is less than or equal to \code{0.0}, then an error message will be displayed.  If fixed is set to \code{TRUE}, then the function fixes the mass at the given input value.  If fixed is changed to \code{FALSE}, then regardless of whether or not an input value is specified for the mass, the function will automatically sample a mass value from the Gamma function on the given shape and rate parameters.
}
\value{
A object of class \code{shallot.mass} or \code{shallot.discount}.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\examples{
mass()
mass(1.0, fixed=TRUE)
mass(1.4, fixed=FALSE)
mass(0.5, 1, fixed=FALSE)
discount()
discount(0.2, fixed=TRUE)
discount(1, 3, fixed=FALSE)
}

