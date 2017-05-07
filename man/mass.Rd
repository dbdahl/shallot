\name{mass}
\alias{mass}
\alias{discount}
\alias{temperature}
\alias{print.shallot.mass}
\alias{print.shallot.discount}
\alias{print.shallot.temperature}
\title{Mass, Discount, and Temperature Parameters}
\description{
Functions setting mass, discount, and temperature parameters and, in the case of them being random, specifying the parameters of their distribution.}
\usage{
mass(..., fixed = TRUE)
discount(..., fixed = TRUE)
temperature(..., fixed = TRUE)
\method{print}{shallot.mass}(x, ...)
\method{print}{shallot.discount}(x, ...)
\method{print}{shallot.temperature}(x, ...)
}
\arguments{
\item{...}{A number greater than \code{0.0} representing the value of the mass, discount, or temperaure parameters.  Or, in the case of them being random, a vector of two numbers representing either: i. the shape and rate parameters of the gamma distribution for the mass or temperature, or ii. the shape parameters of the beta distribution for the discount.  This argument is currently ignored for the associated print functions.}
\item{x}{An object from the \code{\link{mass}}, \code{\link{discount}}, or \code{\link{temperature}} functions.}
\item{fixed}{If \code{TRUE}, the parameter is fixed.  If \code{FALSE}, the parameter value is samples from either: i. a gamma distribution for the mass or temperature, or ii. a beta distribution for the discount.}
}
\details{
If no parameters are specified, the mass parameter defaults to \code{1.2}, the discount parameter defaults to \code{0.05}, the temperature parameter defaults to \code{3.0}.
If the mass parameter is random, the default shape and rate parameters of the gamma distribution are \code{2.5} and \code{2}, respectively.
If the discount parameter is random, the default shape parameters of the beta distribution are \code{1.0} and \code{1.0}.
If the temperature parameter is random, the default shape and rate parameters of the gamma distribution are \code{2} and \code{0.5}, respectively.
}
\value{
An object of class \code{shallot.mass}, \code{shallot.discount}, or \code{shallot.temperature}.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\examples{
mass()
mass(1.0)
mass(1.4, fixed=FALSE)
mass(0.5, 1, fixed=FALSE)
discount()
discount(0.2)
discount(1, 3, fixed=FALSE)
temperature()
temperature(2)
temperature(2, 4, fixed=FALSE)
}
% DBD: 2017-05-02
