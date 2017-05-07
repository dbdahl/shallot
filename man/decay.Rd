\name{decay}
\alias{decay}
\alias{decay.reciprocal}
\alias{decay.exponential}
\alias{decay.subtraction}
\alias{print.shallot.decay}
\title{Decay}
\description{
  Function storing the selected decay method at the given temperature for future use.}
\usage{
decay.reciprocal(..., distance, fixed = TRUE)
decay.exponential(..., distance, fixed = TRUE)
decay.subtraction(..., distance, multiplier = 1.01, fixed = TRUE)
\method{print}{shallot.decay}(x, ...)
}
\arguments{
\item{...}{A numeric vector of length one giving the temperature value, OR a two numeric vectors of length one specifying the shape and rate parameters for the Gamma function that will sample the temperature value at random.  As used in the print functions, an unused argument.}
\item{fixed}{A boolean object that indicates whether the given temperature is fixed or not.  When the temperature is fixed without a temperature specified, the function will default to temperature being \code{3}.  When a two-number vector is specified instead of a single number, \var{fixed} is required to be \code{FALSE}.}
\item{multiplier}{An additional factor that affects the decay rate.  The default is \code{1.01}.}
\item{distance}{An object of class \code{shallot.distance} which will be used to ensure that, in the class of \code{\var{fixed}=TRUE}, the temperatures do not lead to problems with machine precision.}
\item{x}{As used in \code{\link{print.shallot.temperature}}, an object of class \code{shallot.temperature}.  As used in \code{\link{print.shallot.decay}}, an object of class \code{shallot.decay}.}
}
\details{
The shallot program is designed to determine how similar two alternate distributions of the dataset truly are by allowing the user to modify several variables along the way.  The program is built on a Dirichlet process in which each variable has a different distribution.  In fact, the entire shallot program actually generates a Normal distribution of the number of subsets against the proprotion of partitions that have each particular number of subsets.  In this distribution, the temperature variable is a  parameter that dictates how wide or tight the distribution is, although admittedly there is relatively little variation graphically.  However, the temperature does greatly affect the adjusted Rand index, and it has enormous control not only over the concentration, but also over the number of subsets, when the DDCRP distribution method is being used.  Hence, on the DDCRP method, it is a location parameter.  On DDCRP, the smaller and larger temperature values yield distributions that are more tightly concentrated, while the distribution at a temperature value of approximately \code{3} yields the widest distribution.  On the other distribution methods, the mass variable has control over this property, but on DDCRP, it is the temperature.

The temperature variable has the Gamma distribution with default parameters shape of \code{2} and rate of \code{0.5}.  The default temperature is \code{3}.  The temperature function allows the user to establish either the temperature of the distribution or the shape and rate parameters of the Gamma distribution.  This function is an internal function only and cannot be used directly by the user.  The parameters are therefore established in the decay functions described below, and the temperature function is called internally.  It accepts as input either a single number intended for the temperature itself or a two-number vector containing the intended values for the shape and rate parameters of the Gamma function from which the temperature is to be sampled.  When the temperature is specified, the shape and rate parameters will remain at their default values, and vice-versa.  When no numbers are specified, all three internal variables will remain at their default values.  When a numeric vector is specified, \var{fixed} is required to be \code{FALSE}.

The decay functions, having called the temperature function internally to establish the parameters, then indicate how the distribution tails are to decay.  The three choices are exponential, subtraction, and reciprocal.  The temperature indicates by how much they are to decay.  The three decay functions simply accept the given input values of the temperature/parameters and whether the temperature is fixed or not.  Each function then passes these values on to the \code{decay.generic} function, along with specifying the values for the multiplier and the decay type (\dQuote{reciprocal},\dQuote{exponential}, or \dQuote{subtraction}).

The \code{decay.generic} function takes the above inputs, lists them, and calls the temperature function described above to process the temperature prior to listing.
}
\note{
When each decay function is written out with all of its arguments, it is coded internally as in the following example:
decay.exponential(\code{9.0}, \var{fixed} = \code{TRUE}) decay.generic(..., \var{multiplier} = \code{1.01}, type = \dQuote{exponential}, \var{fixed} = \code{fixed})
The decay functions are not coded as such in the usage section because the package program read this format as an error.  However, this is the way it is written in shallot.R.

Also, there is hardly any difference in the adjusted Rand index between the three decay methods.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{mass}}
\code{\link{discount}}
\code{\link{attraction}}
}
\examples{
\dontrun{
## Demonstrate decay.exponential.  (The same for decay.reciprocal and decay.subtraction functions.)
decay.exponential(9.0, fixed = TRUE)
decay.exponential(9.0, fixed = TRUE)$type
decay.exponential(9.0, fixed = TRUE)$multiplier
decay.exponential(9.0, fixed = TRUE)$temperature$value
decay.exponential(9.0, fixed = TRUE)$temperature$shape
decay.exponential(9.0, fixed = TRUE)$temperature$rate
decay.exponential(9.0, fixed = TRUE)$temperature$fixed
decay.exponential(c(2, 3), fixed = FALSE)
decay.exponential(fixed = TRUE)
decay.exponential(fixed = FALSE)
}
}

