\name{mass}
\alias{mass}
\alias{print.shallot.mass}
\title{Mass}
\description{
  Function storing a given value for the mass and whether the mass is fixed.}
\usage{
mass(..., fixed = TRUE)
\method{print}{shallot.mass}(x, ...)
}
\arguments{
\item{...}{A number greater than but not equal to \code{0.0} representing the mass, OR a vector of two numbers representing the shape and rate parameters of the Gamma function.  As used in \code{\link{print.shallot.mass}}, an unused argument.}
\item{x}{An object of class \code{shallot.mass}.}
\item{fixed}{A boolean object that can equal either \code{TRUE} or \code{FALSE}, \code{TRUE} being the default.  If \code{TRUE}, then the function fixes the mass at the given input value of x.  If \code{FALSE}, then the mass is sampled at random from the Gamma function determined by the given shape and rate parameters.}
}
\details{
With the temperature and decay method both specified, the next variable to consider is the mass.  The mass variable is a location parameter for the Normal distribution controlling the entire shallot project.  This variable controls the similarity of the distributions.  That is, the bigger the mass, the closer together different distributions across different mass values are.  For example, on the Ewens Pitman Attraction method, distributions with mass equal to \code{64} are a lot closer to the distributions with mass equal to \code{80} than they are to those distributions that have mass equal to \code{48}.  On the DDCRP method, however, the temperature is the main variable controlling this property, and the mass variable can only shift the distributions slightly closer to the right.

Despite the above, the mass, as it increases, does cause the distributions to become wider until the mass reaches the bottom of the parabola.  The maximum values of the distributions form a loose parabolic curve.  This is much more noticeable in the DDCRP method when the temperature variable is controlling the number of subsets, but the same thing is true in the other distribution methods when the mass variable is in control.  After reaching the bottom of the parabola, as the mass gets even bigger beyond this point, the distributions become tighter again.

The mass variable has the Gamma distribution with default parameters shape equal to \code{2.5} and rate equal to \code{2}.  The mass variable itself defaults to \code{1.2}.  If no input is specified, then the function will fix the mass at the default value of \code{1.2}.  If an input for the mass value is specified, then the mass function sets that number as the mass value instead of \code{1.2}.  If the input value is less than or equal to \code{0.0}, then an error message will be displayed.  If fixed is set to \code{TRUE}, then the function fixes the mass at the given input value.  If fixed is changed to \code{FALSE}, then regardless of whether or not an input value is specified for the mass, the function will automatically sample a mass value from the Gamma function on the given shape and rate parameters.

If, instead of an input value for the mass, a two-number vector is specified as the input, then the function will set the first of the two numbers as the shape parameter and the other as the rate parameter.  If fixed is changed to \code{FALSE}, then the function will sample the mass value from the Gamma function determined by the new values for the shape and rate parameters.  If fixed is left at \code{TRUE}, however, then an error message will be displayed.  An error message will also be displayed if the vector contains more than two numbers.

The \code{\link{print.shallot.mass}} function is an internal function that automatically displays a concise output for the mass function rather than listing the full results.  If the mass is fixed, then this function will print \dQuote{mass fixed at } and the given mass.  If the mass is not fixed, then this function will print \dQuote{mass sampled from Gamma(shape = } and the shape parameter, followed by \dQuote{, rate = } and the rate parameter.
}
\value{
A data frame of class \code{shallot.mass} containing the established values for the mass, shape, rate, and whether the mass is fixed, along with other internal commands referencing the jobs to be called.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{discount}}
\code{\link{decay}}
}
\examples{
\dontrun{
## Demonstrate mass with mass fixed at given mass input.
mass(1.0, fixed = TRUE)
mass(1.0, fixed = TRUE)$value
mass(1.0, fixed = TRUE)$shape
mass(1.0, fixed = TRUE)$rate
mass(1.0, fixed = TRUE)$fixed

## Demonstrate mass with given mass input, but fixed set to FALSE.
mass(1.0, fixed = FALSE)
mass(1.0, fixed = FALSE)$value
mass(1.0, fixed = FALSE)$shape
mass(1.0, fixed = FALSE)$rate
mass(1.0, fixed = FALSE)$fixed

## Demonstrate mass with given shape and rate input, fixed set to FALSE.
mass(c(0.5, 1), fixed = FALSE)
mass(c(0.5, 1), fixed = FALSE)$value
mass(c(0.5, 1), fixed = FALSE)$shape
mass(c(0.5, 1), fixed = FALSE)$rate
mass(c(0.5, 1), fixed = FALSE)$fixed
}
}

