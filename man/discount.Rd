\name{discount}
\alias{discount}
\alias{print.shallot.discount}
\title{Discount}
\description{
  Function storing a given value for the discount and whether the discount is fixed.}
\usage{
discount(..., fixed = TRUE)
\method{print}{shallot.discount}(x, ...)
}
\arguments{
\item{...}{A number greater than or equal to \code{0.0} and less than but not equal to \code{1.0} representing the discount, OR a vector of two numbers representing the shape parameters of the Beta function.  As used in \code{\link{print.shallot.discount}}, an unused argument.}
\item{x}{An object of class \code{shallot.discount}.}
\item{fixed}{A boolean object that can equal either \code{TRUE} or \code{FALSE}, \code{TRUE} being the default.  If \code{TRUE}, then the function fixes the mass at the given input value of x.  If \code{FALSE}, then the mass is sampled at random from the Gamma function determined by the given shape and rate parameters.}
}
\details{
The discount variable is another concentration parameter, in that the bigger the discount, the more it concentrates the distributions closer to the right of the parabola described in the \code{\link{mass}} function instructions.  Furthermore, the discount seems to wield a greater amount of control over the parabolic curve than the mass does when on the Ewens Pitman Attraction method.

This function is largely similar to the \code{\link{mass}} function, except that the input may equal \code{0.0} but may not equal \code{1.0} or higher, and the Beta function is sampled instead of the Gamma function when the input is not fixed.

The \code{\link{discount}} function already has default values of \code{0.05} for the discount value and \code{1} for both of the shape parameters.  If no input is specified, then the function will fix the discount at the default value of \code{0.05}.  If an input for the discount value is specified, then the \code{\link{discount}} function sets that number as the discount value instead of \code{0.05}.  If the input value is less than \code{0.0} or greater than or equal to \code{1.0}, then an error message will be displayed.  If fixed is set to \code{TRUE}, then the function fixes the discount at the given input value.  If fixed is changed to \code{FALSE}, then regardless of whether or not an input value is specified for the discount, the function will automatically sample a discount value from the Beta function on the given shape parameters.

If, instead of an input value for the discount, a two-number vector is specified as the input, then the function will set the two numbers as the shape parameters.  If fixed is changed to \code{FALSE}, then the function will sample the discount value from the Beta function determined by the new values for the shape parameters.  If fixed is left at \code{TRUE}, however, then an error message will be displayed.  An error message will also be displayed if the vector contains more than two numbers.

The \code{\link{print.shallot.discount}} function is an internal function that automatically displays a concise output for the \code{\link{discount}} function rather than listing the full results.  If the discount is fixed, then this function will print \dQuote{discount fixed at } and the given discount.  If the discount is not fixed, then this function will print \dQuote{discount sampled from Beta(shape1=} and the first shape parameter, followed by \dQuote{,shape2=} and the other shape parameter.
}
\value{
A data frame of class \code{shallot.discount} containing the established values for the discount, shape1, shape2, and whether the discount is fixed, along with other internal commands referencing the jobs to be called.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{mass}}
\code{\link{decay}}
}
\examples{
\dontrun{
## Demonstrate discount with discount fixed at given discount input.
discount(0.5, fixed = TRUE)
discount(0.5, fixed = TRUE)$value
discount(0.5, fixed = TRUE)$shape1
discount(0.5, fixed = TRUE)$shape2
discount(0.5, fixed = TRUE)$fixed

## Demonstrate discount with given discount input, but fixed set to FALSE.
discount(0.5, fixed = FALSE)
discount(0.5, fixed = FALSE)$value
discount(0.5, fixed = FALSE)$shape1
discount(0.5, fixed = FALSE)$shape2
discount(0.5, fixed = FALSE)$fixed

## Demonstrate discount with given discount input at 1.0.
discount(0.99, fixed = TRUE)

## Demonstrate discount with given shape1 and shape2 input, fixed set to FALSE.
discount(c(0.5, 1), fixed = FALSE)
discount(c(0.5, 1), fixed = FALSE)$value
discount(c(0.5, 1), fixed = FALSE)$shape1
discount(c(0.5, 1), fixed = FALSE)$shape2
discount(c(0.5, 1), fixed = FALSE)$fixed
}
}

