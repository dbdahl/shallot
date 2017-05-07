\name{distribution}
\alias{distribution}
\alias{ewens}
\alias{print.shallot.distribution.ewens}
\alias{ewens.pitman}
\alias{print.shallot.distribution.ewensPitman}
\alias{ewens.attraction}
\alias{print.shallot.distribution.ewensAttraction}
\alias{ewens.pitman.attraction}
\alias{print.shallot.distribution.ewensPitmanAttraction}
\alias{ddcrp}
\alias{print.shallot.distribution.ddcrp}
\title{Distribution}
\description{
  Functions handling the selection of the distribution method to be used.}
\usage{
ewens(mass,n.items,names = paste0("c", 1:n.items))
\method{print}{shallot.distribution.ewens}(x, ...)
ewens.pitman(mass,discount,n.items,names = paste0("c", 1:n.items))
\method{print}{shallot.distribution.ewensPitman}(x, ...)
ewens.attraction(mass, attraction)
\method{print}{shallot.distribution.ewensAttraction}(x, ...)
ewens.pitman.attraction(mass, discount, attraction)
\method{print}{shallot.distribution.ewensPitmanAttraction}(x, ...)
ddcrp(mass, attraction)
\method{print}{shallot.distribution.ddcrp}(x, ...)
}
\arguments{
\item{mass}{An object of class \code{shallot.mass}.  Please see instructions for the \code{\link{mass}} function.}
\item{n.items}{An integer containing the number of rows in the dataset.  Obtained from \code{attraction$n.items}.  Please see instructions for the \code{\link{attraction}} function.}
\item{names}{An array containing the names of the observations.  The names assigned are of the form \dQuote{c1}, \dQuote{c2}, etc.}
\item{discount}{An object of class \code{shallot.discount}.  Please see instructions for the \code{\link{discount}} function.}
\item{attraction}{An object of class \code{shallot.attraction}.  Please see instructions for the \code{\link{attraction}} function.}
\item{x}{An object of class \code{shallot.distribution} with the subclass varying depending on the \code{\link{distribution}} type.}
\item{...}{An unused argument.}
}
\details{
At this stage, the user can now specify which distribution method to use.  There are five to choose from, but two in particular seem to work the best.  The five methods are ddcrp, ewens, ewensPitman, ewensAttraction, and ewensPitmanAttraction.

The \code{\link{ewens}} function takes as its input the \code{shallot.mass} object and the number of rows in the dataset and stores them for subsequent use.

The \code{\link{ewens.pitman}} function does the same, but also includes the \code{shallot.discount} object.

The \code{\link{ewens.attraction}} function does the same.

The \code{\link{ewens.pitman.attraction}} function does the same, except that it also includes the \code{shallot.discount} object.

For the input values for these functions, please see the instructions for the \code{\link{mass}}, \code{\link{discount}}, and \code{\link{attraction}} functions.

The various print functions simply print the given values for the variables.
}
\value{
An object of class \code{shallot.distribution} with an additional subclass that varies depending on the function.  The five options for the subclass are ewens, ewensPitman, ddcrp, ewensAttraction, and ewensPitmanAttraction.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{mass}}
\code{\link{discount}}
\code{\link{attraction}}
\code{\link{sample.partition}}
}
\examples{
\dontrun{
## Load example dataset.
data <- Theoph[, c("Wt", "Dose", "Time", "conc", "Subject")]
truth <- as.numeric(levels(data[, ncol(data)]))[data[, ncol(data)]]
distance <- dist(scale(data[, -ncol(data)]))

## Prepare attraction.
attract <- attraction(distance, permutation(n.items =
    nrow(data), fixed = FALSE), decay.exponential(9.0, fixed = TRUE))
mass <- mass(1.0, fixed = TRUE)
discount <- discount(0.1, fixed = TRUE)
p <- 0.5

## Demonstrate ewens.
ewens(mass, n.items = attract$n.items)

## Demonstrate ewens.pitman.
ewens.pitman(mass, discount, n.items = attract$n.items)

## Demonstrate ewens.attraction.
ewens.attraction(mass, attract)

## Demonstrate ewens.pitman.attraction.
ewens.pitman.attraction(mass, discount, attract)
}
}

