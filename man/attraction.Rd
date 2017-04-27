\name{attraction}
\alias{attraction}
\alias{print.shallot.attraction}
\alias{as.matrix.shallot.attraction}
\title{Attraction}
\description{
  Function storing the selected decay method at the given temperature for future use.}
\usage{
attraction(distance, permutation, decay)
\method{print}{shallot.attraction}(x, ...)
\method{as.matrix}{shallot.attraction}(x, ...)
}
\arguments{
\item{distance}{A data frame of class \code{shallot.distance} containing the number of observations, the name of each observation, and the established values for the various distances between the observations of the dataset (although the values are not immediately accessible by the user at this stage), along with some internal commands referencing the jobs to be called.}
\item{permutation}{A data frame of class \code{shallot.permutation} containing the value of n.items, whether the permutation is fixed, and the permutation itself (if given; again, this is not accessible by the user at this stage), along with some internal commands referencing the jobs to be called.}
\item{decay}{A data frame of class \code{shallot.decay} containing the values of the temperature, shape, and rate, whether the permutation is fixed, and the decay method type and the multiplier.}
\item{x}{An object of class \code{shallot.attraction} containing the stored variable values described above.}
\item{...}{An additional argument, not used.}
}
\details{
The \code{\link{attraction}} function brings all of the values specified in the arguments section above (specifically, \var{distance}, \var{permutation}, and \var{decay}) together into one data frame for further processing.  At this stage, for the \code{shallot.permutation} function, no permutation is specified, n.items is set to the number of observations in the dataset, and fixed must equal \code{FALSE}.  For the most part, each function in this package simply compiles the work of the previous functions in the process rather than actually doing anything significant.  Yet multiple functions are used instead of compiling everything in one large function to allow the user greater flexibility and convenience in adjusting the values for the variables.  When the \code{\link{attraction}} function is called, for example, the variables that have been adjusted thus far include the temperature, whether the temperature is fixed, and the type of decay method selected.  The function also calls several jobs as internal variable values but does not actually execute these jobs at this time.

The \code{\link{print.shallot.attraction}} function is an internal function similar to the \code{\link{print.shallot.temperature}} function.  It prints a statement that restates the number of observations, the decay method selected, and the value of the temperature (or if the temperature is not fixed, the values of the shape and rate parameters for the Gamma distribution).

The \code{\link{as.matrix.shallot.attraction}} function takes as its input the distance attribute of the \code{\link{attraction}} function result and displays these distances in matrix form.
}
\value{
A data frame of class \code{shallot.attraction} that contains the arguments listed above.  For more details on the contents of these output values, please see the links below.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
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
attraction(distance(as.matrix(distance)), permutation(n.items = nrow(data) ,fixed = FALSE),
    decay.exponential(9.0, fixed = TRUE))
attraction(distance(as.matrix(distance)), permutation(n.items = nrow(data), fixed = FALSE),
    decay.subtraction(9.0, fixed = TRUE))
attraction(distance(as.matrix(distance)), permutation(n.items = nrow(data), fixed = FALSE),
    decay.reciprocal(9.0, fixed = TRUE))
}
}

