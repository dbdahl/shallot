\name{permutation}
\alias{permutation}
\alias{print.shallot.permutation}
\title{Permutation}
\description{
  Function storing the given permutation for future use.}
\usage{
permutation(..., n.items = NULL, fixed = TRUE)
\method{print}{shallot.permutation}(x, ...)
}
\arguments{
\item{...}{A list containing the intended permutation to be stored and used later.  The list must not include anything other than consecutive integers starting at \code{1}, or an error message will appear.  However, these integers do not have to be in consecutive order.  For example, if the integers were \code{1}, \code{2}, and \code{3}, they can be inputted as \code{c(2,3,1)}.  For the shallot process, however, this particular argument usually remains blank.  As used in \code{\link{print.shallot.permutation}}, an unused argument.}
\item{x}{An object of class \code{shallot.permutation}.}
\item{n.items}{An integer stating the number of items being considered.  The default is \code{NULL}, but in the shallot process, n.items is usually set to be the number of observations in the dataset.  If the value assigned to n.items is different from the number of values in the permutation argument described above, then an error message appears.  When the permutation argument described above is left blank, n.items is required to have a specific value.  If it is \code{NULL}, then an error message will appear.}
\item{fixed}{A boolean object that is either \code{TRUE} or \code{FALSE}.  The default is \code{TRUE}, but if the permutation argument described above is left blank, then fixed is required to be \code{FALSE} or an error message will appear.}
}
\details{
The \code{\link{permutation}} function takes as its input a list of numbers representing an intended permutation.  The function will subtract \code{1} from each value in the list to start the series at \code{0} instead of at \code{1}.  However, for the purposes of the shallot process, no permutation is specified.  When this is the case, the function selects a permutation at random from the given dataset.  In order for this to work properly, n.items and fixed cannot remain at their default values.  Fixed has to be set to \code{FALSE} and n.items has to equal the number of observations in the dataset, or error messages will appear.  The permutation is sampled from the uniform distribution, and the internal coding does not specify any default values for the parameters.

The \code{\link{print.shallot.permutation}} function simply prints the phrase \dQuote{permutation of }, n.items, \dQuote{ items fixed at }, and the permutation list.  If no permutation is specified, then it prints \dQuote{permutation of }, n.items, \dQuote{ items sampled uniformly}.
}
\value{
A data frame of class \dQuote{shallot.permutation} containing the value of n.items, whether the permutation is fixed, and the permutation itself (if given; the permutation itself is not accessible by the user at this stage), along with some internal commands referencing the jobs to be called.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{attraction}}
}
\examples{
\dontrun{
## Demonstrate permutation.
permutation(c(1, 2, 3), n.items = 3, fixed = TRUE)
permutation(c(1, 2, 3), n.items = 3, fixed = TRUE)$n.items
permutation(c(1, 2, 3), n.items = 3, fixed = TRUE)$fixed
permutation(c(1, 2, 3) ,n.items = 3, fixed = FALSE)
permutation(c(1, 2, 3), n.items = 3, fixed = FALSE)$n.items
permutation(c(1, 2, 3), n.items = 3, fixed = FALSE)$fixed

## Demonstrate permutation with no permutation specified and with n.items set to a specific value.
permutation(n.items = 3, fixed = FALSE)
}
}

