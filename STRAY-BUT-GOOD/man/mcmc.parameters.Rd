\name{mcmc.parameters}
\alias{mcmc.parameters}
\title{Setting MCMC Parameters}
\description{
  Function that locks some variable parameters into place.}
\usage{
mcmc.parameters(log.density = NULL, sample = NULL, mass.rw.sd = 0.5,
discount.rw.sd = 0.1, permutation.grab.size = 10, temperature.rw.sd = 0.5,
n.iterations.per.sample = 1)
}
\arguments{
\item{log.density}{A function that returns the logorithm of the sampling density function.  Defaults to \code{NULL}.}
\item{sample}{A function that returns a sample from the centering distribution.  Defaults to \code{NULL}.}
\item{mass.rw.sd}{A number containing the specified mass value.  Defaults to \code{0.5}.}
\item{discount.rw.sd}{A number containing the specified discount value.  Defaults to \code{0.1}.}
\item{permutation.grab.size}{An integer representing the number of rows in the dataset.  Defaults to \code{10}.}
\item{temperature.rw.sd}{A number containing the specified temperature value.  Defaults to \code{0.5}.}
\item{n.iterations.per.sample}{An integer specifying the number of iterations per sample.  Defaults to \code{1}.}
}
\details{
This function gives the user the ability to fine-tune the parameters a little bit.  The name of the function is derived from the fact that it relies on a Markov Chain Monte Carlo method.  This function takes the above input, compiles it, and stores it for future use.  The function also creates an internal variable called sampling.model and adjusts its value depending on whether \var{log.density} is \code{NULL}.
}
\value{
An object of class \code{shallot.mcmc.parameters} storing the various input values listed above and sampling.model.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{collect}}
}
\examples{
\dontrun{
## Demonstrate mcmc.parameters.
mcmc.parameters(log.density = NULL, sample = NULL, mass.rw.sd = 0.5,
    discount.rw.sd = 0.1, permutation.grab.size = 10, temperature.rw.sd =
    0.5, n.iterations.per.sample = 1)
}
}

