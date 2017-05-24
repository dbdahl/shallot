\name{decay}
\alias{decay}
\alias{decay.reciprocal}
\alias{decay.exponential}
\alias{decay.subtraction}
\alias{print.shallot.decay}
\title{Decay}
\description{
These functions specify the decay to map distances to attractions.}
\usage{
decay.reciprocal(temperature, distance)
decay.exponential(temperature, distance)
decay.subtraction(temperature, distance, multiplier = 1.01)
\method{print}{shallot.decay}(x, ...)
}
\arguments{
\item{temperature}{An object of class \code{shallot.temperature}.}
\item{distance}{An object of class \code{dist}.}
\item{multiplier}{An scalar greater than \code{1.0} to ensure that attractions from \code{decay.subtraction} are finite.}
\item{x}{An object of class \code{shallot.decay}.}
\item{...}{Currently ignored.}
}
\details{
There are currently three choices for decay functions: reciprocal, exponential, and subtraction.

The reciprocal decay maps a distance \code{d} to an attraction \code{a} as follows: \code{a = 1/d^t}, where \code{t} is the temperature.

The exponential decay maps a distance \code{d} to an attraction \code{a} as follows: \code{a = exp(-t*d)}, where \code{t} is the temperature.

The subtract decay maps a distance \code{d} to an attraction \code{a} as follows: \code{a = (m-d)^t}, where \code{t} is the temperature and \code{m} is the maximum distance in \code{distance} multiplied by the supplied \var{multiplier}.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link[stats]{dist}},
\code{\link{temperature}},
\code{\link{attraction}}
}
\references{
\href{https://doi.org/10.1080/01621459.2016.1165103}{Dahl, D. B., Day, R., and Tsai, J. (2017), Random Partition Distribution Indexed by Pairwise Information, \emph{Journal of the American Statistical Association}, accepted. <DOI:10.1080/01621459.2016.1165103>}
}
\examples{
temp <- temperature(1.0)
distance <- dist(scale(USArrests))
decay1 <- decay.reciprocal(temp,distance)
decay2 <- decay.exponential(temp,distance)
decay3 <- decay.subtraction(temp,distance)
}
