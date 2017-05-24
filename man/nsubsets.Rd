\name{nsubsets}
\alias{nsubsets.random}
\alias{nsubsets.probability}
\alias{nsubsets.average}
\alias{nsubsets.variance}
\title{Number of Subsets}
\description{
These functions either sample the number of subsets for supported partition distributions or computes probabilities, means, and variances of these distributions.
}
\usage{
nsubsets.random(x, n.samples)
nsubsets.probability(x, n.subsets)
nsubsets.average(x)
nsubsets.variance(x)
}
\arguments{
\item{x}{An object of class \code{shallot.distribution}.}
\item{n.samples}{An integer containing the number of samples.}
\item{n.subsets}{An integer containing the number of subsets.}
\item{...}{Currently ignored.}
}
\value{
The \code{\link{nsubsets.random}} function returns a vector of random samples of the number of subsets in the distribution \var{x}.

The \code{\link{nsubsets.probability}} function returns the probability that the number of subsets is \var{n.subsets} in the distribution \var{x}.  Depending on the number of items and the value of \var{n.subsets}, this function can be computationally intensive.

The \code{\link{nsubsets.average}} and \code{\link{nsubsets.variance}} functions return the mean and variances, respectively, of the number of subsets in the distribution \var{x}.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{partition.distribution}}
}
\references{
\href{https://doi.org/10.1080/01621459.2016.1165103}{Dahl, D. B., Day, R., and Tsai, J. (2017), Random Partition Distribution Indexed by Pairwise Information, \emph{Journal of the American Statistical Association}, accepted. <DOI:10.1080/01621459.2016.1165103>}
}
\examples{
\donttest{
pd <- ewens.pitman.attraction(
        mass(1),
        discount(0.05),
        attraction(permutation(n.items=50,fixed=FALSE),
          decay.exponential(temperature(1.0),dist(scale(USArrests)))))
mean(nsubsets.random(pd,1000))
nsubsets.average(pd)

pde <- ewens(mass(1),50)
nsubsets.variance(pde)
nsubsets.probability(pde,4)
}
}
