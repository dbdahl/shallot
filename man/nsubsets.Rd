\name{nsubsets}
\alias{nsubsets.random}
\alias{nsubsets.probability}
\alias{nsubsets.average}
\alias{nsubsets.variance}
\title{Q}
\description{
  Functions whose titles end as \dQuote{.q} and handle several lesser tasks.}
\usage{
nsubsets.random(x, n.samples)
nsubsets.probability(x, n.subsets)
nsubsets.average(x)
nsubsets.variance(x)
}
\arguments{
\item{x}{An object of class \code{shallot.distribution}, varying between any of the six distribution methods provided.}
\item{n.samples}{An integer containing the number of samples.  The number of samples is usually at least \code{500}.}
\item{n.subsets}{An integer containing the number of subsets.  This corresponds to the number of options within the truth column.}
\item{...}{An unused argument.}
}
\details{
The \code{\link{nsubsets.random}} function returns an array of the number of subsets for each of the specified samples.  The number of subsets can vary to much higher than the number of maxSubsets in the truth column, but must be at least \code{1}.  As the name suggests, the number of subsets for each sample is chosen at random each time \code{\link{nsubsets.random}} is run.

The \code{\link{nsubsets.probability}} function returns the probability that a dataset with a given distribution method will obtain the given input number of subsets.  The mass and the discount (if present) have to be fixed or an error message will appear.

The \code{\link{nsubsets.average}} function returns the average number of subsets for the given distribution method.  The mass and the discount (if present) have to be fixed or an error message will appear.

The \code{\link{nsubsets.variance}} function returns the variance of the number of subsets available for the given distribution method.  The mass and the discount (if present) have to be fixed or an error message will appear.  However, at present, this function is really only designed to work with the ewens and the ewensAttraction methods and will redirect the user to the \code{\link{nsubsets.random}} function for everything else.
}
\value{
The \code{\link{nsubsets.random}} function returns an array of integers representing the number of subsets for each of the specified samples.  The other three functions each return a number that, depending on the function, represents the probability of obtaining a certain number of subsets, the average number of subsets, or the variance of the subsets.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{distribution}}
}
\references{
Dahl, D. B., Day, R., and Tsai, J. (2017), Random Partition Distribution Indexed by Pairwise Information, \emph{Journal of the American Statistical Association}, accepted.
}
\examples{
\dontrun{
## Load example dataset.
data <- Theoph[, c("Wt", "Dose", "Time", "conc", "Subject")]
truth <- as.numeric(levels(data[, ncol(data)]))[data[, ncol(data)]]
distance <- dist(scale(data[, -ncol(data)]))

## Prepare attraction.
attract <- attraction(distance,
    permutation(n.items=nrow(data), fixed = FALSE), decay.exponential(9.0, fixed = TRUE))
mass <- mass(1.0, fixed = TRUE)
discount <- discount(0.1, fixed = TRUE)
distribution <- ewens.pitman.attraction(mass, discount, attract)

## Demonstrate nsubsets.random.
nsubsets.random(distribution, 500)

## Demonstrate nsubsets.probability.
nsubsets.probability(distribution, 2)

## Demonstrate nsubsets.average.
nsubsets.average(distribution)

## Demonstrate nsubsets.variance.
distribution <- ewens(mass, 132)
nsubsets.variance(distribution)
}
}

