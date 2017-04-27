\name{q}
\alias{q}
\alias{random.q}
\alias{probability.q}
\alias{average.q}
\alias{variance.q}
\title{Q}
\description{
  Functions whose titles end as \dQuote{.q} and handle several lesser tasks.}
\usage{
random.q(x, n.samples)
probability.q(x, n.subsets)
average.q(x)
variance.q(x)
}
\arguments{
\item{x}{An object of class \code{shallot.distribution}, varying between any of the six distribution methods provided.}
\item{n.samples}{An integer containing the number of samples.  The number of samples is usually at least \code{500}.}
\item{n.subsets}{An integer containing the number of subsets.  This corresponds to the number of options within the truth column.}
\item{...}{An unused argument.}
}
\details{
The \code{\link{random.q}} function returns an array of the number of subsets for each of the specified samples.  The number of subsets can vary to much higher than the number of maxSubsets in the truth column, but must be at least \code{1}.  As the name suggests, the number of subsets for each sample is chosen at random each time \code{\link{random.q}} is run.

The \code{\link{probability.q}} function returns the probability that a dataset with a given distribution method will obtain the given input number of subsets.  The mass and the discount (if present) have to be fixed or an error message will appear.

The \code{\link{average.q}} function returns the average number of subsets for the given distribution method.  The mass and the discount (if present) have to be fixed or an error message will appear.

The \code{\link{variance.q}} function returns the variance of the number of subsets available for the given distribution method.  The mass and the discount (if present) have to be fixed or an error message will appear.  However, at present, this function is really only designed to work with the ewens and the ewensAttraction methods and will redirect the user to the \code{\link{random.q}} function for everything else.
}
\value{
The \code{\link{random.q}} function returns an array of integers representing the number of subsets for each of the specified samples.  The other three functions each return a number that, depending on the function, represents the probability of obtaining a certain number of subsets, the average number of subsets, or the variance of the subsets.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{distribution}}
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

## Demonstrate random.q.
random.q(distribution, 500)

## Demonstrate probability.q.
probability.q(distribution, 2)

## Demonstrate average.q.
average.q(distribution)

## Demonstrate variance.q.
distribution <- ewens(mass, 132)
variance.q(distribution)
}
}

