\name{sample.partitions.posterior}
\alias{sample.partitions.posterior}
\title{Sample Partitions from Posterior Distribution of Partition}
\description{
This function samples partitions from the posterior distribution of a partition based on a user-supplied likelihood and the following prior partition distributions:  Ewens, Ewens-Pitman, Ewens attraction, Ewens-Pitman attraction, and ddCRP distributions.
}
\usage{
sample.partitions.posterior(partition, sampling.model, partition.model,
    n.draws, massRWSD = 0.5, discountRWSD = 0.1,
    k = min(length(partition), 25), temperatureRWSD = 0.5,
    progress.bar = interactive())
}
\arguments{
\item{partition}{An object of class \code{shallot.distribution.data}}
\item{sampling.model}{An object of class \code{shallot.distribution.data} obtained from the \code{\link{sampling.model}} function.}
\item{partition.model}{An object of class \code{shallot.distribution} obtained, for example, from the \code{\link{ewens.pitman.attraction}} function.}
\item{n.draws}{An integer representing the desired number of samples.}
\item{massRWSD}{The standard deviation of the random walk proposal for updating the mass parameter.}
\item{discountRWSD}{The standard deviation of the random walk proposal for updating the discount parameter.}
\item{k}{The number of items to shuffle when proposing an update for the permutation.}
\item{temperatureRWSD}{The standard deviation of the random walk proposal for updating the temperature parameter.}
\item{progress.bar}{Should a progress bar be shown while sampling?}
}
\value{
An object of class \code{shallot.samples.raw} which can be subsequently be used in \code{\link{process.samples}}, \code{\link{pairwise.probabilities}}, \code{\link{estimate.partition}},  
}
\note{
If this function is interrupted by the user, the computation engine will be broken and subsequent calls to package functions may fail until a new \R session is started.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{partition.distribution}},
\code{\link{process.samples}},
\code{\link{pairwise.probabilities}},
\code{\link{estimate.partition}}
\code{\link{sample.partitions}}
}
\examples{
\donttest{
mass <- mass(1.0, fixed=TRUE)
discount <- discount(0.05, fixed=TRUE)
distance <- dist(scale(USArrests[1:9,]))
if ( min(distance[upper.tri(distance)],na.rm=TRUE) == 0 )
  stop("Oops, distances must be strictly positive.")

n.items <- attr(distance,"Size")
permutation <- permutation(n.items=n.items, fixed=FALSE)
temperature <- temperature(2,fixed=TRUE)
attraction <- attraction(permutation,decay.exponential(temperature,distance))
partition.distribution <- ewens.pitman.attraction(mass, discount, attraction)

## Model inputs.
data <- c(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
sigma  <- 0.1
mu0    <- 0.0
sigma0 <- 1.0

## Derived values.
s2 <- sigma * sigma
s02 <- sigma0 * sigma0
s02Inv <- 1.0 / s02
c <- -1.0 / (2.0 * s2)

## Sampling model of Neal (JCGS, 2009)
## Function to perform an MCMC update of the parameter.
sample.parameter <- function(indices=c(), parameter=NULL) {
  sum <- sum(data[indices])
  variance <- 1 / (s02Inv + length(indices) / s2)
  mean <- variance * (mu0 / s02 + sum / s2)
  rnorm(1, mean=mean, sd=sqrt(variance))
}

## Function to evaluate the likelihood contribution for an observation.
log.density <- function(i, indices, parameter) {
  resid <- data[i] - parameter
  c * resid * resid
}

sampling.model <- sampling.model(sample.parameter, log.density)

## Perform posterior sampling.
initial.partition <- rep(1,length(data))
n.draws <- 100
raw <- sample.partitions.posterior(initial.partition,sampling.model,partition.distribution,
                                   massRWSD=3,temperatureRWSD=1,n.draws)
samples.format1 <- process.samples(raw,as.matrix=TRUE)
samples.format2 <- process.samples(raw,as.matrix=TRUE, expand=TRUE)
samples.format3 <- process.samples(raw,as.matrix=FALSE)

tail(samples.format1$hyperparameters)

## Shrinkage to group means?
plot(data,apply(samples.format2$partitions,2,mean))
abline(a=0,b=1)

## Post processing to find the partition estimate.
pp <- pairwise.probabilities(raw)
est <- estimate.partition(raw, pp)
plot(confidence(pp,est))

samples.format1$hyperparameters
}
}
