## Install shallot 0.3.2-1.  Since this version is not on CRAN yet, use the following...
#> library(devtools)
#> install_github("dbdahl/shallot")

## Installing shallot will also install the rscala package.  The rscala package uses
## Scala.  If Scala is not installed on your operating system, you can install it with...
#> rscala::scalaInstall()

## But you must install Java yourself (if it is not already installed).

## Set the 'rscala.heap.maximum' if memory runs out, e.g.,
#> options(rscala.heap.maximum="4g")

library(shallot)

mass <- mass(1.0, fixed=TRUE)
discount <- discount(0.05, fixed=TRUE)

## The number of items may be large.
distance <- dist(rnorm(1000))

## But make the number of items small for the sake of the enumeration below.
distance <- dist(scale(USArrests[1:9,]))

if ( min(distance[upper.tri(distance)],na.rm=TRUE) == 0 ) stop("Oops, distances must be strictly positive.")

n.items <- attr(distance,"Size")
permutation <- permutation(n.items=n.items, fixed=FALSE)
permutation <- permutation(1:9, fixed=TRUE)
temperature <- temperature(2,fixed=FALSE)
temperature <- temperature(2,fixed=TRUE)
attraction <- attraction(permutation,decay.exponential(temperature,distance))

partition.distribution <- ewens.pitman.attraction(mass, discount, attraction)
partition.distribution <- ewens.attraction(mass, attraction)
partition.distribution <- ewens.pitman(mass, discount, n.items=9)
partition.distribution <- ewens(mass, n.items=9)

## This returns the probability mass function of the specified distribution.
pmf <- partition.pmf(partition.distribution)

## For example, the log of the probability that all items are clustered together is...
pmf(rep(1,length=partition.distribution$n.items))


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
raw <- sample.partitions.posterior(initial.partition,sampling.model,partition.distribution,massRWSD=3,temperatureRWSD=1,n.draws)
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


