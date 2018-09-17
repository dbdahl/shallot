library(rscala)
scala("commonsMath")

mass <- mass(1.0, fixed=TRUE)

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

## Sampling model of Neal (JCGS, 2000)
## Function to perform an MCMC update of the parameter.
sample.parameter <- function(indices=scalaType("I1"), parameter=scalaType("D1")) {
  sum <- sum(data[indices])
  variance <- 1 / (s02Inv + length(indices) / s2)
  mean <- variance * (mu0 / s02 + sum / s2)
  rnorm(1, mean=mean, sd=sqrt(variance))
}
s(data=data,s02Inv=s02Inv,s2=s2,mu0=mu0,s02=s02) ^ sample.parameter
s(data,s02Inv,s2,mu0,s02) ^ sample.parameter

## Function to evaluate the likelihood contribution for an observation.
log.density <- function(i, indices, parameter) {
  resid <- data[i] - parameter
  c * resid * resid
}

sampling.model <- sampling.model(sample.parameter, log.density)

partition.distribution <- ewens(mass, n.items=length(data))

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

samples.format1$hyperparameters


