library("shallot")

# Example from Neal (2000)
data <- c(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
sigma <- 0.1
mu0 <- 0.0 
sigma0 <- 1.0

s2 <- sigma*sigma
s02 <- sigma0*sigma0
s02Inv <- 1/s02

jainneal.sample <- function(subset=c()) {
  if ( length(subset) ) rnorm(1,mean=mu0,sd=sigma0)
  else {
    s <- sum(data[subset])
    variance <- 1 / (s02Inv + length(subset) / s2) 
    mean <- variance * (mu0 / s02 + s / s2) 
    rnorm(1,mean=mean,sd=sqrt(variance))
  }
}

jainneal.log.density <- function(i,parameter) {
  dnorm(data[i],mean=parameter,sd=sigma,log=TRUE)
}

mcmc.p <- mcmc.parameters(log.density.NAME="jainneal.log.density",sample.NAME="jainneal.sample")

mass(3.2)
mass(3.2,fixed=FALSE)
mass(2.9,0.3,fixed=FALSE)
dist <- ewens(mass(fixed=FALSE),length(data))

mcmc1 <- collect(dist,n.draws=1000)
current.state(mcmc1)

mcmc2 <- collect(dist,n.draws=30,mcmc.parameters=mcmc.p)
current.state(mcmc2)

mcmc2$n.samples
mcmc2 <- collect(mcmc2,reset=TRUE)
mcmc2$n.samples

for ( i in 1:1000 ) {
  mu0 <- rnorm(1,mean=0,sd=0.1)  # Apply whatever update is needed.
  mcmc2 <- collect(mcmc2,n.draws=1)
}
mcmc2$proc.time

x <- process(mcmc2)

hist(x$n.subsets)
plot(x$entropy,type="l")

pp <- pairwise.probabilities(mcmc2)
pp
as.matrix(pp)[1:5,1:5]

e <- estimate.partition(mcmc2,pp,max.subsets=3)
e
e$partition



