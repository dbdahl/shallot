library("shallot")

mass(3.2)
mass(3.2,fixed=FALSE)
mass(2.9,0.3,fixed=FALSE)

discount(0.2)
discount(2.9,0.3,fixed=FALSE)

d <- distance(as.matrix(dist(iris[,-5])))

permutation(4,3,2,1)
permutation(sample(150))
permutation(n.items=10,fixed=FALSE)
permutation(1,2,3,fixed=FALSE)
permutation(1,2,3,fixed=TRUE)

decay.exponential(3.2)
decay.reciprocal(1,4,fixed=FALSE)
decay.subtraction(fixed=FALSE)
decay <- decay.exponential(fixed=TRUE)

attraction(d,permutation(n.items=d$n.items,fixed=FALSE),decay)
attraction(d,permutation(1:d$n.items,fixed=TRUE),decay)
attraction(d,permutation(n.items=d$n.items,fixed=FALSE),decay)
a <- attraction(d,permutation(sample(d$n.items)),decay)
a
apply(as.matrix(a),1,sum)

ewens(mass(),30)

ewens.pitman(mass(),discount(),30)

ewens.pitman.attraction(mass(),discount(),a)

ea <- ewens.attraction(mass(),a)
ea

mcmc.parameters()

samples <- collect(ea,n.draws=100)
samples <- collect(samples,n.draws=100)
y <- process(samples)
str(y)
plot(y$entropy,type="l")

collect(ewens.pitman(mass(1.0),discount(0.05),10),n.draws=100)
collect(ewens(mass(1.0),10),n.draws=100)

dist <- ewens(mass(1.5,fixed=FALSE),10)
asdf.mcmc <-    collect(dist,n.draws=100,mcmc.parameters=mcmc.parameters())
asdf.forward <- collect(dist,n.draws=100)

y1 <- process(asdf.mcmc)
y2 <- process(asdf.forward)
mean(apply(y1$partitions,1,max))
mean(apply(y2$partitions,1,max))

dist <- ewens.pitman(mass(1.5,fixed=FALSE),discount(0.05,fixed=FALSE),10)
asdf.mcmc <-    collect(dist,n.draws=1000,mcmc.parameters=mcmc.parameters())
asdf.forward <- collect(dist,n.draws=1000)
y1 <- process(asdf.mcmc)
y2 <- process(asdf.forward)
mean(apply(y1$partitions,1,max))
mean(apply(y2$partitions,1,max))

e2 <- ewens(mass(fixed=FALSE),150)
asdf2.mcmc <-    collect(e2,n.draws=300,mcmc.parameters=mcmc.parameters())
asdf2.mcmc <-    collect(asdf2.mcmc,n.draws=300)
asdf2.mcmc <-    collect(asdf2.mcmc,n.draws=300)

ea2 <- ewens.attraction(mass(fixed=TRUE),attraction(d,permutation(1:150,fixed=TRUE),decay))
asdf.mcmc <-    collect(ea2,n.draws=300,mcmc.parameters=mcmc.parameters())
asdf.mcmc <-    collect(asdf.mcmc,n.draws=300)
asdf.mcmc$proc.time
plot(process(asdf.mcmc)$entropy,type="l")
acf(process(asdf.mcmc)$entropy)

ea2 <- ewens.attraction(mass(fixed=FALSE),attraction(d,permutation(1:150,fixed=FALSE),decay))
asdf.mcmc <-    collect(ea2,n.draws=300,mcmc.parameters=mcmc.parameters())
process(asdf.mcmc)$rates







data <- c(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
sigma <- 0.1
mu0 <- 0.0 
sigma0 <- 1.0

s2 <- sigma*sigma
s02 <- sigma0*sigma0
s02Inv <- 1/s02

shallot.sample <- function(subset,parameter) {
  if ( length(subset) ) rnorm(1,mean=mu0,sd=sigma0)
  else {
    s <- sum(data[subset])
    variance <- 1 / (s02Inv + length(subset) / s2) 
    mean <- variance * (mu0 / s02 + s / s2) 
    rnorm(1,mean=mean,sd=sqrt(variance))
  }
}

shallot.log.density <- function(i,parameter) {
  dnorm(data[i],mean=parameter,sd=sigma,log=TRUE)
}

mcmc.p <- mcmc.parameters(log.density=shallot.log.density,sample=shallot.sample)

e2 <- ewens(mass(fixed=FALSE),length(data))
asdf1.mcmc <- collect(e2,n.draws=300)

asdf2.mcmc <- collect(e2,n.draws=30,mcmc.parameters=mcmc.p)
for ( i in 1:100 ) {
  mu0 <- rnorm(1,mean=0,sd=0.1)
  asdf2.mcmc <- collect(asdf2.mcmc,n.draws=1)
}
asdf2.mcmc$proc.time

asdf2.p <- process(asdf2.mcmc)

table(asdf2.p$n.subsets)/length(asdf2.p$n.subsets)

pp <- pairwise.probabilities(samples)
pp
as.matrix(pp)[1:5,1:5]

e <- estimate.partition(samples,pp,max.subsets=3)
e



