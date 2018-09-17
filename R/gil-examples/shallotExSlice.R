library(shallot)
shallot.load(java.heap.max="2G")

padata <- read.csv("pa_daily_to_Gil.csv")
source("SliceForShallot.R")

mvpa <- padata[padata$id==1,3]
mvpa <- list(mvpa)
for (i in 2:115){
  new1 <- padata[padata$id==i,3]
  new1 <- list(new1)
  mvpa <- c(mvpa,new1)
}

male <- padata[padata$id==1,8]
male <- list(male)
for (i in 2:115){
  new1 <- padata[padata$id==i,8]
  new1 <- list(new1)
  male <- c(male,new1)
}


day <- padata[padata$id==1,2]
day <- list(day)
for (i in 2:115){
  new1 <- padata[padata$id==i,2]
  new1 <- list(new1)
  day <- c(day,new1)
}

bdrivec <- padata[padata$id==1,11]
bdrivec <- list(bdrivec)
for (i in 2:115){
  new1 <- padata[padata$id==i,11]
  new1 <- list(new1)
  bdrivec <- c(bdrivec,new1)
}


satis <- padata[padata$id==1,5]
satis <- list(satis)
for (i in 2:115){
  new1 <- padata[padata$id==i,5]
  new1 <- list(new1)
  satis <- c(satis,new1)
}

shape1 <- 2
rate1 <- .4

# parameter list
# parameter[1] -- intercept model for pi
# parameter[2] -- coeficient for male, model for pi
# parameter[3] -- coeficient for day, model for pi
# parameter[4] -- coeficient for bdrivec, model for pi
# parameter[5] -- coeficient for satis, model for pi
# parameter[6] -- intercept, model for mu
# parameter[7] -- coeficient for gender, model for mu
# parameter[8] -- coeficient for bdrivec, model for mu
# parameter[9] -- coeficient for satis, model for mu
# parameter[10] -- shape parameter



uni.slice.evals <- rep(0,10)
jainneal.log.likelihood <- function(i,parameter){
  int1 <- exp(parameter[1]+parameter[2]*male[[i]]+parameter[3]*day[[i]]+parameter[4]*bdrivec[[i]]+parameter[5]*satis[[i]])
  pii <- int1/(1+int1)
  muii <- exp(parameter[6]+parameter[7]*male[[i]]+parameter[8]*bdrivec[[i]]+parameter[9]*satis[[i]])
  ratei <- parameter[10]/muii
  if (all(mvpa[[i]]==0)){
    return(sum(log(pii)))
  }
  if (all(mvpa[[i]]!=0)){
    f1 <- sum(dgamma(mvpa[[i]],shape=parameter[10],rate=ratei,log=TRUE))
    f2 <- sum(log(1-pii))
    return(f1+f2)
  }
  if (any(mvpa[[i]]==0) & any(mvpa[[i]]!=0)) {
    f1 <- sum(dgamma(mvpa[[i]][mvpa[[i]]!=0],shape=parameter[10],rate=ratei,log=TRUE))
    f2 <- sum((mvpa[[i]]!=0)*log(1-pii))
    f3 <- sum((mvpa[[i]]==0)*log(pii))
    return(f1+f2+f3)
  }
}
#accept <- rep(0,107)
#attempt <- rep(0,107)
jainneal.sample <- function(subset,parameter){
  if (length(subset)==0) {
    c(rnorm(1,0,3),rnorm(8,0,.8),rgamma(1,shape=shape1,rate=rate1))
  } 
  else {
    for (j in 1:(length(parameter)-1)){
      out <- uni.slice(parameter,subset,w=1, m=500, lower=-Inf, upper=+Inf, j)
      parameter[j] <- out
  }
    out <- uni.slice(parameter,subset,w=1, m=500, lower=0, upper=+Inf, 10)
    parameter[10] <- out
    parameter
}
}

mcmc.params <- mcmc.parameters(
  log.density.NAME="jainneal.log.likelihood",
  sample.NAME="jainneal.sample")
dist <- ewens(mass(2.0,fixed=TRUE),length(data),names=1:length(data))
dist <- ewens(mass(4,.5,fixed=FALSE),length(data),names=1:length(data))
ptm <- proc.time()
mcmc <- collect(dist,n.draws=10500,mcmc.parameters=mcmc.params)
etm <- proc.time()-ptm

mcmc <- collect(mcmc,n.draws=500)

output <- process(mcmc)

K <- apply(output$partitions,1,function(x) max(unique(x)))


extract <- function(output,index) {
  s <- lapply(output$parameters,function(x) sapply(x,function(y) y[index]))
  samples <- output$partitions
  for ( i in 1:nrow(samples) ) {
    samples[i,] <- s[[i]][output$partitions[i,]]
  }
  samples
}

b0p <- extract(output,1)
bdayp <- extract(output,2)
bmalep <- extract(output,3)
bbdrivecp <- extract(output,4)
bsatisp <- extract(output,5)
b0m <- extract(output,6)
bmalem <- extract(output,7)
bbdrivecm <- extract(output,8)
bsatism <- extract(output,9)
shapes <- extract(output,10)


props <- props[-(1:500),]
shapes <- shapes[-(1:500),]
rates <- rates[-(1:500),]

sim.props <- as.mcmc(props)
sim.shapes <- as.mcmc(shapes)
sim.rates <- as.mcmc(rates)

prop <- as.vector(props)
shape <- as.vector(shapes)
rate <- as.vector(rates)
scale <- as.vector(1/rates)
eval <- shape*scale
f1 <- kde2d(eval,shape,n=50,lims=c(0,110,0,20))
persp(f1, xlab='Exp Value', ylab = 'shape',phi = 30, theta = 30, r=2, d = 5,ticktype='detailed')
persp(f1, xlab='Exp Value', ylab = 'shape',phi = 20, theta = 110, r=2, d = 5,ticktype='detailed')
f2 <- kde2d(eval,scale,n=50,lims=c(0,110,0,50))
persp(f2, xlab='Exp Value', ylab = 'scale',phi = 20, theta = 120, r=2, d = 5,ticktype='detailed')
f3 <- kde2d(shape,scale,n=50,lims=c(0,20,0,50))
persp(f3, xlab='shape', ylab = 'scale',phi = 20, theta = 120, r=2, d = 5,ticktype='detailed')


cntr <- rep(0,107)
for (i in 2:1000){
  for (j in 1:107){
    if (shapes[i-1,j] != shapes[i,j]) cntr[j] <- cntr[j] + 1
                  }
                }
cntr
