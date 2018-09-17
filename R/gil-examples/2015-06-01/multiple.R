library("shallot")
#shallot.load(java.opts="-Xdebug -Xrunjdwp:transport=dt_socket,address=8787,server=y,suspend=n",debug.filename="/tmp/rscala")
shallot.load()
rscala::intpSettings(shallot:::E$s,quiet=FALSE)


source("SliceForShallot.R")
padata <- read.csv('pa_daily_to_Gil.csv')
padata1 <- padata[!is.na(padata$satis),]
mvpa <- padata1[padata1$id==1,3]
mvpa <- list(mvpa)
for (i in 2:115){
  new1 <- padata1[padata1$id==i,3]
  new1 <- list(new1)
  mvpa <- c(mvpa,new1)
}

data <- mvpa

male <- padata1[padata1$id==1,8]
male <- list(male)
for (i in 2:115){
  new1 <- padata1[padata1$id==i,8]
  new1 <- list(new1)
  male <- c(male,new1)
}


day <- padata1[padata1$id==1,2]
day <- list(day)
for (i in 2:115){
  new1 <- padata1[padata1$id==i,2]
  new1 <- list(new1)
  day <- c(day,new1)
}

bdrivec <- padata1[padata1$id==1,11]
bdrivec <- list(bdrivec)
for (i in 2:115){
  new1 <- padata1[padata1$id==i,11]
  new1 <- list(new1)
  bdrivec <- c(bdrivec,new1)
}


satis <- padata1[padata1$id==1,5]
satis <- list(satis)
for (i in 2:115){
  new1 <- padata1[padata1$id==i,5]
  new1 <- list(new1)
  satis <- c(satis,new1)
}

shape1 <- 2
rate1 <- .4
uni.slice.evals <- rep(0,10)

loglike.i  <- function(i, parameter) {   # where parameters is a list/vector for all the parameters
  cat("i = ",i,": ",paste(parameter,collapse=","),"\n")
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
    f1 <- sum(dgamma(mvpa[[i]][mvpa[[i]]!=0],shape=parameter[10],rate=ratei[mvpa[[i]]!=0],log=TRUE))
    f2 <- sum((mvpa[[i]]!=0)*log(1-pii))
    f3 <- sum((mvpa[[i]]==0)*log(pii))
    return(f1+f2+f3)
  }
}  

loglike    <- function()    {
  sum <- 0.0
  for ( i in 1:nrow(data) ) {
    p <- patch(i)
    sum <- sum + loglike.i(i, p)
  }
  sum
}

params <- c(1,0,0,0,0,1,0,0,0,10)   # p1, p2, ..., p10
p1 <- p6 <- p10 <- rep(1,length(data))

patch <- function(i) {
  p <- params
  p[1]  <- p1[i]
  p[6]  <- p6[i]
  p[10] <- p10[i]
  p
}

expand <- function(mcmc) {
  p <- current.state(mcmc)
  unlist(p$parameters)[p$partition]    
}

mkLogLike <- function(j) {
  function(i, parameter) {
    p <- patch(i)
    p[j] <- parameter
    loglike.i(i, p)
  }
}

sample.1 <- function(subset, parameter) {
  if ( length(subset) == 0 ) {
    rnorm(1,0,1)
  } else {
    p1u <- p1
    p1u[subset] <- parameter    
    paraM <- cbind(p1u,params[2],params[3],params[4],params[5],p6,params[7],params[8],params[9],p10)
    uni.slice(paraM,subset,w=1, m=500, lower=-5, upper=+5, 1)
  }
}

sample.6 <- function(subset, parameter) {
  if ( length(subset) == 0 ) {
    rnorm(1,0,.1)
  } else {
    p6u <- p6
    p6u[subset] <- parameter
    paraM <- cbind(p1,params[2],params[3],params[4],params[5],p6u,params[7],params[8],params[9],p10)
    uni.slice(paraM,subset,w=1, m=500, lower=-1, upper=+1, 6)
  }
}

sample.10 <- function(subset, parameter) {
  if ( length(subset) == 0 ) {
    rgamma(1,shape=shape1,rate=rate1)
  } else {
    p10u <- p10
    p10u[subset] <- parameter    
    paraM <- cbind(p1,params[2],params[3],params[4],params[5],p6,params[7],params[8],params[9],p10u)
    uni.slice(paraM,subset,w=1, m=500, lower=0, upper=+Inf, 10)
  }
}

mcmc.parameters.1  <- mcmc.parameters(log.density=mkLogLike(1),  sample=sample.1)
mcmc.parameters.6  <- mcmc.parameters(log.density=mkLogLike(6),  sample=sample.6)
mcmc.parameters.10 <- mcmc.parameters(log.density=mkLogLike(10), sample=sample.10)

mcmc.1  <- collect(ewens(mass(1.0, fixed = TRUE), length(data), names = names(data)),
                   n.draws = 1, mcmc.parameters = mcmc.parameters.1)

mcmc.6  <- collect(ewens(mass(1.0, fixed = TRUE), length(data), names = names(data)),
                   n.draws = 1, mcmc.parameters = mcmc.parameters.6)

mcmc.10 <- collect(ewens(mass(1.0, fixed = TRUE), length(data), names = names(data)),
                   n.draws = 1, mcmc.parameters = mcmc.parameters.10)

nSamples <- 10
results <- matrix(NA, nrow = nSamples, ncol = length(params))  # Ignore columns 1, 6, 10
results[1, ] <- params


for ( r in 2:nSamples ) {
  paraM <- cbind(p1,params[2],params[3],params[4],params[5],p6,params[7],params[8],params[9],p10)

  # Update p1
  mcmc.1 <-  collect(mcmc.1,  n.draws = 1)
  paraM[,1] <- p1 <- expand(mcmc.1)
  
  # Update p2,...,p5
  params[2] <- uni.slice(paraM,1:115,w=1, m=500, lower=-.5, upper=+.5, 2)
  paraM[,2] <- rep(params[2],115)
  params[3] <- uni.slice(paraM,1:115,w=1, m=100, lower=-.1, upper=+.1, 3)
  paraM[,3] <- rep(params[3],115)
  params[4] <- uni.slice(paraM,1:115,w=1, m=500, lower=-1, upper=+1, 4)
  paraM[,4] <- rep(params[4],115)
  params[5] <- uni.slice(paraM,1:115,w=1, m=500, lower=-1, upper=+1, 5)
  paraM[,5] <- rep(params[5],115)
  
  # Update p6
  mcmc.6 <-  collect(mcmc.6,  n.draws = 1)  
  paraM[,6] <- p6 <- expand(mcmc.6)
  
  # Update p7,...,p9
  params[7] <- uni.slice(paraM,1:115,w=1, m=500, lower=-1, upper=+1, 7)
  paraM[,7] <- rep(params[7],115)
  params[8] <- uni.slice(paraM,1:115,w=1, m=500, lower=-1, upper=+1, 8)
  paraM[,8] <- rep(param[8],115)
  params[9] <- uni.slice(paraM,1:115,w=1, m=500, lower=-1, upper=+1, 9)
  
  # Update p10
  mcmc.10 <- collect(mcmc.10, n.draws = 1)
  paraM[,10] <- p10 <- expand(mcmc.10)
  
  results[r, ] <- params
}
output.1  <- process(mcmc.1)
output.6  <- process(mcmc.6)
output.10 <- process(mcmc.10)



