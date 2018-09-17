loglike.i  <- function(i, parameters) {   # where parameters is a list/vector for all the parameters
  
}

loglike    <- function(parameters)    {   # where parameters is a list/vector for all the parameters
  sum <- 0.0
  for ( i in 1:nrow(data) ) {
    sum <- sum + loglike.i(i, parameters)
  }
  sum
}

params <- rep(1,10)   # p1, p2, ..., p10

patch <- function(i, parameters) {
  p <- parameters
  p1 <- current.state(mcmc.1)
  p[1] <- p1$parameters[[p1$partition[i]]]
  p6 <- current.state(mcmc.6)
  p[6] <- p6$parameters[[p6$partition[i]]]
  p10 <- current.state(mcmc.10)
  p[10] <- p10$parameters[[p10$partition[i]]]
  p
}

mkLogLike <- function(j) {
  function(i, parameter) {
    p <- patch(i, params)
    p[j] <- parameter
    loglike.i(i, p)
  }
}

sample.1 <- function(subset, parameter) {
  if ( length(subset) == 0 ) {
    ##### Do stuff
  } else {
    ##### Do other stuff
  }
}

sample.6 <- function(subset, parameter) {
  if ( length(subset) == 0 ) {
    ##### Do stuff
  } else {
    ##### Do other stuff
  }
}

sample.10 <- function(subset, parameter) {
  if ( length(subset) == 0 ) {
    ##### Do stuff
  } else {
    ##### Do other stuff
  }
}

mcmc.1  <- collect(ewens(mass(1.0, fixed = TRUE), length(data), names = names(data)),
                   n.draws = 1,
                   mcmc.parameters = mcmc.parameters(log.density=mkLogLike(1), sample=sample.1))

mcmc.6  <- collect(ewens(mass(1.0, fixed = TRUE), length(data), names = names(data)),
                   n.draws = 1,
                   mcmc.parameters = mcmc.parameters(log.density=mkLogLike(6), sample=sample.6))

mcmc.10 <- collect(ewens(mass(1.0, fixed = TRUE), length(data), names = names(data)),
                   n.draws = 1,
                   mcmc.parameters = mcmc.parameters(log.density=mkLogLike(10), sample=sample.10))

results <- matrix(NA, nrow = nSamples, ncol = length(params))  # Ignore columns 1, 6, 10
results[1, ] <- params
for ( r in 2:nSamples ) {
  mcmc.1 <-  collect(mcmc.1,  n.draws = 1)
  # Update p2,...,p5
  mcmc.6 <-  collect(mcmc.6,  n.draws = 1)
  # Update p7,...,p9
  mcmc.10 <- collect(mcmc.10, n.draws = 1)
  results[r, ] <- params
}
output.1  <- process(mcmc.1)
output.6  <- process(mcmc.6)
output.10 <- process(mcmc.10)


