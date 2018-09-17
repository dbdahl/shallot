#### DBD: How about the probability of a partition?  That should be an R function.

























##########################
#### Old stuff
##########################

# PPMx similarities
ppmx.similarity.categorical <- function(covariates, levels=NULL, alpha=0.5) {
  if ( ! is.matrix(covariates) ) stop("Covariates must be a matrix.")
  storage.mode(covariates) <- "character"
  n.items <- nrow(covariates)
  names <- rownames(covariates)
  covariates <- t(covariates)
  levels <- if ( is.null(levels) ) sort(unique(covariates))
  else {
    l <- as.character(levels)
    if ( ! all(unique(covariates) %in% l) ) stop("covariates contains some values that are not in levels.")
    l
  }
  alpha <- if ( length(alpha) == 1 ) as.numeric(rep(alpha,times=length(levels)))
  else {
    a <- as.numeric(alpha)
    if ( length(a) != length(levels) ) stop("alpha must have length 1 or have length equal to that of levels.") 
  }
  jobjRef <- s$.org.ddahl.shallot.distribution.PPMxSimilarityCategorical$new(covariates,levels,alpha)
  result <- list(jobjRef=jobjRef,n.items=n.items,names=names)
  class(result) <- "shallot.ppmx.similarity"
  result
}

ppmx.similarity.continuous <- function(covariates, mu0=NULL, n0=NULL, alpha=NULL, beta=NULL) {
  if ( ! is.matrix(covariates) ) stop("Covariates must be a matrix.")
  storage.mode(covariates) <- "numeric"
  n.items <- nrow(covariates)
  names <- rownames(covariates)
  covariates <- scale(covariates)
  covariates <- t(covariates)
  hyperparametersAreNull <- c(is.null(mu0),is.null(n0),is.null(alpha),is.null(beta))
  if ( any(hyperparametersAreNull) && ! all(hyperparametersAreNull) ) stop("If any hyperparameter is null, all must be null.")
  jobjRef <- if ( all(hyperparametersAreNull) ) {
    s$.org.ddahl.shallot.distribution.PPMxSimilarityContinuousMean$new(covariates)
  } else {
    mu0 <- as.numeric(mu0)[1]
    n0 <- as.numeric(n0)[1]
    alpha <- as.numeric(alpha)[1]
    beta <- as.numeric(beta)[1]
    s$.org.ddahl.shallot.distribution.PPMxSimilarityContinuousMeanVariance$new(covariates,mu0,n0,alpha,beta)
  }
  result <- list(jobjRef=jobjRef,n.items=n.items,names=names)
  class(result) <- "shallot.ppmx.similarity"
  result
}

ppmx.similarity.composite <- function(...) {
  similarities <- list(...)
  if ( length(similarities) == 0 ) stop("At least one argument must be provided.")
  x <- s$.Seq$"empty[org.ddahl.shallot.distribution.PPMxSimilarity]"()
  n.items <- similarities[[1]]$n.items
  names <- similarities[[1]]$names
  for ( i in 1:length(similarities) ) {
    if ( class(similarities[[i]]) != "shallot.ppmx.similarity" ) stop("All arguments must be PPMx similarity functions.")
    if ( similarities[[i]]$n.items != n.items ) stop("Inconsistent number of items.")
    if ( !all(similarities[[i]]$names == names)) stop("Inconsistent names.")
    x <- x$":+"(similarities[[i]]$jobjRef)
  }
  jobjRef <- s$.org.ddahl.shallot.distribution.PPMxSimilarityComposite$new(x)
  result <- list(jobjRef=jobjRef,n.items=n.items,names=names)
  class(result) <- "shallot.ppmx.similarity"
  result
}

# Product partition model with covariates (PPMx)
ppmx <- function(mass, similarity) {
  result <- list(mass=mass,attraction=similarity,n.items=nrow(covariates),names=rownames(covariates))
  class(result) <- "shallot.distribution.ppmx"
  result
}

print.shallot.distribution.ppmx <- function(x, ...) {
  cat("PPMx distribution with\n")
  cat("  ")
  print(x$mass)
  cat(paste("  similarity: ",x$attraction$jobjRef[['type']],"\n",sep=""))
}

# ppmx(mass(2.3),USArrests,"continuous")





# MCMC tuning
mcmc.parameters <- function(log.density=NULL,sample=NULL,mass.rw.sd=0.5,discount.rw.sd=0.1,permutation.grab.size=10,temperature.rw.sd=0.5,n.iterations.per.sample=1) {
  usesPredictiveDensity <- FALSE
  sampling.model <- if ( is.null(log.density) ) NULL
  else {
    callBackEngine <- s$.R
    if ( is.null(sample) ) {
      usesPredictiveDensity <- TRUE
      sample <- NULL
    }
    argsNames <- names(formals(log.density))
    if (usesPredictiveDensity) {
      if ( ( length(argsNames) != 2 ) || ( ! all( argsNames == c("i","subset") ) ) )
        stop("Signature of the log density function should be exactly c('i','subset').")
    } else {
      if ( ( length(argsNames) != 2 ) || ( ! all( argsNames == c("i","parameter") ) ) )
        stop("Signature of the log density function should be exactly c('i','parameter').")
      argsNames <- names(formals(sample))
      if ( ( length(argsNames) != 2 ) || ( ! all( argsNames == c("subset","parameter") ) ) )
        stop("Signature of the sample function should be exactly c('subset','parameter').")
    }
    log.density.REFERENCE <- scalaWrap(s,log.density)
    sample.REFERENCE <- scalaWrap(s,sample)
    s$.org.ddahl.shallot.r.RAdapter$apply(callBackEngine,as.logical(usesPredictiveDensity),log.density.REFERENCE,sample.REFERENCE)
  }
  result <- list(
    mass.rw.sd=as.double(mass.rw.sd),
    discount.rw.sd=as.double(discount.rw.sd),
    permutation.grab.size=as.integer(permutation.grab.size),
    temperature.rw.sd=as.double(temperature.rw.sd),
    n.iterations.per.sample=as.integer(n.iterations.per.sample),
    uses.predictive.density=as.logical(usesPredictiveDensity),
    sampling.model=sampling.model)
  class(result) <- "shallot.mcmc.parameters"
  result
}

# Collect
collect <- function(x,n.draws=1000,mcmc.parameters=NULL,parallel=TRUE,seed=NULL,reset=FALSE) {
  start.time <- proc.time()
  if ( n.draws < 1 ) stop("n.draws must be at least 1.")
  if ( substr(class(x),1,nchar("shallot.distribution.")) == "shallot.distribution." ) {
    pm <- x
    rdg <- .rdg()
    if ( ! is.null(seed) ) rdg$reSeed(seed)
    distribution.name <- class(pm)
    mass <- as.double(pm$mass$value)
    massShape <- as.double(pm$mass$shape)
    massRate <- as.double(pm$mass$rate)
    massFixed <- as.logical(pm$mass$fixed)
    if ( distribution.name %in% c("shallot.distribution.ewensPitman","shallot.distribution.ewensPitmanAttraction") ) {
      discount <- as.double(pm$discount$value)
      discountShape1 <- as.double(pm$discount$shape1)
      discountShape2 <- as.double(pm$discount$shape2)
      discountFixed <- as.logical(pm$discount$fixed)
    } else {
      discount <- 0.0
      discountShape1 <- 0.0
      discountShape2 <- 0.0
      discountFixed <- TRUE
    }
    if ( distribution.name %in% c("shallot.distribution.ewensAttraction","shallot.distribution.ewensPitmanAttraction","shallot.distribution.ddcrp") ) {
      distance <- pm$attraction$distance$jobjRef
      permutation <- pm$attraction$permutation$jobjRef
      permutationFixed <- as.logical(pm$attraction$permutation$fixed)
      decayString <- as.character(pm$attraction$decay$type)
      decayMultiplier <- as.double(pm$attraction$decay$multiplier)
      temperature <- as.double(pm$attraction$decay$temperature$value)
      temperatureShape <- as.double(pm$attraction$decay$temperature$shape)
      temperatureRate <- as.double(pm$attraction$decay$temperature$rate)
      temperatureFixed <- as.logical(pm$attraction$decay$temperature$fixed)
    } else {
      distance <- s$.null
      permutation <- s$.null
      permutationFixed <- TRUE
      decayString <- ""
      decayMultiplier <- 1.0
      temperature <- 0.0
      temperatureShape <- 0.0
      temperatureRate <- 0.0
      temperatureFixed <- TRUE
    }
    if ( distribution.name == "shallot.distribution.ppmx" ) {
      ppmxSimilarity <- pm$attraction$jobjRef
    } else {
      ppmxSimilarity <- s$.null
    }
    if ( is.null(mcmc.parameters) ) {
      sampler <- s$.org.ddahl.shallot.r.RInterface$makeForwardSamplerForR(
        distribution.name,as.integer(pm$n.items),distance,
        mass,massShape,massRate,massFixed,
        discount,discountShape1,discountShape2,discountFixed,
        permutation,permutationFixed,
        decayString,decayMultiplier,
        temperature,temperatureShape,temperatureRate,temperatureFixed,ppmxSimilarity)
      augmentedSamples <- NULL
      samples <- s$.org.ddahl.shallot.r.RInterface$sampleForward(as.integer(n.draws),rdg,sampler,as.logical(parallel))
    } else {
      sampling.model <- if ( is.null(mcmc.parameters$sampling.model) ) {
        s$.null
      }
      else mcmc.parameters$sampling.model
      sampler <- s$.org.ddahl.shallot.r.RInterface$makeMCMCSamplerForR(
        distribution.name,as.integer(pm$n.items),distance,
        s$.Tuple5$apply(mass,massShape,massRate,massFixed,mcmc.parameters$mass.rw.sd),
        discount,discountShape1,discountShape2,discountFixed,mcmc.parameters$discount.rw.sd,
        permutation,permutationFixed,as.integer(min(mcmc.parameters$permutation.grab.size,pm$n.items)),
        decayString,decayMultiplier,
        temperature,temperatureShape,temperatureRate,temperatureFixed,mcmc.parameters$temperature.rw.sd,
        sampling.model)
      augmentedSamples <- s$.org.ddahl.shallot.r.RInterface$sampleMCMC(as.integer(n.draws),as.integer(mcmc.parameters$n.iterations.per.sample),rdg,sampler)
      samples <- s$.org.ddahl.shallot.r.RInterface$extractSamples(augmentedSamples)
    }
    n.samples <- samples$length()
    result <- list(proc.time=proc.time()-start.time,n.items=pm$n.items,n.samples=n.samples,names=x$names,mcmc.parameters=mcmc.parameters,jobjRef=list(rdg=rdg,samples=samples,sampler=sampler,augmentedSamples=augmentedSamples))
  } else if ( class(x) == "shallot.samples.raw" ) {
    if ( ! is.null(mcmc.parameters) ) stop("mcmc.parameters must be null when continuing.")
    if ( ! is.null(seed) ) stop("seed must be null when continuing.")
    if ( is.null(x$mcmc.parameters) ) {
      augmentedSamples <- NULL
      if ( reset ) {
        samples <- x$jobjRef$samples$take(as.integer(1))
      } else {
        samples <- s$.org.ddahl.shallot.r.RInterface$sampleForward(as.integer(n.draws),x$jobjRef$rdg,x$jobjRef$sampler,as.logical(parallel))
        samples <- s$.org.ddahl.shallot.r.RInterface$mergeSamples(samples,x$jobjRef$samples)
      }
    } else {
      if ( reset ) {
        augmentedSamples <- x$jobjRef$augmentedSamples$take(as.integer(1))
        samples <- x$jobjRef$samples$take(as.integer(1))
      } else {
        augmentedSamples <- s$.org.ddahl.shallot.r.RInterface$sampleMCMC(as.integer(n.draws),x$mcmc.parameters$n.iterations.per.sample,x$jobjRef$rdg,x$jobjRef$sampler)
        samples <- s$.org.ddahl.shallot.r.RInterface$extractSamples(augmentedSamples)
        samples <- s$.org.ddahl.shallot.r.RInterface$mergeSamples(samples,x$jobjRef$samples)
        augmentedSamples <- s$.org.ddahl.shallot.r.RInterface$mergeAugmentedSamples(augmentedSamples,x$jobjRef$augmentedSamples)
      }
    }
    if ( reset ) x$jobjRef$sampler$reset()
    n.samples <- samples$length()
    result <- list(proc.time=x$proc.time+(proc.time()-start.time),n.items=x$n.items,n.samples=n.samples,names=x$names,mcmc.parameters=x$mcmc.parameters,jobjRef=list(rdg=x$jobjRef$rdg,samples=samples,sampler=x$jobjRef$sampler,augmentedSamples=augmentedSamples))
  } else stop("First argument is not recognized.")
  class(result) <- "shallot.samples.raw"
  result
}

# Process
process <- function(x, ...) {
  start.time <- proc.time()
  withParameters <- ! is.null(x$jobjRef$augmentedSamples) && ! is.null(x$mcmc.parameters) && ! x$mcmc.parameters$uses.predictive.density
  samples <- x$jobjRef$samples$reverse()
  labelsWithParameters <- samples$map(s[['env']]$mapper.toLabelsWithParameters)
  partitions <- labelsWithParameters$map(s[['env']]$mapper.toLabels)$toArray()
  colnames(partitions) <- x$names
  n.subsets <- samples$map(s[['env']]$mapper.nSubsets)$toArray()
  entropy <- samples$map(s[['env']]$mapper.entropy)$toArray()
  result <- list(partitions=partitions,n.subsets=n.subsets,entropy=entropy)
  if ( ! is.null(x$jobjRef$augmentedSamples) ) {
    recordMass <- x$jobjRef$sampler$recordMass()
    recordDiscount <- x$jobjRef$sampler$recordDiscount()
    recordPermutation <- x$jobjRef$sampler$recordPermutation()
    recordTemperature <- x$jobjRef$sampler$recordTemperature()
    augmentedSamples <- x$jobjRef$augmentedSamples$reverse()
    if ( recordMass )        result <- c(result,list(mass=       augmentedSamples$map(s[['env']]$mapper.mass       )$toArray()))
    if ( recordDiscount )    result <- c(result,list(discount=   augmentedSamples$map(s[['env']]$mapper.discount   )$toArray()))
    if ( recordTemperature ) result <- c(result,list(temperature=augmentedSamples$map(s[['env']]$mapper.temperature)$toArray()))
    if ( ! is.null(x$mcmc.parameters) ) {
      rates <- x$jobjRef$sampler$rates()
      names(rates) <- c("mass","discount","permutation","temperature")
      rates <- rates[c(recordMass,recordDiscount,recordPermutation,recordTemperature)]
      if ( withParameters ) {
        parameterReferences <- labelsWithParameters$map(s[['env']]$mapper.toParameters)$toArray()
        parameters <- list()
        length(parameters) <- length(n.subsets)
        for ( i in 1:length(parameters) ) {
          ri <- parameterReferences$apply(as.integer(i-1))
          tl <- list()
          nParams <- n.subsets[i]
          length(tl) <- nParams
          for ( j in 1:nParams ) {
            tl[[j]] <- scalaUnwrap(s,ri$apply(as.integer(j-1)))
          }
          parameters[[i]] <- tl
        }
      } else parameters <- NULL
      result <- c(result,list(parameters=parameters,rates=rates))
    }
  }
  result <- c(result,list(proc.time=proc.time()-start.time))
  class(result) <- "shallot.samples"
  result
}

print.shallot.samples <- function(x, ...) {
  cat("samples --- list with elements: ",paste(names(x),collapse=", "),".\n",sep="")
}

# Current state
current.state <- function(x) {
  withParameters <- ! is.null(x$jobjRef$augmentedSamples) && ! is.null(x$mcmc.parameters) && ! x$mcmc.parameters$uses.predictive.density
  sample <- x$jobjRef$samples$head()
  labelsWithParameters <- s[['env']]$mapper.toLabelsWithParameters$apply(sample)
  partition <- labelsWithParameters$"_1"()
  names(partition) <- x$names
  if ( ! is.null(x$jobjRef$augmentedSamples) && ! is.null(x$mcmc.parameters) ) {
    if ( withParameters ) {
      parameterReferences <- labelsWithParameters$"_2"()
      parameters <- list()
      nParams <- parameterReferences$length()
      length(parameters) <- nParams
      for ( j in 1:nParams ) {
        parameters[[j]] <- scalaUnwrap(s,parameterReferences$apply(as.integer(j-1)))
      }
    } else parameters <- NULL
    result <- list(partition=partition,parameters=parameters)
  } else result <- list(partition=partition)
  class(result) <- "shallot.state"
  result
}

print.shallot.state <- function(x, ...) {
  cat("partition\n")
  print(x$partition)
  if ( ! is.null(x$parameters) ) {
    cat("\n")
    cat("parameters\n")
    print(x$parameters)
  }
}

