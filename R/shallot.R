### Development
# library(rscala); library(shallot); s <- shallot:::s

##########################
#### Simplified
##########################

# Random Data Generator
.rdg <- function() s %.!% 'rdg()'



# Adjusted Rand Index
adj.rand.index <- function(c1,c2) {
  n <- length(c1)
  if ( length(c2) != n ) stop("Clusterings must be the same length.")
  t1  <- table(c1)
  t2  <- table(c2)
  t12 <- table(c1,c2)
  expected <- sum(choose(t1,2)) * sum(choose(t2,2)) / choose(n,2)
  numerator <- sum(choose(t12,2)) - expected
  denominator <- 0.5*(sum(choose(t1,2))+sum(choose(t2,2))) - expected
  numerator / denominator
}

# adj.rand.index(c(1,1,1,2,2,2,3,3,3),c(1,1,1,1,1,1,2,2,2))


#### To be added: "variation of information" by Meila (2007).


# Mass
mass <- function(...,fixed=TRUE) {
  x <- c(...)
  fixed <- as.logical(fixed)
  value <- 1.2
  shape <- 2.5
  rate <- 2
  if ( length(x) == 0 ) {
  } else if ( length(x) == 1 ) {
    value <- as.double(x[1])
    if ( value <= 0 ) stop("'value' must be positive.")
  } else if ( length(x) == 2 ) {
    if ( fixed ) stop("'fixed' should be FALSE if distribution parameters are specified.")
    shape <- as.double(x[1])
    if ( shape <= 0.0 ) stop("'shape' must be positive.")
    rate <- as.double(x[2])
    if ( rate <= 0.0 ) stop("'rate' must be positive.")
    value <- shape/rate
  } else stop("Incorrect number of arguments.")
  result <- list(value=value,shape=shape,rate=rate,fixed=fixed)
  class(result) <- "shallot.mass"
  result
}

# Not exported:  Mass wrapper
.mass <- function(mass=NULL) s %.!% '
  parameter.Mass(R.evalD0(mass+"$value"))
'

.massFactory <- function(mass=NULL) s %.!% '
  if ( R.evalL0(mass+"$fixed")  ) {
    parameter.Mass.factory(R.evalD0(mass+"$value"))
  } else {
    parameter.Mass.factory(R.evalD0(mass+"$shape"),R.evalD0(mass+"$rate"),rdg())
  }
'

print.shallot.mass <- function(x, ...) {
  if ( x$fixed ) cat("mass fixed at ",x$value,"\n",sep="")
  else cat("mass at ",x$value," and distributed Gamma(shape=",x$shape,",rate=",x$rate,")\n",sep="")
}

# mass()
# mass(2.3)
# mass(1,fixed=FALSE)
# mass(3,4,fixed=FALSE)



# Discount
discount <- function(...,fixed=TRUE) {
  x <- c(...)
  fixed <- as.logical(fixed)
  value <- 0.05
  shape1 <- 1.0
  shape2 <- 1.0
  if ( length(x) == 0 ) {
  } else if ( length(x) == 1 ) {
    value <- as.double(x[1])
    if ( ( value < 0 ) || ( value >= 1.0 ) ) stop("'value' must be in [0,1).")
  } else if ( length(x) == 2 ) {
    if ( fixed ) stop("'fixed' should be FALSE if distribution parameters are specified.")
    shape1 <- as.double(x[1])
    if ( shape1 <= 0.0 ) stop("'shape1' must be positive.")
    shape2 <- as.double(x[2])
    if ( shape2 <= 0.0 ) stop("'shape2' must be positive.")
    value <- shape1/(shape1+shape2)
  } else stop("Incorrect number of arguments.")
  result <- list(value=value,shape1=shape1,shape2=shape2,fixed=fixed)
  class(result) <- "shallot.discount"
  result
}

# Not exported:  Discount wrapper 
.discount <- function(discount=NULL) s %.!% '
  parameter.Discount(R.evalD0(discount+"$value"))
'

.discountFactory <- function(discount=NULL) s %.!% '
  if ( R.evalL0(discount+"$fixed")  ) {
    parameter.Discount.factory(R.evalD0(discount+"$value"))
  } else {
    parameter.Discount.factory(R.evalD0(discount+"$shape1"),R.evalD0(discount+"$shape2"),rdg())
  }
'

print.shallot.discount <- function(x, ...) {
  if ( x$fixed ) cat("discount fixed at ",x$value,"\n",sep="")
  else cat("discount at ",x$value," and distributed Beta(shape1=",x$shape1,",shape2=",x$shape2,")\n",sep="")
}

# discount()
# discount(0.1)
# discount(0.1,fixed=FALSE)
# discount(3,4,fixed=FALSE)



# Permutation
permutation <- function(...,n.items=NULL,fixed=TRUE) {
  x <- c(...)
  fixed <- as.logical(fixed)
  if ( length(x) == 0 ) {
    if ( is.null(n.items) ) stop("'n.items' must be specified if permutation is not given.")
    n.items <- as.integer(n.items[1])
    if ( fixed ) stop("'fixed' must be FALSE if permutation is not given.")
    x <- sample(1:n.items)
  } else {
    if ( is.null(n.items) ) n.items <- length(x)
    n.items <- as.integer(n.items[1])
    if ( n.items != length(x) ) stop("'n.items' is not equal to the length of the permutation.")
    if ( length(unique(x)) != n.items ) stop("A permutation cannot have repeated values.")
    if ( min(x) < 1 ) stop("The smallest value in a permutation should be 1.")
    if ( max(x) > n.items ) stop("The largest value in a permutation should be its length.")
  }
  result <- list(value=as.integer(x),n.items=n.items,fixed=fixed)
  class(result) <- "shallot.permutation"
  result
}

.permutation <- function(permutation) {
  s$.parameter.Permutation$apply(permutation$value-1L)
}

.permutationFactory <- function(permutation) {
  if ( permutation$fixed ) {
    s$.parameter.Permutation$factory(permutation$value-1L)
  } else {
    s$.parameter.Permutation$factory(permutation$n.items,.rdg())
  }
}

print.shallot.permutation <- function(x, ...) {
  if ( x$fixed ) {
    cat("permutation of ",x$n.items," items fixed at\n",sep="")
    print(x$value)
  }
  else {
    cat("permutation of ",x$n.items," items distributed uniformly whose value is\n",sep="")
    print(x$value)
  }
}

# permutation(n.items=30,fixed=FALSE)
# permutation(2,3,4,1)
# permutation(2,3,4,1,fixed=FALSE)



# Temperature
temperature <- function(...,fixed=TRUE) {
  x <- c(...)
  fixed <- as.logical(fixed)
  value <- 3
  shape <- 2
  rate <- 0.5
  if ( length(x) == 0 ) {
  } else if ( length(x) == 1 ) {
    if ( ! fixed ) stop("'fixed' should be TRUE if value is specified.")
    value <- as.double(x[1])
  } else if ( length(x) == 2 ) {
    if ( fixed ) stop("'fixed' should be FALSE if distribution parameters are specified.")
    shape <- as.double(x[1])
    if ( shape <= 0.0 ) stop("'shape' must be positive.")
    rate <- as.double(x[2])
    if ( rate <= 0.0 ) stop("'rate' must be positive.")
    value <- shape/rate
  } else stop("Incorrect number of arguments.")
  result <- list(value=value,shape=shape,rate=rate,fixed=fixed)
  class(result) <- "shallot.temperature"
  result
}

print.shallot.temperature <- function(x, ...) {
  if ( x$fixed ) cat("temperature fixed at ",x$value,"\n",sep="")
  else cat("temperature at ",x$value," and distributed Gamma(shape=",x$shape,",rate=",x$rate,")\n",sep="")
}

# temperature()
# temperature(2,4,fixed=FALSE)
# temperature(2)



# Decay functions
decay.reciprocal <- function(temperature,distance) {
  if ( ! inherits(distance,"dist") ) stop("'distance' must be of class 'dist'")
  x1 <- -log(.Machine$double.xmin)/log(max(distance))
  x2 <- -(log(.Machine$double.xmax)-log(attr(distance,"Size")))/log(min(distance))
  x <- Inf
  if ( x1 > 0 ) x <- min(x,x1)
  if ( x2 > 0 ) x <- min(x,x2)
  max.temperature <- x
  decay.generic(temperature,distance,type="reciprocal",max.temperature=max.temperature)
}

decay.exponential <- function(temperature,distance) {
  if ( ! inherits(distance,"dist") ) stop("'distance' must be of class 'dist'")
  x1 <- -log(.Machine$double.xmin)/max(distance)
  x2 <- -(log(.Machine$double.xmax)-log(attr(distance,"Size")))/min(distance)
  x <- Inf
  if ( x1 > 0 ) x <- min(x,x1)
  if ( x2 > 0 ) x <- min(x,x2)
  max.temperature <- x
  decay.generic(temperature,distance,type="exponential",max.temperature=max.temperature)
}

decay.subtraction <- function(temperature,distance,multiplier=1.01) {
  if ( ! inherits(distance,"dist") ) stop("'distance' must be of class 'dist'")
  if ( ( length(multiplier) != 1 ) || ( ! is.numeric(multiplier) ) || ( multiplier <= 1.0 ) ) stop("'multiplier' must be a scalar greater than 1.")
  max <- max(distance)
  min <- min(distance)
  max.dist <- multiplier*max
  d1 <- ( max.dist - min )
  d2 <- ( max.dist - max )
  x1 <- log(.Machine$double.xmin)/log(d2)
  x2 <- (log(.Machine$double.xmax)-log(attr(distance,"Size")))/log(d1)
  x <- Inf
  if ( x1 > 0 ) x <- min(x,x1)
  if ( x2 > 0 ) x <- min(x,x2)
  max.temperature <- x
  decay.generic(temperature,distance,type="subtraction",max.temperature=max.temperature,max.distance=max.dist)
}

decay.generic <- function(temperature,distance,type,max.temperature,max.distance=NULL) {
  if ( ! inherits(temperature,"shallot.temperature") ) stop("'temperature' must be of class 'shallot.temperature'")
  result <- list(temperature=temperature,distance=distance,type=type,max.temperature=max.temperature,max.distance=max.distance)
  class(result) <- "shallot.decay"
  result
}

.decay <- function(decay) {
  temp <- min(decay$temperature$value,decay$max.temperature)
       if ( decay$type == "reciprocal" )  s$.parameter.decay.ReciprocalDecay$new(temp)
  else if ( decay$type == "exponential" ) s$.parameter.decay.ExponentialDecay$new(temp)
  else if ( decay$type == "subtraction" ) s$.parameter.decay.SubtractionDecay$new(temp,decay$max.distance)
}

.decayFactory <- function(decay) {
  if ( decay$temperature$fixed ) {
    temp <- min(decay$temperature$value,decay$max.temperature)
         if ( decay$type == "reciprocal" )  s$.parameter.decay.ReciprocalDecayFactory$factory(temp)
    else if ( decay$type == "exponential" ) s$.parameter.decay.ExponentialDecayFactory$factory(temp)
    else if ( decay$type == "subtraction" ) s$.parameter.decay.SubtractionDecayFactory$new(decay$max.distance)$factory(temp)
  } else {
    shape <- decay$temperature$shape
    rate <- decay$temperature$rate
         if ( decay$type == "reciprocal" )  s$.parameter.decay.ReciprocalDecayFactory$factory(shape,rate,.rdg())
    else if ( decay$type == "exponential" ) s$.parameter.decay.ExponentialDecayFactory$factory(shape,rate,.rdg())
    else if ( decay$type == "subtraction" ) s$.parameter.decay.SubtractionDecayFactory$new(decay$max.distance)$factory(shape,rate,.rdg())
  }
}

print.shallot.decay <- function(x, ...) {
  cat(x$type,"decay function with ")
  if ( x$temperature$fixed ) cat("temperature fixed at ",x$temperature$value,"\n",sep="")
  else cat("temperature sampled from Gamma(shape=",x$temperature$shape,",rate=",x$temperature$rate,")\n",sep="")
}

# decay.reciprocal()
# decay.reciprocal(4)
# decay.reciprocal(2,4,fixed=FALSE)
# decay.reciprocal(Inf,distance=d)
# decay.exponential()
# decay.exponential(4)
# decay.exponential(2,4,fixed=FALSE)
# decay.exponential(Inf,distance=d)
# decay.subtraction(distance=d)
# decay.subtraction(4,distance=d)
# decay.subtraction(2,4,fixed=FALSE,distance=d)
# decay.subtraction(Inf,distance=d)



# Attraction
attraction <- function(permutation, decay) {
  if ( ( is.vector(permutation) ) && ( length(permutation) == 1 ) && ( missing(decay) ) ) {
    result <- list(constant=TRUE,n.items=as.integer(permutation))
  } else {
    if ( ! inherits(permutation,"shallot.permutation") ) stop("'permutation' must be of class 'shallot.permutation'")
    if ( ! inherits(decay,"shallot.decay") ) stop("'decay' must be of class 'shallot.decay'")
    result <- list(constant=FALSE,n.items=permutation$n.items, permutation=permutation, decay=decay, names=attr(decay$distance,"Labels"))
  }
  class(result) <- "shallot.attraction"
  result
}

.distance <- function(distance) {
  s$.parameter.Distance$apply(as.matrix(distance),FALSE)
}

.attraction <- function(attraction) {
  if ( attraction$constant ) {
    s$.distribution.Attraction$constant(attraction$n.items)
  } else {
    tryCatch(
      s$.distribution.Attraction$apply(.distance(as.matrix(attraction$decay$distance)),.permutation(attraction$permutation),.decay(attraction$decay))
    ,error = function(e) stop("Attraction is invalid because 'distance' and 'decay' appear to be incompatible.  Perhaps lower the 'temperature'."))
  }
}

.attractionFactory <- function(attraction) {
  if ( attraction$constant ) {
    s$.distribution.Attraction$factory(attraction$n.items)
  } else if ( attraction$permutation$fixed && attraction$decay$temperature$fixed ) {
    tryCatch(
      s$.distribution.Attraction$factory(.distance(as.matrix(attraction$decay$distance)),.permutation(attraction$permutation),.decay(attraction$decay))
    ,error = function(e) stop("Attraction is invalid because 'distance' and 'decay' appear to be incompatible.  Perhaps lower the 'temperature'."))
  } else if ( attraction$permutation$fixed ) {
    tryCatch(
      s$.distribution.Attraction$factory(.distance(as.matrix(attraction$decay$distance)),.permutation(attraction$permutation),.decayFactory(attraction$decay))
    ,error = function(e) stop("Attraction is invalid because 'distance' and 'decay' appear to be incompatible.  Perhaps lower the 'temperature'."))
  } else if ( attraction$decay$temperature$fixed ) {
    tryCatch(
      s$.distribution.Attraction$factory(.distance(as.matrix(attraction$decay$distance)),.permutationFactory(attraction$permutation),.decay(attraction$decay))
    ,error = function(e) stop("Attraction is invalid because 'distance' and 'decay' appear to be incompatible.  Perhaps lower the 'temperature'."))
  } else {
    tryCatch(
      s$.distribution.Attraction$factory(.distance(as.matrix(attraction$decay$distance)),.permutationFactory(attraction$permutation),.decayFactory(attraction$decay))
    ,error = function(e) stop("Attraction is invalid because 'distance' and 'decay' appear to be incompatible.  Perhaps lower the 'temperature'."))
  }
}

print.shallot.attraction <- function(x, ...) {
  cat("attraction for ",x$n.items," items\n",sep="")
  if ( ! is.null(x$distance) ) {
    cat("  ")
    print(x$distance)
  }
  if ( ! is.null(x$permutation) ) {
    cat("  ")
    print(x$permutation)
  }
  if ( ! is.null(x$decay) ) {
    cat("  ")
    print(x$decay)
  }
}

as.matrix.shallot.attraction <- function(x, ...) {
  stop("No yet implemented.")
  y <- x$ref$toArray()
  structure(y, dimnames=list(x$names,x$names))
}

# attraction(10)
# attraction(d,permutation(n.items=attr(d,"Size"),fixed=FALSE),decay.exponential(fixed=FALSE))
# attraction(d,permutation(n.items=attr(d,"Size"),fixed=FALSE),decay.exponential(fixed=TRUE))
# attraction(d,permutation(1:attr(d,"Size"),fixed=TRUE),decay.exponential(fixed=FALSE))
# attraction(d,permutation(1:attr(d,"Size"),fixed=TRUE),decay.exponential(fixed=TRUE))

# d <- dist(scale(USArrests))
# a <- attraction(d,permutation(n.items=attr(d,"Size"),fixed=FALSE),decay.exponential(fixed=FALSE))



# Ewens
ewens <- function(mass, n.items, names=paste0("c",1:n.items)) {
  n.items <- as.integer(n.items[1])
  if ( n.items < 0 ) stop("'n.items' must be nonnegative.")
  result <- list(mass=mass,n.items=n.items,names=names)
  class(result) <- "shallot.distribution.ewens"
  result
}

print.shallot.distribution.ewens <- function(x, ...) {
  cat("Ewens distribution with\n")
  cat("  ",x$n.items," items\n",sep="")
  cat("  ")
  print(x$mass)
}

# ewens(mass(fixed=FALSE),50)



# Ewens Pitman
ewens.pitman <- function(mass, discount, n.items, names=paste0("c",1:n.items)) {
  n.items <- as.integer(n.items[1])
  if ( n.items < 0 ) stop("'n.items' must be nonnegative.")
  result <- list(mass=mass,discount=discount,n.items=n.items,names=names)
  class(result) <- "shallot.distribution.ewensPitman"
  result
}

print.shallot.distribution.ewensPitman <- function(x, ...) {
  cat("Ewens Pitman distribution with\n")
  cat("  ",x$n.items," items\n",sep="")
  cat("  ")
  print(x$mass)
  cat("  ")
  print(x$discount)
}

# ewens.pitman(mass(fixed=FALSE),discount(fixed=FALSE),50)



# Ewens Attraction
ewens.attraction <- function(mass, attraction) {
  result <- list(mass=mass,attraction=attraction,n.items=attraction$n.items,names=attraction$names)
  class(result) <- "shallot.distribution.ewensAttraction"
  result
}

print.shallot.distribution.ewensAttraction <- function(x, ...) {
  cat("Ewens Attraction distribution with\n")
  cat("  ")
  print(x$mass)
  x <- capture.output(print(x$attraction))
  cat(paste("  ",x,"\n",sep=""))
}

# ewens.attraction(mass(fixed=FALSE),a)



# Ewens Pitman Attraction
ewens.pitman.attraction <- function(mass, discount, attraction) {
  result <- list(mass=mass,discount=discount,attraction=attraction,n.items=attraction$n.items,names=attraction$names)
  class(result) <- "shallot.distribution.ewensPitmanAttraction"
  result
}

print.shallot.distribution.ewensPitmanAttraction <- function(x, ...) {
  cat("Ewens Pitman Attraction distribution with\n")
  cat("  ")
  print(x$mass)
  cat("  ")
  print(x$discount)
  x <- capture.output(print(x$attraction))
  cat(paste("  ",x,"\n",sep=""))
}

# ewens.pitman.attraction(mass(fixed=FALSE),discount(fixed=FALSE),a)



# Distance dependent Chinese restaurant process (ddCRP)
ddcrp <- function(mass, attraction) {
  result <- list(mass=mass,attraction=attraction,n.items=attraction$n.items,names=attraction$names)
  class(result) <- "shallot.distribution.ddcrp"
  result
}

print.shallot.distribution.ddcrp <- function(x, ...) {
  cat("ddCRP distribution with\n")
  cat("  ")
  print(x$mass)
  x <- capture.output(print(x$attraction))
  cat(paste("  ",x,"\n",sep=""))
}

# ddcrp(mass(2.3),a)



# Distribution of the number of subsets
nsubsets.random <- function(x,n.samples) {
  n.samples <- as.integer(n.samples)
  if ( n.samples < 0 ) stop("'n.samples' must be positive.")
  if ( any( ! ( c("mass","n.items") %in% names(x) ) ) ) stop("Unrecognized distribution.")
  mass <- .massFactory(x$mass)
  n.items <- x$n.items
  if ( inherits(x,"shallot.distribution.ewens") ) {
    s$.distribution.Ewens$sampleNumberOfSubsets(n.items,mass,n.samples)
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    s$.distribution.EwensAttraction$sampleNumberOfSubsets(n.items,mass,n.samples)
  } else {
    discount <- .discountFactory(x$discount)
    if ( inherits(x,"shallot.distribution.ewensPitman") ) {
      s$.distribution.EwensPitman$sampleNumberOfSubsets(n.items,mass,discount,n.samples)
    } else if ( inherits(x,"shallot.distribution.ewensPitmanAttraction" ) ) {
      s$.distribution.EwensPitmanAttraction$sampleNumberOfSubsets(n.items,mass,discount,n.samples)
    } else stop("Unrecognized distribution.")
  }
}

# nsubsets.random(ewens.pitman.attraction(mass(fixed=FALSE),discount(fixed=FALSE),a),1000)



# Probability of the number of subsets
nsubsets.probability <- function(x,n.subsets) {
  n.subsets <- as.integer(n.subsets)
  if ( any( ! ( c("mass","n.items") %in% names(x) ) ) ) stop("Unrecognized distribution.")
  if ( ! x$mass$fixed ) stop("'mass' must be fixed for this function, but emperical estimates are available through the 'nsubsets.random' function.")
  mass <- .mass(x$mass)
  n.items <- x$n.items
  if ( inherits(x,"shallot.distribution.ewens") ) {
    s$.distribution.Ewens$probabilityNumberOfSubsets(n.items,n.subsets,mass)
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    s$.distribution.EwensAttraction$probabilityNumberOfSubsets(n.items,n.subsets,mass)
  } else {
    if ( ! x$discount$fixed ) stop("'discount' must be fixed for this function, but emperical estimates are available through the 'nsubsets.random' function.")
    discount <- .discount(x$discount)
    if ( inherits(x,"shallot.distribution.ewensPitman") ) {
      s$.distribution.EwensPitman$probabilityNumberOfSubsets(n.items,n.subsets,mass,discount)
    } else if ( inherits(x,"shallot.distribution.ewensPitmanAttraction") ) {
      s$.distribution.EwensPitmanAttraction$probabilityNumberOfSubsets(n.items,n.subsets,mass,discount)
    } else stop("Unrecognized distribution.")
  }
}

# nsubsets.probability(ewens.pitman.attraction(mass(n.items=a$n.items,fixed=TRUE),discount(0.1,fixed=TRUE),a),4)



# Expected number of subsets
nsubsets.average <- function(x) {
  if ( any( ! ( c("mass","n.items") %in% names(x) ) ) ) stop("Unrecognized distribution.")
  if ( ! x$mass$fixed ) stop("'mass' must be fixed for this function, but emperical estimates are available through the 'nsubsets.random' function.")
  mass <- .mass(x$mass)
  n.items <- x$n.items
  if ( inherits(x,"shallot.distribution.ewens") ) {
    s$.distribution.Ewens$meanNumberOfSubsets(n.items,mass)
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    s$.distribution.EwensAttraction$meanNumberOfSubsets(n.items,mass)
  } else {
    if ( ! x$discount$fixed ) stop("'discount' must be fixed for this function, but emperical estimates are available through the 'nsubsets.random' function.")
    discount <- .discount(x$discount)
    if ( inherits(x,"shallot.distribution.ewensPitman") ) {
      s$.distribution.EwensPitman$meanNumberOfSubsets(n.items,mass,discount)
    } else if ( inherits(x,"shallot.distribution.ewensPitmanAttraction") ) {
      s$.distribution.EwensPitmanAttraction$meanNumberOfSubsets(n.items,mass,discount)
    } else stop("Unrecognized distribution.")
  }
}

# nsubsets.average(ewens.pitman.attraction(mass(2,fixed=TRUE),discount(0.1,fixed=TRUE),a))



# Variance of the number of subsets
nsubsets.variance <- function(x) {
  if ( any( ! ( c("mass","n.items") %in% names(x) ) ) ) stop("Unrecognized distribution.")
  if ( ! x$mass$fixed ) stop("'mass' must be fixed for this function, but emperical estimates are available through the 'nsubsets.random' function.")
  mass <- .mass(x$mass)
  n.items <- x$n.items
  if ( inherits(x,"shallot.distribution.ewens") ) {
    s$.distribution.Ewens$varianceNumberOfSubsets(n.items,mass)
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    s$.distribution.EwensAttraction$varianceNumberOfSubsets(n.items,mass)
  } else {
    if ( ! x$discount$fixed ) stop("'discount' must be fixed, but emperical estimates are available through the 'nsubsets.random' function.")
    discount <- .discount(x$discount)
    if ( inherits(x,"shallot.distribution.ewensPitman") ) {
      stop("Unsupported distribution, but emperical estimates are available through the 'nsubsets.random' function.")
    } else if ( inherits(x,"shallot.distribution.ewensPitmanAttraction") ) {
      stop("Unsupported distribution, but emperical estimates are available through the 'nsubsets.random' function.")
    } else stop("Unrecognized distribution.")
  }
}

# nsubsets.variance(ewens.attraction(mass(2,fixed=TRUE),a))



# Sample for partition distributions.
.ewens <- function(x, samplingModel=.nullModel()) {
  mass <- .mass(x$mass)
  s$.distribution.Ewens$apply(samplingModel,mass,.AS.REFERENCE=TRUE)
}

.sample.ewens <- function(nItems=0L, massFactory=scalaNull('() => parameter.Mass')) s %.!% '
  val samplingModel = parameter.NullSamplingModel
  val partitionModelFactory = distribution.Ewens.factory(samplingModel,massFactory)
  distribution.PartitionModel.forwardSampler(nItems,partitionModelFactory)
'

.ewensPitman <- function(x, samplingModel=.nullModel()) {
  mass <- .mass(x$mass)
  discount <- .discount(x$discount)
  s$.distribution.EwensPitman$apply(samplingModel,mass,discount,.AS.REFERENCE=TRUE)
}

.sample.ewensPitman <- function(nItems=0L, massFactory=scalaNull('() => parameter.Mass'), discountFactory=scalaNull('() => parameter.Discount')) s %.!% '
  val samplingModel = parameter.NullSamplingModel
  val partitionModelFactory = distribution.EwensPitman.factory(samplingModel,massFactory,discountFactory)
  distribution.PartitionModel.forwardSampler(nItems,partitionModelFactory)
'

.ewensAttraction <- function(x, samplingModel=.nullModel()) {
  mass <- .mass(x$mass)
  attraction <- .attraction(x$attraction)
  s$.distribution.EwensAttraction$apply(samplingModel,mass,attraction,.AS.REFERENCE=TRUE)
}

.sample.ewensAttraction <- function(nItems=0L, massFactory=scalaNull('() => parameter.Mass'), attractionFactory=scalaNull('() => distribution.Attraction')) s %.!% '
  val samplingModel = parameter.NullSamplingModel
  val partitionModelFactory = distribution.EwensAttraction.factory(samplingModel,massFactory,attractionFactory)
  distribution.PartitionModel.forwardSampler(nItems,partitionModelFactory)
'

.ewensPitmanAttraction <- function(x, samplingModel=.nullModel()) {
  mass <- .mass(x$mass)
  discount <- .discount(x$discount)
  attraction <- .attraction(x$attraction)
  s$.distribution.EwensPitmanAttraction$apply(samplingModel,mass,discount,attraction,.AS.REFERENCE=TRUE)
}

.sample.ewensPitmanAttraction <- function(nItems=0L, massFactory=scalaNull('() => parameter.Mass'), discountFactory=scalaNull('() => parameter.Discount'), attractionFactory=scalaNull('() => distribution.Attraction')) s %.!% '
  val samplingModel = parameter.NullSamplingModel
  val partitionModelFactory = distribution.EwensPitmanAttraction.factory(samplingModel,massFactory,discountFactory,attractionFactory)
  distribution.PartitionModel.forwardSampler(nItems,partitionModelFactory)
'

.partitionsToMatrix <- function(x=scalaNull('List[parameter.partition.Partition[PersistentReference]]')) s %!% '
  x.map(_.toLabels).toArray
'

.partitionsToMatrixWithParameters <- function(x=scalaNull('List[parameter.partition.Partition[PersistentReference]]')) s %!% '
  val labelsWithParameters = x.map(_.toLabelsWithParameters)
  val labels = labelsWithParameters.map(_._1).toArray
  val parameters = labelsWithParameters.map(_._2.map(_.toString)).toArray
  (labels, parameters)
'

.sampleForward <- function(nSamples=0L, rdg=scalaNull('RDG'), sampler=scalaNull('Function2[Int, RDG, List[parameter.partition.Partition[Null]]]'), parallel=TRUE) s %.!% '
  if (!parallel) sampler(nSamples, rdg)
  else {
    val nCores = Runtime.getRuntime.availableProcessors
    val nSamplesPerCore = (nSamples / nCores) + 1
    val randomGenerator = rdg.getRandomGenerator
    val rdgList = List.fill(nCores) { new RDG(randomGenerator) }
    rdgList.par.map(r => sampler(nSamplesPerCore, r)).toList.flatten
  }
'

# Sample partitions.
sample.partitions <- function(x, n.draws, parallel=TRUE) {
  forwardSampler <- if ( inherits(x,"shallot.distribution.ewens") ) {
    .sample.ewens(x$n.items,.massFactory(x$mass))
  } else if ( inherits(x,"shallot.distribution.ewensPitman") ) {
    .sample.ewensPitman(x$n.items,.massFactory(x$mass),.discountFactory(x$discount))
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    .sample.ewensAttraction(x$n.items,.massFactory(x$mass),.attractionFactory(x$attraction))
  } else if ( inherits(x,"shallot.distribution.ewensPitmanAttraction") ) {
    .sample.ewensPitmanAttraction(x$n.items,.massFactory(x$mass),.discountFactory(x$discount),.attractionFactory(x$attraction))
  } else stop("Unrecognized partition distribution.")
  ref <- .sampleForward(n.draws,.rdg(),forwardSampler,parallel)
  structure(list(ref=ref, names=x$names), class="shallot.samples.raw")
}

print.shallot.samples.raw <- function(x, ...) {
  cat("raw partition samples --- use the 'process.partition' function to extract information\n")
}



# Posterior simulation via MCMC.
sample.partitions.posterior <- function(partition, sampling.model, partition.model, n.draws, progress.bar=interactive()) {
  sampler <- function(p=scalaNull("parameter.partition.Partition[PersistentReference]"),
                      sm=scalaNull("parameter.SamplingModel[PersistentReference]"),
                      pm=scalaNull("distribution.EwensPitmanAttraction[PersistentReference]"),
                      rdg=scalaNull("RDG"),
                      progressBar=NULL, showProgressBar=TRUE) s %.!% '
    val nDraws = R.getI0("n.draws")
    val monitor = mcmc.AcceptanceRateMonitor()
    var samples = new scala.collection.mutable.ListBuffer[parameter.partition.Partition[PersistentReference]]()
    var partition = p
    var counter = 0
    for ( i <- 1 to nDraws ) {
      partition = mcmc.AuxiliaryGibbsSampler(partition, sm, pm, rdg)._1
      samples += partition
      if ( showProgressBar && ( 100*i % nDraws == 0 ) ) {
        counter += 1
        R.invoke("setTxtProgressBar",progressBar,counter)
      }
    }
    samples.toList
  '
  sm <- .samplingModel(sampling.model)
  pm <- .partitionModel(partition.model, sm)
  p <- .labels2partition(partition, sm)
  rdg <- .rdg()
  pb <- if ( progress.bar ) txtProgressBar(min=0, max=100, style=3) else NULL
  ref <- sampler(p,sm,pm,rdg,pb,progress.bar)
  if ( progress.bar ) close(pb)
  structure(list(ref=ref, names=partition.model$names), class="shallot.samples.raw")
}



.partitionModel <- function(x, samplingModel=.nullModel()) {
  if ( inherits(x,"shallot.distribution.ewens") ) {
    .ewens(x, samplingModel)
  } else if ( inherits(x,"shallot.distribution.ewensPitman") ) {
    .ewensPitman(x, samplingModel)
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    .ewensAttraction(x, samplingModel)
  } else if ( inherits(x,"shallot.distribution.ewensPitmanAttraction") ) {
    .ewensPitmanAttraction(x, samplingModel)
  } else stop("Unrecognized partition distribution.")
}



# Make probability mass function.
partition.pmf <- function(x) {
  distribution <- .partitionModel(x)
  pmf <- distribution$logProbability(scalaNull("parameter.partition.Partition[PersistentReference]"),.EVALUATE=FALSE)
  function(x, log=TRUE) {
    partition <- if ( is.vector(x) ) .labels2partition(x,.nullModel())
    else if ( is.list(x) ) .partition2partition(x)
    else stop("'x' should be: i. a vector of cluster labels, or ii. a list containing partitions.")
    v <- pmf(partition)
    if ( log ) v else exp(v)
  }
}



# Serialize partitions to R.
serializePartitions <- function(ref, as.matrix, sample.parameter) {
  withParameters <- ( ! is.null(sample.parameter) ) && ( ref[['type']] != "List[org.ddahl.shallot.parameter.partition.Partition[Null]]" )
  if ( withParameters ) {
    zandp <- .partitionsToMatrixWithParameters(ref)
    z <- zandp$"_1"()
    p <- if ( withParameters && is.function(sample.parameter) ) apply(z,1,function(zz) lapply(1:max(zz),function(i) sample.parameter()))
    else {
      extractor <- function(x=scalaNull("Array[Array[String]]"),i=0L) s %!% 'x(i-1)'
      pp <- zandp$"_2"()
      ppp <- vector(mode="list", length=n.draws)
      for ( i in 1:n.draws ) {
        ppp[[i]] <- sapply(extractor(pp,i), function(x) s$var(x), USE.NAMES=FALSE)
      }
      ppp
    }
  } else {
    z <- .partitionsToMatrix(ref)
    p <- NULL
  }
  if ( as.matrix) {
    z <- z+1L
    r <- list(labels=z,parameters=p)
    structure(r, class="shallot.samples.labels")
  } else {
    n.draws <- nrow(z)
    r <- vector(mode="list", length=n.draws)
    for ( i in 1:n.draws ) {
      zz <- z[i,]
      n.subsets <- max(zz)+1
      rr <- vector(mode="list", length=n.subsets)
      for ( j in 1:n.subsets ) {
        items <- which(zz==j-1)
        rr[[j]] <- if ( is.null(p) ) list(items=items, parameter=NULL)
        else list(items=items,parameter=p[[i]][[j]])
      }
      r[[i]] <- rr
    }
    structure(r, class="shallot.samples.partition")
  }
}



# Sampling model
sampling.model <- function(sample.parameter, log.density) {
  if ( ! is.function(sample.parameter) ) stop("'sample.parameter' should be a function.")
  if ( length(formals(sample.parameter)) != 2 ) stop("'sample.parameter' should take two arguments named 'indices' and 'parameter'.")
  if ( ! identical(names(formals(sample.parameter)),c("indices","parameter")) )
    stop("'sample.parameter' should take two arguments named 'indices' and 'parameter'.")
  if ( ! is.function(log.density) ) stop("'log.density' should be a function.")
  if ( length(formals(log.density)) != 3 ) stop("'log.density' should take three arguments named 'i', 'indices', and 'parameter'.")
  if ( ! identical(names(formals(log.density)),c("i","indices","parameter")) )
    stop("'log.density' should take three arguments named 'i', 'indices', and 'parameter'.")
  tryCatch(p <- sample.parameter(), error=function(e) stop("'sample.parameter' should sample from the centering distribution when no arguments are provided."))
  tryCatch(d <- log.density(1,1,p), error=function(e) stop("'log.density' should doesn't pass sanity check log.density(1,1,p)."))
  if ( ( ! is.vector(d) ) || ( length(d) != 1 ) ) stop("'log.density' should return a numeric vector of length one.")
  r <- list(sample.parameter=sample.parameter, log.density=log.density)
  structure(r, class="shallot.distribution.data")
}

.samplingModel <- function(samplingModel=NULL) {
  s %.!% '
    val sp = R.evalReference(samplingModel+"$sample.parameter")
    val ld = R.evalReference(samplingModel+"$log.density")

    new parameter.SamplingModel[PersistentReference] {

      def logDensity(i: Int, subset: parameter.partition.Subset[PersistentReference]): Double = {
        R.invokeD0(ld, i+1, subset.toArray.map(_+1), subset.parameter)
      }

      def sample(subset: parameter.partition.Subset[PersistentReference]): PersistentReference = {
        R.invokeReference(sp, subset.toArray.map(_+1), subset.parameter)
      }

      def sample(): PersistentReference = {
        R.invokeReference(sp)
      }

    }
  '
}



# Process partitions that were sampled.
process.partitions <- function(x, as.matrix=TRUE, expand=FALSE, sample.parameter=FALSE) {
  if ( ! inherits(x,"shallot.samples.raw") ) stop("'x' should be a result from the 'sample.partition' function.")
  result <- serializePartitions(x$ref, as.matrix, sample.parameter)
  if ( expand ) {
    if ( identical(sample.parameter,NULL) ) stop("'sample.parameter' may not be null when 'expand=TRUE'.")
    if ( ! as.matrix ) stop("'expand=TRUE' is only valid when 'as.matrix=TRUE'.")
    mega <- matrix(NA,nrow=nrow(result$labels),ncol=ncol(result$labels))
    for ( i in 1:nrow(mega) ) {
      rhs <- result$parameters[[i]]
      if ( ! is.vector(rhs) ) stop("'expand=TRUE' is only valid when all the parameters are scalars.")
      mega[i, ] <- rhs[result$labels[i,]]
    }
    mega
  } else result
}



# Null sampling model
.nullModel <- function() s %.!% 'new parameter.GeneralNullSamplingModel[PersistentReference]()'



# Pairwise Probabilities
enumerate.partitions <- function(n.items, as.matrix=TRUE) {
  ref <- s$.parameter.partition.Partition$enumerate(.nullModel(),as.integer(n.items)[1])
  serializePartitions(ref, as.matrix=as.matrix, sample.parameter=NULL)
}



# Pairwise Probabilities
pairwise.probabilities <- function(x, parallel=TRUE) {
  if ( ! inherits(x,"shallot.samples.raw") ) stop("'x' should be a result from the 'sample.partition' function.")
  start.time <- proc.time()
  ref <- s$.parameter.partition.PairwiseProbability$apply(x$ref,as.logical(parallel))
  result <- list(ref=ref,n.items=ref$nItems(),names=x$names,proc.time=proc.time()-start.time)
  structure(result, class="shallot.pairwiseProbability")
}

print.shallot.pairwiseProbability <- function(x, ...) {
  cat("pairwise probabilities for ",x$n.items," items --- use 'as.matrix' function to obtain matrix\n",sep="")
}

as.matrix.shallot.pairwiseProbability <- function(x, ...) {
  structure(x$ref$toArray(), dimnames=list(x$names,x$names))
}



# Estimate partition
estimate.partition <- function(x, pairwise.probabilities=NULL, max.subsets=0, max.scans=0, parallel=TRUE) {
  if ( ! inherits(x,"shallot.samples.raw") ) stop("'x' should be a result from the 'sample.partition' function.")
  if ( is.null(pairwise.probabilities) ) pairwise.probabilities <- pairwise.probabilities(x)
  ref <- s$.parameter.partition.MinBinder$apply(
      pairwise.probabilities$ref, x$ref, as.integer(max.subsets), as.integer(max.scans), as.logical(parallel))
  structure(ref$toLabels(), names=x$names)
}



# Confidence
.labels2partition <- function(partition=integer(), samplingModel=scalaNull("parameter.SamplingModel[PersistentReference]")) s %.!% '
  parameter.partition.Partition(samplingModel,partition)
'

.partition2partition <- function(partition) {
  stop("Not yet implemented.")
}

confidence <- function(pairwise.probabilities, partition) {
  tmpObj <- pairwise.probabilities$ref$confidenceComputations(.labels2partition(partition,.nullModel()))
  partition <- tmpObj$"_1"() + 1
  names(partition) <- pairwise.probabilities$names
  confidence <- tmpObj$"_2"()
  names(confidence) <- names(partition)
  confidence.matrix <- tmpObj$"_3"()
  dimnames(confidence.matrix) <- list(1:nrow(confidence.matrix),1:ncol(confidence.matrix))
  order <- tmpObj$"_4"() + 1L
  names(order) <- names(partition)
  exemplar <- tmpObj$"_5"() + 1L
  names(exemplar) <- 1:length(exemplar)
  result <- list(partition=partition,confidence=confidence,confidence.matrix=confidence.matrix,exemplar=exemplar,order=order,pairwise.probabilities=pairwise.probabilities)
  class(result) <- "shallot.confidence"
  result
}



# Confidence or pairs plot
.rotateForConfidencePlot <- function(pp=scalaNull("parameter.partition.PairwiseProbability"), order=integer()) s %!% '
  val nItems = pp.nItems
  val xx = Array.ofDim[Double](nItems, nItems)
  for (i <- 0 until nItems) {
    for (j <- 0 until nItems) {
      xx(i)(nItems - j - 1) = pp(order(i) - 1, order(j) - 1)
    }
  }
  xx
'

plot.shallot.confidence <- function(x, partition=NULL, data=NULL, show.labels=length(x$partition)<=50, ...) {
  if ( ! is.null(data) ) {
    if ( ! is.null(partition) ) stop("'partition' must be 'NULL' for pairs plot.")
    i <- x$exemplar[x$partition]
    c <- rainbow(length(x$exemplar))[x$partition]
    panelFnc <- function(x0,y0,...) {
      points(x0,y0,col=c,pch=19,...)
      segments(x0,y0,x0[i],y0[i],col=c,...)
      points(x0[x$exemplar],y0[x$exemplar],pch=22,bg="white",cex=2,...)
    }
    pairs(data,panel=panelFnc)
    return(invisible())
  }
  if ( is.null(partition) ) {
    partition <- x$partition
    o <- x$order
  } else {
    o <- order(partition)
  }
  pm <- .rotateForConfidencePlot(x$pairwise.probabilities$ref,o)
  n <- nrow(pm)
  sizes <- rle(partition[o])$lengths
  cuts <- cumsum(sizes)
  centers <- ( c(0,cuts[-length(cuts)]) + cuts ) / 2
  cuts <- cuts[-length(cuts)]
  labels <- rle(partition[o])$values
  if ( show.labels ) {
    mymai <- c(1.5,0.5,0.5,1.5)
    cexscale <- 0.85 * 50 / length(partition)
  } else {
    mymai <- c(0,0,0,0)
    cexscale <- 1 * 50 / length(partition)
  }
  opar <- par(pty="s",mai=mymai)
  colors <- topo.colors(200)
  colors <- rev(heat.colors(200))
  image(x=1:n,y=1:n,z=pm,axes=FALSE,xlab="",ylab="",col=colors)
  box()
  abline(v=cuts+0.5,lwd=3)
  abline(h=n-cuts+0.5,lwd=3)
  text(centers+0.5,n-centers+0.5,labels,cex=0.8*cexscale*sizes)
  if ( show.labels ) {
    axisLabels <- if ( is.null(names(partition)) ) o
    else names(partition[o])
    axis(4,1:length(partition),rev(axisLabels),las=2,cex.axis=0.8*cexscale)
    axis(1,1:length(partition),axisLabels,las=2,cex.axis=0.8*cexscale)
    nn <- length(colors)
    bx <- par("usr")
    bx.cx <- c(bx[1] - 1.6 * (bx[2] - bx[1]) / 50, bx[1] - 0.3 * (bx[2] - bx[1]) / 50)
    bx.cy <- c(bx[3], bx[3])
    bx.sy <- (bx[4] - bx[3]) / nn
    xx <- rep(bx.cx, each=2)
    for ( i in 1:nn ) {
      yy <- c(bx.cy[1] + (bx.sy * (i - 1)),
              bx.cy[1] + (bx.sy * (i)),
              bx.cy[1] + (bx.sy * (i)),
              bx.cy[1] + (bx.sy * (i - 1)))
      polygon(xx,yy,col=colors[i],border=colors[i],xpd=TRUE)
    }
  }
  par(opar)
  invisible()
}

