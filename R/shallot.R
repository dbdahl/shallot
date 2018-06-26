### Development
# library(rscala); library(shallot); s <- shallot:::s

##########################
#### Simplified
##########################


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
.mass <- function(mass) s$Mass(mass$value)

.massFactory <- function(mass) {
  if ( mass$fixed ) s$Mass.factory(mass$value)
  else s$Mass.factory(mass$shape,mass$rate,s$rdg())
}

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
.discount <- function(discount) s$Discount(discount$value)

.discountFactory <- function(discount) {
  if ( discount$fixed ) s$Discount.factory(discount$value)
  else s$Discount.factory(discount$shape,discount$rate,s$rdg())
}

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
  s$.Permutation.apply(permutation$value-1L)
}

.permutationFactory <- function(permutation) {
  if ( permutation$fixed ) {
    s$.Permutation.factory(permutation$value-1L)
  } else {
    s$.Permutation.factory(permutation$n.items,s$rdg())
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
       if ( decay$type == "reciprocal" )  s$.new_decay.ReciprocalDecay(temp)
  else if ( decay$type == "exponential" ) s$.new_decay.ExponentialDecay(temp)
  else if ( decay$type == "subtraction" ) s$.new_decay.SubtractionDecay(temp,decay$max.distance)
}

.decayFactory <- function(decay) {
  if ( decay$temperature$fixed ) {
    temp <- min(decay$temperature$value,decay$max.temperature)
         if ( decay$type == "reciprocal" )  s$.decay.ReciprocalDecayFactory.factory(temp)
    else if ( decay$type == "exponential" ) s$.decay.ExponentialDecayFactory.factory(temp)
    else if ( decay$type == "subtraction" ) s$.decay.SubtractionDecayFactory.new(decay$max.distance)$factory(temp)
  } else {
    shape <- decay$temperature$shape
    rate <- decay$temperature$rate
         if ( decay$type == "reciprocal" )  s$.decay.ReciprocalDecayFactory.factory(shape,rate,s$rdg())
    else if ( decay$type == "exponential" ) s$.decay.ExponentialDecayFactory.factory(shape,rate,s$rdg())
    else if ( decay$type == "subtraction" ) s$.new_decay.SubtractionDecayFactory(decay$max.distance)$factory(shape,rate,s$rdg())
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
  s$.Distance.apply(as.matrix(distance),FALSE)
}

.attraction <- function(attraction) {
  if ( attraction$constant ) {
    s$.Attraction.constant(attraction$n.items)
  } else {
    tryCatch(
      s$.Attraction.apply(.distance(as.matrix(attraction$decay$distance)),.permutation(attraction$permutation),.decay(attraction$decay))
    ,error = function(e) stop("Attraction is invalid because 'distance' and 'decay' appear to be incompatible.  Perhaps lower the 'temperature'."))
  }
}

.attractionFactory <- function(attraction) {
  if ( attraction$constant ) {
    s$.Attraction.factory(attraction$n.items)
  } else if ( attraction$permutation$fixed && attraction$decay$temperature$fixed ) {
    tryCatch(
      s$.Attraction.factory(.distance(as.matrix(attraction$decay$distance)),.permutation(attraction$permutation),.decay(attraction$decay))
    ,error = function(e) stop("Attraction is invalid because 'distance' and 'decay' appear to be incompatible.  Perhaps lower the 'temperature'."))
  } else if ( attraction$permutation$fixed ) {
    tryCatch(
      s$.Attraction.factory(.distance(as.matrix(attraction$decay$distance)),.permutation(attraction$permutation),.decayFactory(attraction$decay))
    ,error = function(e) stop("Attraction is invalid because 'distance' and 'decay' appear to be incompatible.  Perhaps lower the 'temperature'."))
  } else if ( attraction$decay$temperature$fixed ) {
    tryCatch(
      s$.Attraction.factory(.distance(as.matrix(attraction$decay$distance)),.permutationFactory(attraction$permutation),.decay(attraction$decay))
    ,error = function(e) stop("Attraction is invalid because 'distance' and 'decay' appear to be incompatible.  Perhaps lower the 'temperature'."))
  } else {
    tryCatch(
      s$.Attraction.factory(.distance(as.matrix(attraction$decay$distance)),.permutationFactory(attraction$permutation),.decayFactory(attraction$decay))
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
    s$Ewens.sampleNumberOfSubsets(n.items,mass,n.samples)
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    s$EwensAttraction.sampleNumberOfSubsets(n.items,mass,n.samples)
  } else {
    discount <- .discountFactory(x$discount)
    if ( inherits(x,"shallot.distribution.ewensPitman") ) {
      s$EwensPitman.sampleNumberOfSubsets(n.items,mass,discount,n.samples)
    } else if ( inherits(x,"shallot.distribution.ewensPitmanAttraction" ) ) {
      s$EwensPitmanAttraction.sampleNumberOfSubsets(n.items,mass,discount,n.samples)
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
    s$Ewens.probabilityNumberOfSubsets(n.items,n.subsets,mass)
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    s$EwensAttraction.probabilityNumberOfSubsets(n.items,n.subsets,mass)
  } else {
    if ( ! x$discount$fixed ) stop("'discount' must be fixed for this function, but emperical estimates are available through the 'nsubsets.random' function.")
    discount <- .discount(x$discount)
    if ( inherits(x,"shallot.distribution.ewensPitman") ) {
      s$EwensPitman.probabilityNumberOfSubsets(n.items,n.subsets,mass,discount)
    } else if ( inherits(x,"shallot.distribution.ewensPitmanAttraction") ) {
      s$EwensPitmanAttraction.probabilityNumberOfSubsets(n.items,n.subsets,mass,discount)
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
    s$Ewens.meanNumberOfSubsets(n.items,mass)
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    s$EwensAttraction.meanNumberOfSubsets(n.items,mass)
  } else {
    if ( ! x$discount$fixed ) stop("'discount' must be fixed for this function, but emperical estimates are available through the 'nsubsets.random' function.")
    discount <- .discount(x$discount)
    if ( inherits(x,"shallot.distribution.ewensPitman") ) {
      s$EwensPitman.meanNumberOfSubsets(n.items,mass,discount)
    } else if ( inherits(x,"shallot.distribution.ewensPitmanAttraction") ) {
      s$EwensPitmanAttraction.meanNumberOfSubsets(n.items,mass,discount)
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
    s$Ewens.varianceNumberOfSubsets(n.items,mass)
  } else if ( inherits(x,"shallot.distribution.ewensAttraction") ) {
    s$EwensAttraction.varianceNumberOfSubsets(n.items,mass)
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
  s$.Ewens.apply(samplingModel,mass)
}

.sample.ewens <- function(nItems, massFactory) s(nItems=as.integer(nItems[1]),massFactory=massFactory) ^ '
  val samplingModel = NullSamplingModel
  val partitionModelFactory = Ewens.factory(samplingModel,massFactory)
  PartitionModel.forwardSampler(nItems,partitionModelFactory)
'

.ewensPitman <- function(x, samplingModel=.nullModel()) {
  mass <- .mass(x$mass)
  discount <- .discount(x$discount)
  s$.EwensPitman.apply(samplingModel,mass,discount)
}

.sample.ewensPitman <- function(nItems, massFactory, discountFactory) s(nItems=as.integer(nItems[1]),massFactory=massFactory,discountFactory=discountFactory) ^ '
  val samplingModel = NullSamplingModel
  val partitionModelFactory = EwensPitman.factory(samplingModel,massFactory,discountFactory)
  PartitionModel.forwardSampler(nItems,partitionModelFactory)
'

.ewensAttraction <- function(x, samplingModel=.nullModel()) {
  mass <- .mass(x$mass)
  attraction <- .attraction(x$attraction)
  s$.EwensAttraction.apply(samplingModel,mass,attraction)
}

.sample.ewensAttraction <- function(nItems, massFactory, attractionFactory) s(nItems=as.integer(nItems[1]),massFactory=massFactory,attractionFactory=attractionFactory) ^ '
  val samplingModel = NullSamplingModel
  val partitionModelFactory = EwensAttraction.factory(samplingModel,massFactory,attractionFactory)
  PartitionModel.forwardSampler(nItems,partitionModelFactory)
'

.ewensPitmanAttraction <- function(x, samplingModel=.nullModel()) {
  mass <- .mass(x$mass)
  discount <- .discount(x$discount)
  attraction <- .attraction(x$attraction)
  s$.EwensPitmanAttraction.apply(samplingModel,mass,discount,attraction)
}

.sample.ewensPitmanAttraction <- function(nItems, massFactory, discountFactory, attractionFactory) s(nItems=as.integer(nItems[1]), massFactory=massFactory, discountFactory=discountFactory, attractionFactory=attractionFactory) ^ '
  val samplingModel = NullSamplingModel
  val partitionModelFactory = EwensPitmanAttraction.factory(samplingModel,massFactory,discountFactory,attractionFactory)
  PartitionModel.forwardSampler(nItems,partitionModelFactory)
'

.partitionsToMatrix <- function(x) s(x=x) * '
  x.map(_.toLabels).toArray
'

.partitionsToMatrixWithParameters <- function(x) s(x=x) ^ '
  val labelsWithParameters = x.map(_.toLabelsWithParameters)
  val labels = labelsWithParameters.map(_._1).toArray
  val parameters = labelsWithParameters.map(_._2.map(_.toString)).toArray
  (labels, parameters)
'

.sampleForward <- function(nSamples, rdg, sampler, parallel=TRUE) {
  s(nSamples=as.integer(nSamples), sampler=sampler, parallel=parallel) ^ '
    if (!parallel) sampler(nSamples, rdg)
    else {
      val nCores = Runtime.getRuntime.availableProcessors
      val nSamplesPerCore = (nSamples / nCores) + 1
      val randomGenerator = rdg.getRandomGenerator
      val rdgList = List.fill(nCores) { new RDG(randomGenerator) }
      rdgList.par.map(r => sampler(nSamplesPerCore, r)).toList.flatten
    }
  '
}

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
  ref <- .sampleForward(n.draws,s$rdg(),forwardSampler,parallel)
  structure(list(ref=ref, names=x$names), class="shallot.samples.raw")
}

print.shallot.samples.raw <- function(x, ...) {
  cat("raw partition samples --- use the 'process.samples' function to extract information\n")
}

print.shallot.samples.full <- function(x, ...) {
  cat("raw partition samples with hyperparameter values --- use the 'process.samples' function to extract information\n")
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
  pmf <- distribution$logProbability
  function(x, log=TRUE) {
    partition <- if ( is.vector(x) ) .labels2partition(x,.nullModel())
    else if ( is.list(x) ) .partition2partition(x)
    else stop("'x' should be: i. a vector of cluster labels, or ii. a list containing partitions.")
    v <- pmf(partition)
    if ( log ) v else exp(v)
  }
}



# Serialize partitions to R.
serializePartitions <- function(ref, names, as.matrix, sample.parameter) {
  withParameters <- ( ! is.null(sample.parameter) ) && ( ref[['type']] != "List[org.ddahl.shallot.parameter.partition.Partition[Null]]" )
  if ( withParameters ) {
    zandp <- .partitionsToMatrixWithParameters(ref)
    z <- zandp$"_1"()
    colnames(z) <- names
    p <- if ( withParameters && is.function(sample.parameter) ) apply(z,1,function(zz) lapply(1:max(zz),function(i) sample.parameter()))
    else {
      extractor <- function(x,i=0L) s(x=x) * 'x(i-1)'
      pp <- zandp$"_2"()
      n.draws <- nrow(z)
      ppp <- vector(mode="list", length=n.draws)
      for ( i in 1:n.draws ) {
        ppp[[i]] <- sapply(extractor(pp,i), function(x) s$var(x), USE.NAMES=FALSE)
      }
      ppp
    }
  } else {
    z <- .partitionsToMatrix(ref)
    colnames(z) <- names
    p <- NULL
  }
  if ( as.matrix) {
    z <- z+1L
    r <- list(labels=z,parameters=p)
    list(partitions=structure(r, class="shallot.samples.labels"))
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
    list(partitons=structure(r, class="shallot.samples.partition"))
  }
}




# Process partitions that were sampled.
process.samples <- function(x, as.matrix=TRUE, expand=FALSE, sample.parameter=FALSE) {
  if ( inherits(x,"shallot.samples.full") ) {
    result <- process.samples(x$raw, as.matrix=as.matrix, expand=expand, sample.parameter=sample.parameter)
    result[['hyperparameters']] <- x$hyperparameters
    return(result)
  }
  if ( ! inherits(x,"shallot.samples.raw") ) stop("'x' should be a result from the functions 'sample.partition' or 'sample.partition.posterior'.")
  result <- serializePartitions(x$ref, x$names, as.matrix, sample.parameter)
  if ( expand ) {
    if ( identical(sample.parameter,NULL) ) stop("'sample.parameter' may not be null when 'expand=TRUE'.")
    if ( ! as.matrix ) stop("'expand=TRUE' is only valid when 'as.matrix=TRUE'.")
    mega <- matrix(NA,nrow=nrow(result$partitions$labels),ncol=ncol(result$partitions$labels))
    for ( i in 1:nrow(mega) ) {
      rhs <- result$partitions$parameters[[i]]
      if ( ! is.vector(rhs) ) stop("'expand=TRUE' is only valid when all the parameters are scalars.")
      mega[i, ] <- rhs[result$partitions$labels[i,]]
    }
    colnames(mega) <- x$names
    list(partitions=mega)
  } else result
}



# Null sampling model
.nullModel <- function() s$'.new_GeneralNullSamplingModel[PersistentReference]'()



# Pairwise Probabilities
enumerate.partitions <- function(n.items, as.matrix=TRUE) {
  ref <- s$.Partition.enumerate(.nullModel(),as.integer(n.items)[1])
  serializePartitions(ref, NULL, as.matrix=as.matrix, sample.parameter=NULL)
}



# Pairwise Probabilities
pairwise.probabilities <- function(x, parallel=TRUE) {
  if ( inherits(x,"shallot.samples.full") ) x <- x$raw
  if ( ! inherits(x,"shallot.samples.raw") ) stop("'x' should be a result from the functions 'sample.partition' or 'sample.partition.posterior'.")
  start.time <- proc.time()
  ref <- s$.PairwiseProbability.apply(x$ref,as.logical(parallel))
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
  if ( inherits(x,"shallot.samples.full") ) x <- x$raw
  if ( ! inherits(x,"shallot.samples.raw") ) stop("'x' should be a result from the functions 'sample.partition' or 'sample.partition.posterior'.")
  if ( is.null(pairwise.probabilities) ) pairwise.probabilities <- pairwise.probabilities(x)
  ref <- s$.MinBinder.apply(
      pairwise.probabilities$ref, x$ref, as.integer(max.subsets), as.integer(max.scans), as.logical(parallel))
  structure(ref$toLabels(), names=x$names)
}


# Misc
.labels2partition <- function(partition, samplingModel) {
  s(partition=as.integer(partition),samplingModel=samplingModel) ^ '
    Partition((i: Int) => samplingModel.sample(),partition)
  '
}

.partition2partition <- function(partition) {
  stop("Not yet implemented.")
}

latest <- function() {
  install.packages('https://dahl.byu.edu/public/rscala_LATEST.tar.gz',repos=NULL,type='source')
  install.packages('https://dahl.byu.edu/public/commonsMath_LATEST.tar.gz',repos=NULL,type='source')
  install.packages('https://dahl.byu.edu/public/sdols_LATEST.tar.gz',repos=NULL,type='source')
  install.packages('https://dahl.byu.edu/public/shallot_LATEST.tar.gz',repos=NULL,type='source')
  cat("\n\n*******\n\nPlease restart your R session.\n\n*******\n")

}

