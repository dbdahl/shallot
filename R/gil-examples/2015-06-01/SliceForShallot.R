
# R FUNCTIONS FOR PERFORMING UNIVARIATE SLICE SAMPLING.
#
# Radford M. Neal, 17 March 2008.
#
# Fellingham modified to interface with Shallot, 18 February 2015, 23 March 2015
#
# Fellingham modified 29 May,2015 to deal with some covariates the same across
# all individuals, and some covriates clustered
# This modification assumes 115 individuals
#
#
# Implements, with slight modifications and extensions, the algorithm described
# in Figures 3 and 5 of the following paper:
#
#   Neal, R. M (2003) "Slice sampling" (with discussion), Annals of Statistics,
#      vol. 31, no. 3, pp. 705-767.
#
# See the documentation for the function uni.slice below for how to use it.
# The function uni.slice.test was used to test the uni.slice function.


# GLOBAL VARIABLES FOR RECORDING PERFORMANCE.

# uni.slice.calls <- 0  # Number of calls of the slice sampling function
# uni.slice.evals <- rep(0,length(parameter))	# Number of density evaluations done in these calls


# UNIVARIATE SLICE SAMPLING WITH STEPPING OUT AND SHRINKAGE.
#
# Performs a slice sampling update from an initial point to a new point that 
# leaves invariant the distribution with the specified log density function.
#
# Arguments:
#
#   x0    Initial point
#   g     Function returning the log of the probability density (plus constant)
#   w     Size of the steps for creating interval (default 1)
#   m     Limit on steps (default infinite)
#   lower Lower bound on support of the distribution (default -Inf)
#   upper Upper bound on support of the distribution (default +Inf)
#   gx0   Value of g(x0), if known (default is not known)
#
# The log density function may return -Inf for points outside the support 
# of the distribution.  If a lower and/or upper bound is specified for the
# support, the log density function will not be called outside such limits.
#
# The value of this function is the new point sampled, with an attribute
# of "log.density" giving the value of the log density function, g, at this
# point.  Depending on the context, this log density might be passed as the 
# gx0 argument of a future call of uni.slice. 
#
# The global variable uni.slice.calls is incremented by one for each call
# of uni.slice.  The global variable uni.slice.evals is incremented by the
# number of calls made to the g function passed.
#
# WARNING:  If you provide a value for g(x0), it must of course be correct!
# In addition to giving wrong answers, wrong values for gx0 may result in
# the uni.slice function going into an infinite loop.

uni.slice <- function (parameter ,subset, w=1, m=500, lower=0, upper=+Inf, parmNumber)
{
  # Check the validity of the arguments.
  
#  if (!is.numeric(x0) || length(x0)!=1
#      || !is.function(g) 
#      || !is.numeric(w) || length(w)!=1 || w<=0 
#      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
#      || !is.numeric(lower) || length(lower)!=1 || x0<lower
#      || !is.numeric(upper) || length(upper)!=1 || x0>upper
#      || upper<=lower 
#      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
#  { 
#    stop ("Invalid slice sampling argument")
#  }
  
  # Keep track of the number of calls made to this function.
  
#  uni.slice.calls <<- uni.slice.calls + 1
  
  # Find the log density at the initial point, if not already known.
  

  
   uni.slice.evals[parmNumber] <<- uni.slice.evals[parmNumber] + 1
    gxs <- 0
    for (i in 1:length(subset)){
      gxs <- gxs + loglike.i(subset[i],parameter[subset[i],])
    }
  
  
  # Determine the slice level, in log terms.
  
  logy <- gxs - rexp(1)
  
  # Find the initial interval to sample from.
  
  u <- runif(1,0,w)
  L <- parameter[1:115,parmNumber] - u
  R <- parameter[1:115,parmNumber] + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
  
  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.
  
  if (is.infinite(m))  # no limit on number of steps
  { 
    repeat
    { if (L<=lower) break
      uni.slice.evals[parmNumber] <<- uni.slice.evals[parmNumber] + 1

      newParm <- newParameter(parameter,parmNumber,L)
      gL <- 0
      for (i in 1:length(subset)){
        gL <- gL + loglike.i(subset[i],newParm[subset[i]])
      }

      if (gL<=logy) break
      L <- L - w
    }
    
    repeat
    { if (R>=upper) break
      uni.slice.evals[parmNumber] <<- uni.slice.evals[parmNumber] + 1
      newParm <- newParameter(parameter,parmNumber,R)
      gR <- 0
      for (i in 1:length(subset)){
        gR <- gR + loglike.i(subset[i],newParm[subset[i]])
      }
      if (gR<=logy) break
      R <- R + w
    }

  
  } 
  else # syntax question if (m>1)  # limit on steps, bigger than one
  { 
    J <- floor(runif(1,0,m))
    K <- (m-1) - J
    
    while (J>0)
    { if (L<=lower) break
      uni.slice.evals[parmNumber] <<- uni.slice.evals[parmNumber] + 1
      newParm <- newParameter(parameter,parmNumber,L)
      gL <- 0
      for (i in 1:length(subset)){
        gL <- gL + loglike.i(subset[i],newParm[subset[i]])
      }      

      if (gL<=logy) break
      L <- L - w
      J <- J - 1
    }
    
    while (K>0)
    { if (R>=upper) break
      uni.slice.evals[parmNumber] <<- uni.slice.evals[parmNumber] + 1
      newParm <- newParameter(parameter,parmNumber,R)
      gR <- 0
      for (i in 1:length(subset)){
        gR <- gR + loglike.i(subset[i],newParm[subset[i]])
      }      
      if (gR<=logy) break
      R <- R + w
      K <- K - 1
    }
  }
  
  # Shrink interval to lower and upper bounds.
  
  if (L<lower) 
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }
  
  # Sample from the interval, shrinking it on each rejection.
  
  repeat
  { 
    x1 <- runif(1,L,R)
    
    uni.slice.evals[parmNumber] <<- uni.slice.evals[parmNumber] + 1
    newParm <- newParameter(parameter,parmNumber,x1)
    gx1 <- 0
    for (i in 1:length(subset)){
      gx1 <- gx1 + loglike.i(subset[i],newParm[subset[i]])
    }    
    if (gx1>=logy) break
    
    if (x1>parameter[subset[1],parmNumber]) 
    { R <- x1
    }
    else 
    { L <- x1
    }
  }
  
  # Return the point sampled, with its log density attached as an attribute.

#  attr(x1,"log.density") <- gx1
  return (x1)
  
}

newParameter <- function(parameter,parmNumber,sub){
if (dim(parameter)[2] == 1){
  newParm[,1] <- rep(sub,115)
}
else if (parmNumber == 1){
  newParm <- c(sub,parameter[1:115,2:dim(parameter)[2]])
}
else if (parmNumber == dim(parameter)[2]){
  newParm <- c(parameter[1:115,1:(dim(parameter)[2]-1)],sub)  
}
else {
  newParm <- c(parameter[1:115,1:(parmNumber-1)],sub,parameter[1:115,(parmNumber+1):dim(parameter)[2]])
}
return(newParm)
}



