# compute variance ratio
var.ratio <- function(conf) {
  cl.am <- assoc.matrix(conf$clustering)
  epam <- conf$expectedPairwiseAllocationMatrix
  cl.vars <- sapply(unique(conf$clustering), function(i) {
    if (sum(conf$clustering == i) < 3) return(c(0,1))
    cl.am[conf$clustering != i,conf$clustering != i] <- 0
    x <- epam[upper.tri(epam) & cl.am]
    c(sum((mean(x) - x)^2)/(length(x)-1),length(x))
  })
  
  cl.var <- sum(apply(cl.vars,2,prod)/sum(cl.vars[2,]))
  y <- epam[upper.tri(epam)]
  tot.var <- sum((mean(y) - y)^2)/(length(y) - 1)
  
  cl.var/tot.var
}

part.conf <- function(conf){
  if (length(unique(conf$clustering)) == length(conf$clustering)) return(0)
  cl.am <- assoc.matrix(conf$clustering)
  epam <- conf$expectedPairwiseAllocationMatrix
  mean(epam[upper.tri(epam) & cl.am])
}

# algorithm to select mass value given:
# mass values
# parition confidence (pc) 
# variance ratio      (vr)
# number of clusters  (n)
# weights correspond to pc, vr, n
mass.algorithm <- function(mass,pc,vr,n,w=c(1,1,1),two.stage=TRUE) {
  if (length(w) != 3) stop("Weights must be assigned to: PC VR N")
  
  compute.rankings <- function(subset=1:length(pc)) {
    rankings <- rank(1-pc[subset])*w[1] + rank(vr[subset])*w[2] + n[subset]*w[3]
    if (length(rankings) > 3 && FALSE) {
      best <- which.min(sapply(1:(length(rankings)-2), function(i) sum(rankings[i:(i+2)]))) + 1
    } else {
      best <- which.min(rankings)
    }
    
    c(n=n[subset][best], mass = mass[subset][best], pc=pc[subset][best], vr=vr[subset][best])
  }
  
  if (two.stage) {
    all <- sapply(unique(n[n>1]), function(i) compute.rankings(n == i))
    rankings <- rank(1-all["pc",])*w[1] + rank(all["vr",])*w[2] + all["n",]*w[3]
    index <- which.min(rankings)
    out <- t(all[,(index-1):(index+1)])
  } else {
    best <- compute.rankings(n>1)
    best.nm1 <- if (best["n"] > 2) compute.rankings(n == (best["n"]-1)) else best
    best.np1 <- compute.rankings(n == (best["n"]+1))
    
    out <- rbind(best.nm1,best,best.np1)
    rownames(out) <- NULL
  }
  
  plot(mass,pc,type="l",ylim=c(0,1),col="dodgerblue")
  lines(vr[n>1]~mass[n>1],col="forestgreen")
  if (nrow(out) == 2) {
    abline(v=out[1,"mass"])
    abline(v=out[2,"mass"],lty=2)
  } else {
    abline(v=out[1,"mass"],lty=2)
    abline(v=out[2,"mass"])
    abline(v=out[3,"mass"],lty=2)
  }
  
  out
}

# list.epam, list of EPAMs
# new.draws - generates draws at specified mass values
default.mass <- function(mass, list.epam, new.draws = TRUE, w=c(1,1,1), discount=0, dis, temp=10, loss="binder", m=100L, parallel=TRUE) {
  if (missing(mass)) {
    if (loss == "lowerBoundVariationOfInformation") {
      mass <- seq(0.1,10,0.2)
    } else {
      mass <- seq(0.1,5,by=0.05)
    }
  } 
  
  single.epam <- function(mass,discount,dis,m,temp=10,parallel=parallel) {
    n <- nrow(as.matrix(dis))
    epa <- ewens.pitman.attraction(mass(mass), discount(discount), attraction(permutation(n.items=n, fixed=FALSE), decay.exponential(temperature(temp), dis)))
    draws <- sample.partitions(epa, m, parallel = parallel)
    as.matrix(pairwise.probabilities(draws))
  }
  
  if(missing(list.epam) && new.draws == FALSE) stop("must provide list.epam or obtain new draws")
  if (missing(list.epam) || new.draws == TRUE) {
    out <- lapply(mass, single.epam,  m=m, discount=discount, dis=dis,  temp=temp, parallel=parallel)
    sapply(1:length(out), function(i) {
      attr(out[[i]],"mass") <<- mass[i]
      attr(out[[i]],"m") <<- m
      attr(out[[i]],"discount") <<- discount
    })
    if (!missing(list.epam)) {
      mass0 <- sapply(list.epam,attr,"mass")
      combine <- which(mass0 %in% mass)
      add <- which(!(mass %in% mass0))
      if (length(combine) > 0) {
        sapply(combine, function(i) {
          m0 <- attr(list.epam[[i]],"m")
          which.out <- which(mass0[i] == mass)
          list.epam[[i]] <<- (list.epam[[i]]*m0 + out[[which.out]]*m)/(m + m0)
          attr(list.epam[[i]],"m") <<- m0 + m
        })
      }
      if (length(add) > 0) {
        list.epam <- c(list.epam, out[add])
      }
      list.epam <- list.epam[order(sapply(list.epam,attr,"mass"))]
    } else {
      list.epam <- out
    }
  }
  
  out.conf <- lapply(list.epam, function(x) {
    attr(x,"mass") <- NULL
    attr(x,"discount") <- NULL
    attr(x,"m") <- NULL
    cl <- salso(x,loss=loss,multicore=parallel)
    confidence(cl,x)
  })
  
  mass <- sapply(list.epam, attr, "mass")
  weight <- sapply(list.epam, attr, "weight")
  pc <- sapply(out.conf,part.conf)
  vr <- sapply(out.conf,var.ratio)
  ncl <- sapply(out.conf,function(X) length(unique(X$clustering)))
  
  best.mass <- mass.algorithm(mass,pc,vr,ncl,w)
  
  list(selection=best.mass, list.epam=list.epam, conf= out.conf[[which(mass==best.mass[2,"mass"])]])
}
