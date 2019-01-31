#' Association Matrix
#' 
#' This function creates an association matrix for a clustering/partition. 
#' The (i,j) element of the matrix is 1 if item i and j are in the same 
#' cluster/subset and 0 otherwise.
#' 
#' @param cl A vector containing cluster labels for a clustering/partition.
#' @return A matrix of 0s and 1s indicating whether items i and j are in
#' the same cluster/subset.
#' @examples
#' 
#' cl <- rep(1:3,times=c(2,4,3))
#' association.matrix(cl)
#' 
#' @export association.matrix
association.matrix <- function(cl){
  n <- length(cl)
  cl1 <- rep(cl,n)
  cl2 <- rep(cl,each=n)
  matrix(cl1==cl2, ncol=n)*1
}


#' Variance Ratio
#'
#' This function calculates the variance of the expected pairwise allocation
#' matrix (EPAM) within clusters/subsets over the total variance of the expected
#' pairwise allocation matrix.
#'
#' The \code{\link{variance.ratio}} function takes as input an object of class
#' \code{sdols.confidence} and calculates the variance ratio for the estimated
#' partition from the corresponding expected pairwise allocation matrix (EPAM).
#'
#' The variance ratio is the weighted average of the within cluster variances of
#' the EPAM, weighted by the number of pairwise EPAM values per cluster, over
#' the total variance of the EPAM.
#'
#' @param x,y If \code{y} is not specified then \code{x} must be an object of
#'   class \code{sdols.confidence}. Otherwise, \code{x} is a vector of cluster
#'   labels and \code{y} is an expected pairwise allocation matrix.
#'
#' @examples
#' x <- rep(c(1,2,3), times=c(2,3,5))
#' y <- diag(10)
#' y[upper.tri(y)] <- runif(45)
#' variance.ratio(x,y)
#' 
#' @family Default Mass Selection
#' @export
variance.ratio <- function(x,y) {
  if (inherits(x,"sdols.confidence")) {
    cl.am <- association.matrix(x$clustering)
    y <- x$expectedPairwiseAllocationMatrix
    x <- x$clustering
  } else {
    cl.am <- association.matrix(x)
  }
  cl.vars <- sapply(unique(x), function(i) {
    if (sum(x == i) < 3) return(c(0,1))
    cl.am[x != i, x != i] <- 0
    z <- y[upper.tri(y) & cl.am]
    c(sum((mean(z) - z)^2)/(length(z)-1),length(z))
  })

  cl.var <- sum(apply(cl.vars,2,prod)/sum(cl.vars[2,]))
  z <- y[upper.tri(y)]
  tot.var <- sum((mean(z) - z)^2)/(length(z) - 1)

  cl.var/tot.var
}

#' Partition Confidence
#'
#' This function calculates the partition confidence of a partition estimate
#' from the corresponding expected pairwise allocation matrix (EPAM).
#'
#' The \code{\link{partition.confidence}} takes as input an object of class
#' \code{sdols.confidence} and then calculates the partition confidence from the
#' expected pairwise allocation matrix.
#'
#' The partition confidence is the average values of the EPAM for items that are
#' clustered together. Items which are in their own subset do not contribute to
#' partition confidence.
#'
#' @inheritParams variance.ratio
#'
#' @examples
#' x <- rep(c(1,2,3), times=c(2,3,5))
#' y <- diag(10)
#' y[upper.tri(y)] <- runif(45)
#' partition.confidence(x,y)
#'
#' @family Default Mass Selection
#' @export
partition.confidence <- function(x,y){
  if (inherits(x,"sdols.confidence")) {
    if (length(unique(x$clustering)) == length(x$clustering)) return(0)
    cl.am <- association.matrix(x$clustering)
    y <- x$expectedPairwiseAllocationMatrix
  } else {
    cl.am <- association.matrix(x)
  }
  mean(y[upper.tri(y) & cl.am])
}


#' Mass Selection Algorithm
#'
#' This function selects the optimal mass value for Cluster Analysis via Random
#' Partition distributions using the Ewens-Pitman attraction distribution.
#'
#' The \code{\link{mass.algorithm}} function is used internally in the
#' \code{\link{default.mass}} function.
#'
#' The default value for \code{w} is \code{c(1,1,1)}.
#'
#' The general algorithm is as follows: \enumerate{ \item Rank the partition
#' confidence (\code{pc}) and variance ratio (\code{vr}). Select the
#' \code{mass_i} value which minimizes the weigthed sum of \eqn{w_1 pc_i + w_2
#' vr_i + w_3 n_i}. }
#'
#' The two stage algorithm proceeds as follows: \enumerate{ \item Rank the
#' partition confidence (\code{pc}) and variance ratio (\code{vr}). For each
#' number of clusters \code{n} select the index which minimizes the weigthed sum
#' of \eqn{w_1 pc_i + w_2 vr_i}. \item Rerank the \code{pc} and \code{vr} of the
#' selected indices and select the \code{mass_i} value which minimizes the
#' weigthed sum of \eqn{w_1 pc_i + w_2 vr_i + w_3 n_i} from among the selected
#' indices. }
#'
#'
#' @param mass a vector of mass values
#' @param pc a vector of partition confidences for the partition estimates at
#'   the corresponding mass values
#' @param vr a vector of variance ratios for the partition estimates at the
#'   corresponding mass values
#' @param n a vector of the number of subsets in the partition estimates at the
#'   correpsonding mass values
#' @param w a vector of length 3 specifying the weights of \code{pc}, \code{vr},
#'   and \code{n}
#' @param two.stage logical; if \code{TRUE}, the two stage algorithm is
#'   implemented
#'
#' @return A matrix containing the `best' \code{mass} value and corresponding
#'   values for \code{pc}, \code{vr}, and \code{n}. The matrix also contains the
#'   mass values for the partitions estimate with more one more and one less
#'   subset that the selected mass value.
#' @examples
#'
#' @family Default Mass Selection
#' @importFrom graphics abline lines plot
#' @export
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
    out <- t(all[,c(best=index, index-1, index+1)])
  } else {
    best <- compute.rankings(n>1)
    index <- which(unique(n[n>1]) == best["n"])
    others <- sapply(unique(n[n>1])[index + c(-1,1)], function(i) compute.rankings(n == i))
    
    out <- rbind(best,t(others))
  }

  plot(mass,pc,type="l",ylim=c(0,1),col="dodgerblue")
  lines(vr[n>1]~mass[n>1],col="forestgreen")
  abline(v=out[1,"mass"])
  for (i in 2:nrow(out)) {
    abline(v=out[i,"mass"], lty=2)
  }

  out
}

#' Default Mass Selection
#'
#' This function selects an optimal mass value for Cluster Analysis via Random
#' Partition Distribtuions, using the Ewens-Pitman Attraction distribution.
#'
#' The function draws \code{n.draws} partitions at each specified mass value. If
#' a vector of mass values is not given, then the default of
#' \code{seq(0.1,10,0.2)} is used for loss
#' \code{lowerBoundVariationOfInformation} and \code{seq(0.1,5,0.05)} used for
#' the other loss functions.
#'
#' If a list of expected pairwise allocation matrices (EPAM) is provided,
#' additional draws at matching mass values are added to the corresponding
#' matrix. Additionally, no new draws are needed for estimation, if a list of
#' EPAMs is provided.
#'
#' A partition/clustering estimate from each EPAM is obtained using the SALSO
#' method in \code{\link[sdols]{salso}}. The estimate given minimizes the
#' specified \code{loss} function with respect to the EPAM.
#'
#' The function then uses the \code{\link{mass.algorithm}} to select the optimal
#' mass value for clustering estimation.
#'
#'
#' @param mass optional, a vector of mass values.
#' @param list.epam optional, a list of expected pairwise allocation matrices.
#'   Each matrix in the list needs the attributes "\code{mass}" and
#'   "\code{n.draws}".
#' @param dis a dissimilarity structure of class \code{dist}.
#' @param new.draws logical; if \code{TRUE} then new draws are obtained at each
#'   mass value.
#' @param w a vector of length 3 of the weights to be used in the
#'   \code{\link{mass.algorithm}}.
#' @param discount parameter of the Ewens-Pitman Attraction distribution.
#' @param temp temperature parameter of the Ewens-Pitman Attraction
#'   distribution.
#' @param loss One of "\code{squaredError}", "\code{absoluteError}",
#'   "\code{binder}", or "\code{lowerBoundVariationOfInformation}" to indicate
#'   the optimization should seek to minimize squared error loss, absolute error
#'   loss, Binder loss (Binder 1978), or the lower bound of the variation of
#'   information loss (Wade & Ghahramani 2017), respectively.
#' @param n.draws number of draws of partitions to be obtained at each mass
#'   value.
#' @param two.stage logical; if \code{TRUE}, the two stage algorithm is
#'   implemented in \code{\link{mass.algorithm}}.
#' @param parallel logical; if \code{TRUE} computations will take advantage
#'   multiple CPU cores.
#' @param x An object from the \code{\link{default.mass}} function.
#' @param ... currently ignored
#'
#' @return An object of class \code{shallot.default.mass}. This object is a list
#'   containing a matrix of `best' possible mass values to maximize partition
#'   confidence and minimize the variance ratio, the clustering estimate, the
#'   expected pairwise allocation matrix, parameters used for optimization and
#'   the EPA distribution, and the list of expected pairwise allocation matrices
#'   for each mass value.
#' @examples
#'
#' @family Default Mass Selection
#' @export
default.mass <- function(mass, list.epam, dis, new.draws = TRUE, w=c(1,1,1), discount=0, temp=10, loss="binder", n.draws=100L, two.stage=TRUE, parallel=TRUE) {
  if (missing(mass)) {
    if (loss == "lowerBoundVariationOfInformation") {
      mass <- seq(0.1,10,0.2)
    } else {
      mass <- seq(0.1,5,by=0.05)
    }
  }

  single.epam <- function(mass,discount,dis,n.draws,temp=10,parallel=parallel) {
    n <- nrow(as.matrix(dis))
    epa <- ewens.pitman.attraction(mass(mass), discount(discount), attraction(permutation(n.items=n, fixed=FALSE), decay.exponential(temperature(temp), dis)))
    draws <- sample.partitions(epa, n.draws, parallel = parallel)
    as.matrix(pairwise.probabilities(draws))
  }

  if(missing(list.epam) && new.draws == FALSE) stop("must provide list.epam or obtain new draws")
  if (missing(list.epam) || new.draws == TRUE) {
    out <- lapply(mass, single.epam,  n.draws=n.draws, discount=discount, dis=dis,  temp=temp, parallel=parallel)
    sapply(1:length(out), function(i) {
      attr(out[[i]],"mass") <<- mass[i]
      attr(out[[i]],"n.draws") <<- n.draws
      attr(out[[i]],"discount") <<- discount
      attr(out[[i]],"temperature") <<- temp
    })
    if (!missing(list.epam)) {
      mass0 <- sapply(list.epam,attr,"mass")
      combine <- which(mass0 %in% mass)
      add <- which(!(mass %in% mass0))
      if (length(combine) > 0) {
        sapply(combine, function(i) {
          n.draws0 <- attr(list.epam[[i]],"n.draws")
          which.out <- which(mass0[i] == mass)
          list.epam[[i]] <<- (list.epam[[i]]*n.draws0 + out[[which.out]]*n.draws)/(n.draws + n.draws0)
          attr(list.epam[[i]],"n.draws") <<- n.draws0 + n.draws
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

  out.cl <- lapply(list.epam, function(x) {
    attr(x, "mass") <- NULL
    attr(x, "discount") <- NULL
    attr(x, "n.draws") <- NULL
    attr(x, "temperature") <- NULL
    salso(x, loss=loss, multicore=parallel)
  })

  mass <- sapply(list.epam, attr, "mass")
  weight <- sapply(list.epam, attr, "n.draws")
  pc <- mapply(partition.confidence, out.cl, list.epam)
  vr <- mapply(variance.ratio, out.cl, list.epam)
  ncl <- sapply(out.cl,function(X) length(unique(X)))

  best.mass <- mass.algorithm(mass,pc,vr,ncl,w,two.stage)
  best.index <- which(mass == best.mass[1,"mass"])

  output <- list(criteria=best.mass, 
                 clustering = out.cl[[best.index]],
                 expectedPairwiseAllocationMatrix = list.epam[[best.index]],
                 params = list(discount = discount, temperature = temp, loss = loss),
                 list.epam = list.epam)
  class(output) <- "shallot.default.mass"
  output
}

#' @rdname default.mass
#' @export
print.shallot.default.mass <- function(x, ...) {
  cat("Clustering estimated using loss: ",x$params$loss,", with discount ",
      x$params$discount," and temperature ", x$params$temperature,"\n\n",sep="")
  cat("The selected mass value is", x$criteria[1,"mass"],"\n")
  cat("number of clusters: ",x$criteria[1,"n"],"\npartition confidence: ",x$criteria[1,"pc"],"\nvariance ratio: ", x$criteria[1,"vr"], "\n", sep="")
  
  cat("\nestimated clustering:\n", x$clustering)
}





