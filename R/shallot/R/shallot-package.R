#' Random Partition Distribution Indexed by Pairwise Information
#' 
#' This package implements models described in the paper
#' \href{https://doi.org/10.1080/01621459.2016.1165103}{Dahl, D. B., Day, R., and
#' Tsai, J. (2017), Random Partition Distribution Indexed by Pairwise
#' Information, \emph{Journal of the American Statistical Association},
#' 112, 721-732.} The Ewens, Ewens-Pitman, Ewens attraction, Ewens-Pitman
#' attraction, and ddCRP distributions are available for prior simulation.  We
#' hope in the future to add posterior simulation with a user-supplied
#' likelihood.  Supporting functions for partition estimation and plotting are
#' also planned.
#' 
#' 
#' @name shallot-package
#' @aliases shallot-package shallot
#' @docType package
#' @author David B. Dahl \email{dahl@@stat.byu.edu}
#' @seealso \code{\link{ewens.pitman.attraction}},
#' \code{\link{sample.partitions}}
#' @references \href{https://doi.org/10.1080/01621459.2016.1165103}{Dahl, D. B.,
#' Day, R., and Tsai, J. (2017), Random Partition Distribution Indexed by
#' Pairwise Information, \emph{Journal of the American Statistical
#' Association}, 112, 721-732. <DOI:10.1080/01621459.2016.1165103>}
#' @keywords package
#' @examples
#' 
#' \donttest{
#' data <- iris[,-ncol(iris)]
#' truth <- as.integer(iris[,ncol(iris)])
#' distance <- as.dist(as.matrix(dist(scale(data))+0.001))
#' 
#' decay <- decay.exponential(temperature(9.0, fixed=TRUE), distance)
#' permutation <- permutation(n.items=nrow(data), fixed = FALSE)
#' attraction <- attraction(permutation, decay)
#' mass <- mass(1.0, fixed = TRUE)
#' discount <- discount(0.2, fixed = TRUE)
#' distribution <- ewens.pitman.attraction(mass, discount, attraction)
#' 
#' raw <- sample.partitions(distribution, 500, parallel=FALSE)
#' samples <- process.samples(raw)
#' 
#' \dontshow{
#' rscala::scalaDisconnect(shallot:::s)
#' }
#' }
#' 
NULL



