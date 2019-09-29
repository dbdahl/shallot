#' @import rscala
scalaPush.clustering <- function(x, bridge, withParameters=TRUE) {
  singleton <- is.vector(x)
  if ( singleton ) x <- matrix(x,nrow=1)
  if ( ! is.matrix(x) ) stop("Object is not a matrix.")
  if ( nrow(x) == 0 ) stop("Object should not be empty.")
  result <- s(x=x) ^ 'x.map(Clustering.apply).toList'
  if ( singleton ) result$head() else result
}

#' @import rscala
scalaPull.clustering <- function(reference, bridge, names=NULL, withParameters=TRUE) {
  singleton <- FALSE 
  if ( reference$"isInstanceOf[org.ddahl.sdols.clustering.Clustering[_]]"() ) {
    reference <- s(x=reference) ^ 'List(x)'
    singleton <- TRUE
  }
  type <- scalaType(reference)
  if ( ! reference$"isInstanceOf[List[org.ddahl.sdols.clustering.Clustering[_]]]"() ) stop(paste0("Unexpected reference type: ",type))
  if ( withParameters && ( substr(type,nchar(type)-26,nchar(type)) == "[org.ddahl.rscala.RObject]]" ) ) {
    r <- s(ref=reference) ^ '
      val list = ref.map { partition =>
        (partition.toLabelsWithParameters, partition.nClusters)
      }
      val labels = list.map(_._1._1).toArray
      val parameters = list.map(_._1._2).flatten
      val sizes = list.map(_._2).toArray
      (labels, parameters, sizes)
    '
    labels <- r$"_1"()
    parameters <- scalaPull(r$"_2"(),"generic",bridge)
    sizes <- r$"_3"()
    p <- vector(length(sizes), mode="list")
    j <- 1
    for ( i in seq_along(sizes) ) {
      p[[i]] <- parameters[j:(j+sizes[i]-1)]
      j <- j + sizes[i]
    }
    parameters <- p
  } else {
    labels <- s(ref=reference) * 'ref.map(_.toLabels).toArray'
    parameters <- NULL
  }
  labels <- labels + 1L
  colnames(labels) <- names
  if ( singleton ) {
    r <- labels[1,]
    if ( ! is.null(parameters) ) attr(r,"parameters") <- parameters[1,]
    r
  } else list(labels=labels, parameters=parameters)
}

