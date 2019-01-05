#' @import rscala
.onLoad <- function(libname, pkgname) {
  s <- sdols:::s
  scalaJARs(pkgname,s)
  scalaLazy(function(s) s + '
    import org.ddahl.shallot._
    import parameter._
    import parameter.decay._
    import parameter.partition._
    import distribution._
    import mcmc._

    // This circumvents a bug in the class loader on some versions of Scala/JVM.
    scala.util.Try {
      new org.apache.commons.math3.random.EmpiricalDistribution()
    }
  ')
  assign("s",s,envir=parent.env(environment()))  
}

