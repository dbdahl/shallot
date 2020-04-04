#' @import utils
#'
.onLoad <- function(libname, pkgname) {
  assign("s", NULL, envir=parent.env(environment()))
  globalVariables("s")
}

#' @import rscala
#' @import commonsMath
#'
scalaEnsure <- function() {
  if ( ! is.null(s) ) return()
  commonsMath:::.packageName   # So CRAN checks recognize that its being used.
  s <- scala("commonsMath")
  scalaLazy(function(s) s + '
    import org.ddahl.sdols._
    import org.ddahl.sdols.clustering._
    import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
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

    def rdg() = {
      val ints = R.evalI1("runif(2,-.Machine$integer.max,.Machine$integer.max)")
      val seed = ((ints(0).asInstanceOf[Long]) << 32) | (ints(1) & 0xffffffffL)
      val r = new RDG()
      r.reSeed(seed)
      r
    }
  ')
  scalaPushRegister(scalaPush.clustering,"clustering",s)
  scalaPullRegister(scalaPull.clustering,"clustering",s)
  env <- parent.env(environment())
  # unlockBinding("s", env)
  eval(parse(text=paste0('unlockBinding("s",env)')))
  assign("s", s, envir=env)
  lockBinding("s", env)
}

.onUnload <- function(libpath) {
  if ( ! is.null(s) ) close(s)
}
