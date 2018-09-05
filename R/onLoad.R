.onLoad <- function(libname, pkgname) {
  s <- sdols:::s
  scalaJARs(pkgname)
  scalaLazy(function(s) s + '
    import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
    import org.ddahl.shallot._
    import parameter._
    import parameter.decay._
    import parameter.partition._
    import distribution._
    import mcmc._

    def rdg() = {
      val ints = R.evalI1("as.integer(runif(2,-.Machine$integer.max,.Machine$integer.max))")
      val seed = ((ints(0).asInstanceOf[Long]) << 32) | (ints(1) & 0xffffffffL)
      val r = new RDG()
      r.reSeed(seed)
      r
    }

    // This circumvents a bug in the class loader on some versions of Scala/JVM.
    scala.util.Try {
      new org.apache.commons.math3.random.EmpiricalDistribution()
    }
  ')
  assign("s",s,envir=parent.env(environment()))  
}
