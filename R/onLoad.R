.onLoad <- function(libname, pkgname) {
  snippet <- '
    import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
    import org.ddahl.shallot._
    import parameter._
    import parameter.decay._
    import parameter.partition._
    import distribution._
    import mcmc._

    def rdg() = {
      val ints = R.evalI1("runif(2,-.Machine$integer.max,.Machine$integer.max)")
      val seed = ((ints(0).asInstanceOf[Long]) << 32) | (ints(1) & 0xffffffffL)
      val r = new RDG()
      r.reSeed(seed)
      r
    }
  '
  ## Users may want to use 'options(rscala.heap.maximum="2G")'.
  .rscalaPackage(pkgname,classpath.packages="commonsMath",snippet=snippet)
  ## This circumvents a weird bug in the class loader of Scala 2.11.x.
  sInfo <- scalaInfo()
  if ( ( ! is.null(sInfo) ) && ( sInfo$major.version == "2.11" ) ) {
    s$.org.apache.commons.math3.random.EmpiricalDistribution$new()
  }
}

.onUnload <- function(libpath) {
  .rscalaPackageUnload()
}

