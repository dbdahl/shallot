.onLoad <- function(libname, pkgname) {
  snippet <- '
    import org.ddahl.shallot._
    import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }

    def rdg() = {
      val r = new org.apache.commons.math3.random.RandomDataGenerator()
      val ints = R.evalI1("runif(2,-.Machine$integer.max,.Machine$integer.max)")
      val seed = ((ints(0).asInstanceOf[Long]) << 32) | (ints(1) & 0xffffffffL)
      r.reSeed(seed)
      r
    }
  '
  ## Tell users to set heap.maximum using options("rscala.command.line.options"="-J-Xmx4g")
  .rscalaPackage(pkgname,classpath.packages="commonsMath",snippet=snippet)
  ## This circumvents a weird bug in the class loader of Scala 2.11.x.
  sInfo <- scalaInfo()
  if ( ( ! is.null(sInfo) ) && ( sInfo$major.version == "2.11" ) ) {
    s$.org.apache.commons.math3.random.EmpiricalDistribution$new()
  }
#  s[['env']]$mapper.nSubsets     <- s %.~% '(p: org.ddahl.shallot.parameter.partition.Partition[org.ddahl.shallot.REXP]) => p.nSubsets' 
#  s[['env']]$mapper.entropy      <- s %.~% '(p: org.ddahl.shallot.parameter.partition.Partition[org.ddahl.shallot.REXP]) => p.entropy'
#  s[['env']]$mapper.mass         <- s %.~% '(p: org.ddahl.shallot.AugmentedSample) => p._2'
#  s[['env']]$mapper.discount     <- s %.~% '(p: org.ddahl.shallot.AugmentedSample) => p._3'
#  s[['env']]$mapper.temperature  <- s %.~% '(p: org.ddahl.shallot.AugmentedSample) => p._5'
#  s[['env']]$mapper.toLabelsWithParameters <- s %.~% '(p: org.ddahl.shallot.parameter.partition.Partition[org.ddahl.shallot.REXP]) => org.ddahl.shallot.r.RInterface.toLabelsWithParameters(p)'
#  s[['env']]$mapper.toLabels     <- s %.~% '(x: (Array[Int], Array[org.ddahl.shallot.REXP])) => x._1'
#  s[['env']]$mapper.toParameters <- s %.~% '(x: (Array[Int], Array[org.ddahl.shallot.REXP])) => x._2'
  invisible()
}

.onUnload <- function(libpath) {
  .rscalaPackageUnload()
}

