package org.ddahl.shallot
package distribution

import org.apache.commons.math3.special.Gamma.logGamma
import org.apache.commons.math3.util.FastMath.{log, exp, sqrt}
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import org.apache.commons.math3.distribution.NormalDistribution
import parameter.Similarity

trait PPMxSimilarity {
  
  def logValue(i: Int, subset: Array[Int]): Double = {
    logValue(subset :+ i) - logValue(subset)
  }

  def logValue(subset: Array[Int]): Double
  
}

/*
import org.ddahl.rscala.RClient
class PPMxSimilarityInR(R: RClient) extends PPMxSimilarity {
  
  def logValue(subset: Array[Int]): Double = {
    R.evalD0(s"similarityFunction(c("+subset.map(_+1).mkString(",")+"))")
  }

}
*/

class PPMxSimilarityCategorical(val covariates: Array[Array[String]], val levels: Array[String], val alpha: Array[Double]) extends PPMxSimilarity {

  if ( levels.length != alpha.length ) throw new IllegalArgumentException("Lengths of levels and alpha must be the same.")
  locally {
    val n = covariates(0).length
    for ( i <- 1 until covariates.length ) {
      if ( covariates(i).length != n ) throw new IllegalArgumentException("Number of covariates is not consistent.")
    }
  }

  private val map = levels.zipWithIndex.toMap
  private val sumAlpha = alpha.sum
  private val const = logGamma(sumAlpha) - alpha.map(logGamma).sum
  private val k = covariates.length
  
  override def logValue(i: Int, subset: Array[Int]): Double = {
    covariates.map( cov => {
      val t = cov(i)
      val nt = subset.count(cov(_) == t)
      log(alpha(map(t)) + nt)
    }).sum - k*log(subset.length+sumAlpha)
  }
  
  def logValue(subset: Array[Int]): Double = {
    covariates.map( cov => {
      val beta = new Array[Int](alpha.length)
      subset.foreach( i => beta(map(cov(i))) += 1 )
      val alphaBeta = Range(0,alpha.length).map( i => alpha(i) + beta(i) )
      const - logGamma(subset.length+sumAlpha) + alphaBeta.map(logGamma).sum
    }).sum
  }

}

class PPMxSimilarityContinuousMeanVariance(val covariates: Array[Array[Double]], val mu0: Double, val n0: Double, val alpha: Double, val beta: Double) extends PPMxSimilarity {
  
  // Student t distribution with Bernardo & Smith (2000) parameterization
  private def evalLogDensity(y: Double, df: Double, location: Double, invScale: Double) = {
    logGamma((df+1)/2.0) - logGamma(df/2.0) + 0.5*log(invScale/(df*math.Pi)) - (df+1)/2.0 * log(1+invScale*(y-location)*(y-location)/df)
  }
  
  override def logValue(i: Int, subset: Array[Int]): Double = {
    val n = subset.length.toDouble
    covariates.map( cov => {
      if ( n == 0 ) {
        val df = 2*alpha
        val invScale = n0*alpha/((n0+1)*beta)
        evalLogDensity(cov(i),df,mu0,invScale)
      } else {
        val stats = new SummaryStatistics()
        subset.foreach(j => stats.addValue(cov(j)))
        val xbar = stats.getMean
        val s2 = stats.getVariance
        val df = 2*alpha+n
        val location = (n0*mu0+n*xbar)/(n0+n)
        val betan = beta + n*s2/2 + (n0*n*(mu0-xbar)*(mu0-xbar))/(2*(n0+n))
        val invScale = (n+n0)*(alpha+n/2)/((n+n0+1)*betan)
        evalLogDensity(cov(i),df,location,invScale)
      }
    }).sum
  }
  
  def logValue(subset: Array[Int]): Double = {
    var subset = List[Int]()
    subset.map( i => {
      val result = logValue(i,subset.toArray)
      subset = i +: subset
      result
    }).sum
  }

}

class PPMxSimilarityContinuousMean(val covariates: Array[Array[Double]]) extends PPMxSimilarity {
  
  val mu0: Double = 0      // Because covariates are standardized.
  val S: Double   = 1      // Because covariates are standardized.
  val c1: Double  = 0.5
  val c2: Double  = 10
  
  val v = c1 * S     // variance
  val B = c2 * S     // variance
  val lambda = 1/v   // precision
  val lambda0 = 1/B  // precision
  
  override def logValue(i: Int, subset: Array[Int]): Double = {
    val n = subset.length.toDouble
    covariates.map( cov => {
      val xbar = if ( n != 0 ) {
        val stats = new SummaryStatistics()
        subset.foreach(j => stats.addValue(cov(j)))
        stats.getMean
      } else 0.0
      val lambdan = lambda0 + n*lambda
      val mun = (lambda0*mu0 + n*lambda*xbar)/lambdan
      val sdn = sqrt(1/lambda + 1/lambdan)
      val normal = new NormalDistribution(null, mun, sdn)
      normal.logDensity(cov(i))
    }).sum
  }
  
  def logValue(subset: Array[Int]): Double = {
    var subset = List[Int]()
    subset.map( i => {
      val result = logValue(i,subset.toArray)
      subset = i +: subset
      result
    }).sum
  }

}

class PPMxSimilarityComposite(val ppmxSimilarities: Seq[PPMxSimilarity]) extends PPMxSimilarity {
  
  override def logValue(i: Int, subset: Array[Int]): Double = ppmxSimilarities.map(_.logValue(i,subset)).sum

  def logValue(subset: Array[Int]): Double = ppmxSimilarities.map(_.logValue(subset)).sum
  
}
