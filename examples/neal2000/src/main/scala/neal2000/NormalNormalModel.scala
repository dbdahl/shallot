package neal2000

import org.ddahl.shallot.parameter.SamplingModel
import org.ddahl.shallot.parameter.partition.Subset
import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.{ sqrt, log }

class NormalNormalModel(data: Array[Double], sigma: Double, mu0: Double, sigma0: Double, rdg: RDG) extends SamplingModel[Double] {

  private val s2 = sigma * sigma
  private val s02 = sigma0 * sigma0
  private val s02Inv = 1.0 / s02
  private val c = -1.0 / (2.0 * s2)

  def sample = rdg.nextGaussian(mu0, sigma0)

  def sample(subset: Subset[Double]) = {
    val sum = subset.foldLeft(0.0)((sum, i) => sum + data(i))
    val variance = 1 / (s02Inv + subset.nItems / s2)
    val mean = variance * (mu0 / s02 + sum / s2)
    rdg.nextGaussian(mean, sqrt(variance))
  }

  def copy(x: Double) = x

  def logDensity(i: Int, subset: Subset[Double]): Double = {
    val resid = data(i) - subset.parameter
    c * resid * resid
  }

}

