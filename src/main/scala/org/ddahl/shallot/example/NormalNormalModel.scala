package org.ddahl.shallot
package example

import parameter._
import parameter.partition._
import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.sqrt
import org.apache.commons.math3.util.FastMath.log

class NormalNormalModel(data: Array[Double], sigma: Double, mu0: Double, sigma0: Double, rdg: RDG) extends SamplingModel[Double] {

  private val s2 = sigma * sigma
  private val s02 = sigma0 * sigma0
  private val s02Inv = 1.0 / s02
  private val c = -1.0 / (2.0 * s2)

  def sample = rdg.nextGaussian(mu0, sigma0)

  def sample(subset: Subset[Double]) = {
    val sum = subset.foldLeft(0.0)((sum, i) => sum + data(i))
    val variance = 1 / (s02Inv + subset.size / s2)
    val mean = variance * (mu0 / s02 + sum / s2)
    rdg.nextGaussian(mean, sqrt(variance))
  }

  def logDensity(i: Int, subset: Subset[Double]): Double = {
    val resid = data(i) - subset.parameter
    c * resid * resid
  }

}

class IntegratedNormalNormalModel(data: Array[Double], sigma: Double, mu0: Double, sigma0: Double, rdg: RDG) extends SamplingModel[Null] {

  private val lambda = 1.0 / (sigma * sigma)
  private val lambda0 = 1.0 / (sigma0 * sigma0)
  private val log2pi = log(2 * math.Pi)
  private val lambda0mu0 = lambda0 * mu0

  def sample = null

  def sample(subset: Subset[Null]) = null

  def logDensity(i: Int, subset: Subset[Null]): Double = {
    val sum = subset.foldLeft(0.0)((sum, i) => sum + data(i))
    val precision = lambda0 + subset.size * lambda
    val mean = (lambda0mu0 + lambda * sum) / precision
    val precision2 = lambda * precision / (precision + lambda)
    val resid = data(i) - mean
    0.5 * (log(precision2) - log2pi - precision2 * resid * resid)
  }

}

