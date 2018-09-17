package org.ddahl.shallot
package distribution

import parameter.Mass
import parameter.Discount

trait OneParameterNumberOfSubset {

  def sampleNumberOfSubsets(nItems: Int, massFactory: () => Mass, nSamples: Int): Array[Int] = {
    val samples = new Array[Int](nSamples)
    for (r <- 0 until nSamples) {
      var k = 1
      val mass = massFactory()
      for (t <- 1 until nItems) {
        val p = mass.value / (mass.value + t)
        if (math.random <= p) k += 1
      }
      samples(r) = k
    }
    samples
  }

  def meanNumberOfSubsets(nItems: Int, mass: Mass) = {
    val w = new Array[Double](nItems)
    w(0) = 1.0 / mass.value
    for (t <- 1 until nItems) {
      w(t) = 1.0 / (mass.value + t)
    }
    mass.value * w.sum
  }

  def varianceNumberOfSubsets(nItems: Int, mass: Mass) = {
    val w = new Array[Double](nItems)
    w(0) = 0.0
    for (t <- 1 until nItems) {
      w(t) = t / ((mass.value + t) * (mass.value + t))
    }
    mass.value * w.sum
  }

  def probabilityNumberOfSubsets(nItems: Int, nSubsets: Int, mass: Mass): Double = {
    def snfk(n: Int, k: Int): Double = {
      if (n == 0 || k == 0) if (n == k) 1.0 else 0.0
      else (n - 1) * snfk(n - 1, k) + snfk(n - 1, k - 1)
    }
    if (1 <= nSubsets && nSubsets <= nItems && nItems > 0) math.pow(mass.value, nSubsets) * snfk(nItems, nSubsets) / (0 until nItems).map(_ + mass.value).product
    else 0.0
  }

}
