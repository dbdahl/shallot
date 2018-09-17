package org.ddahl.shallot
package distribution

import parameter.Mass
import parameter.Discount

trait TwoParameterNumberOfSubset {

  def sampleNumberOfSubsets(nItems: Int, massFactory: () => Mass, discountFactory: () => Discount, nSamples: Int): Array[Int] = {
    val samples = new Array[Int](nSamples)
    for (r <- 0 until nSamples) {
      var k = 1
      val mass = massFactory()
      val discount = discountFactory()
      for (t <- 1 until nItems) {
        val p = (mass.value + discount.value * k) / (mass.value + t)
        if (math.random <= p) k += 1
      }
      samples(r) = k
    }
    samples
  }

  def meanNumberOfSubsets(nItems: Int, mass: Mass, discount: Discount) = {
    val w = new Array[Double](nItems)
    w(0) = 1.0
    for (t <- 1 until nItems) {
      w(t) = (mass.value + discount.value * w.view(0, t).sum) / (mass.value + t)
    }
    w.sum
  }

  def probabilityNumberOfSubsets(nItems: Int, nSubsets: Int, mass: Mass, discount: Discount): Double = {
    val m = mass.value
    val d = discount.value
    def numerator(nNewCustomers: Int, nNewTables: Int, nSeatedCustomers: Int, nPreviousTables: Int): Double = {
      if (nNewTables == 0) (0 until nNewCustomers).map(x => nSeatedCustomers + x - nPreviousTables * d).product
      else if (nNewCustomers == nNewTables) (0 until nNewCustomers).map(x => (x + nPreviousTables) * d + m).product
      else (m + nPreviousTables * d) * numerator(nNewCustomers - 1, nNewTables - 1, nSeatedCustomers + 1, nPreviousTables + 1) + (nSeatedCustomers - nPreviousTables * d) * numerator(nNewCustomers - 1, nNewTables, nSeatedCustomers + 1, nPreviousTables)
    }
    m * numerator(nItems - 1, nSubsets - 1, 1, 1) / (0 until nItems).map(_ + m).product
  }

}
