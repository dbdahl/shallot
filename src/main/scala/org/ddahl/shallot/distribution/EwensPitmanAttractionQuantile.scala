package org.ddahl.shallot
package distribution

import org.apache.commons.math3.special.Gamma.logGamma
import org.apache.commons.math3.util.FastMath.log
import parameter._
import parameter.partition._

class EwensPitmanAttractionQuantile[A] private (val samplingModel: SamplingModel[A], val mass: Mass, val discount: Discount, val attraction: Attraction, val p: Double, val timesSize: Boolean) extends PartitionModel[A] with HasSamplingModel[EwensPitmanAttractionQuantile[A], A] with HasMass[EwensPitmanAttractionQuantile[A], A] with HasDiscount[EwensPitmanAttractionQuantile[A], A] with HasAttraction[EwensPitmanAttractionQuantile[A], A] {

  def self = this

  def replaceSamplingModel(newSamplingModel: SamplingModel[A]) = new EwensPitmanAttractionQuantile(newSamplingModel, mass, discount, attraction, p, timesSize)

  def replaceMass(newMass: Mass) = new EwensPitmanAttractionQuantile(samplingModel, newMass, discount, attraction, p, timesSize)

  def replaceDiscount(newDiscount: Discount) = new EwensPitmanAttractionQuantile(samplingModel, mass, newDiscount, attraction, p, timesSize)

  def replaceAttraction(newAttraction: Attraction) = new EwensPitmanAttractionQuantile(samplingModel, mass, discount, newAttraction, p, timesSize)

  override def iteratorForForwardSampling(nItems: Int) = {
    if (nItems != attraction.nItems) throw new IllegalArgumentException("Number of items is inconsistent.")
    attraction.permutation.iterator
  }

  override def logFullConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    val p2 = partition.add(i, subset)
    var p = Partition.empty[A]()
    var sum = 0.0
    val perm = attraction.permutation
    iteratorForForwardSampling(p2.nItems).foreach(k => {
      val sOption = p.find(s => p2.paired(k, s.head))
      val s = if (sOption.isEmpty) Subset.empty(samplingModel.sample()) else sOption.get
      if (perm.nPreceeding(k) >= perm.nPreceeding(i)) sum += logCascadingConditional(k, s, p)
      p = p.add(k, s)
    })
    sum
  }

  private def quantile(x: Iterable[Double]): Double = {
    if (p >= 1.0) return x.max
    if (p <= 0.0) return x.min
    val y = x.toArray.sortWith(_ < _)
    if (y.length == 1) return y(0)
    val h = (y.length - 1) * p
    val hFloor = h.floor.toInt
    y(hFloor) + (h - hFloor) * (y(hFloor + 1) - y(hFloor)) // Type 7 quantile in R
  }

  private def func(x: Iterable[Double]) = {
    if (p == 2.0) {
      x.sum / (if (timesSize) 1.0 else x.size)
    } else {
      quantile(x) * (if (timesSize) x.size else 1.0)
    }
  }

  def logCascadingConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    val numerator = if (subset.size == 0) {
      mass.value + partition.nClusters * discount.value
    } else {
      (partition.nItems - partition.nClusters * discount.value) * func(subset.map(j => attraction.similarity(i, j))) / partition.map(s => func(s.map(j => attraction.similarity(i, j)))).sum
    }
    log(numerator / (mass.value + partition.nItems))
  }

}

object EwensPitmanAttractionQuantile extends TwoParameterNumberOfSubset {

  def apply[A](samplingModel: SamplingModel[A], mass: Mass, discount: Discount, attraction: Attraction, p: Double, timesSize: Boolean) = new EwensPitmanAttractionQuantile(samplingModel, mass, discount, attraction, p, timesSize)

  def factory[A](samplingModel: SamplingModel[A], massFactory: () => Mass, discountFactory: () => Discount, attractionFactory: () => Attraction, p: Double, timesSize: Boolean) = () => apply(samplingModel, massFactory(), discountFactory(), attractionFactory(), p, timesSize)

}

