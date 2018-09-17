package org.ddahl.shallot
package distribution

import org.apache.commons.math3.special.Gamma.logGamma
import org.apache.commons.math3.util.FastMath.log
import parameter._
import parameter.partition._

class EwensAttraction[A] private (val samplingModel: SamplingModel[A], val mass: Mass, val attraction: Attraction) extends PartitionModel[A] with HasSamplingModel[EwensAttraction[A], A] with HasMassEscobarWest[EwensAttraction[A], A] with HasAttraction[EwensAttraction[A], A] {

  def self = this

  def replaceSamplingModel(newSamplingModel: SamplingModel[A]) = new EwensAttraction(newSamplingModel, mass, attraction)

  def replaceMass(newMass: Mass) = new EwensAttraction(samplingModel, newMass, attraction)

  def replaceAttraction(newAttraction: Attraction) = new EwensAttraction(samplingModel, mass, newAttraction)

  override def iteratorForForwardSampling(nItems: Int) = {
    if (nItems != attraction.nItems) throw new IllegalArgumentException("Number of items is inconsistent.")
    attraction.permutation.iterator
  }

  override def logFullConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    val a1 = if (subset.size == 0) mass.logValue else 0.0
    val attrOfiToSubset = attraction(i, subset)
    val a2 = if (attrOfiToSubset > 0.0) log(attrOfiToSubset) else 0.0
    attraction.permutation.thin(subset, i).foldLeft(a1 + a2)((sum, k) => {
      val b1 = attraction(k, subset)
      val b2 = b1 + attraction(k, i)
      sum + (if (b2 > 0.0) log(b2) else 0.0) - (if (b1 > 0.0) log(b1) else 0.0)
    })
  }

  override def logProbability(partition: Partition[A]) = {
    val v1 = mass.logValue * partition.nClusters + mass.logGammaValue - logGamma(mass.value + partition.nItems)
    val v2 = partition.foldLeft(0.0)((sum, subset) => {
      val thinned = attraction.permutation.thin(subset)
      thinned.next // Discard the first which has no attraction to itself
      sum + thinned.foldLeft(0.0)((sum, i) => sum + log(attraction(i, subset)))
    })
    v1 + v2
  }

  def logCascadingConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    val numerator = if (subset.size == 0) mass.value
    else attraction(i, subset)
    log(numerator / (mass.value + partition.nItems))
  }

  def logPredictive(distanceToHeldOut: Array[Double], partition: Partition[A]): Iterable[(Subset[A],Double)] = {
    val s = distanceToHeldOut.map(attraction.decay(_)).toArray
    if (!s.forall(_ > 0.0)) throw new IllegalArgumentException("Not all the elements are strictly positive.")
    if (!s.forall(_ < Double.PositiveInfinity)) throw new IllegalArgumentException("Not all the elements are strictly less than infinity.")
    val p = partition ++ Iterable(Subset.empty(samplingModel.sample()))
    p.map(subset => {
      val numerator = if (subset.size == 0) mass.value
      else subset.foldLeft(0.0)((sum, k) => sum + s(k))
      (subset, log(numerator / (mass.value + partition.nItems)))
    })
  }

}

object EwensAttraction extends OneParameterNumberOfSubset {

  def apply[A](samplingModel: SamplingModel[A], mass: Mass, attraction: Attraction) = new EwensAttraction(samplingModel, mass, attraction)

  def factory[A](samplingModel: SamplingModel[A], massFactory: () => Mass, attractionFactory: () => Attraction) = () => apply(samplingModel, massFactory(), attractionFactory())

}

