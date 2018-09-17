package org.ddahl.shallot
package distribution

import org.apache.commons.math3.special.Gamma.logGamma
import org.apache.commons.math3.util.FastMath.log
import parameter._
import parameter.partition._

class EwensPitman[A] private (val samplingModel: SamplingModel[A], val mass: Mass, val discount: Discount) extends PartitionModel[A] with HasSamplingModel[EwensPitman[A], A] with HasMass[EwensPitman[A], A] with HasDiscount[EwensPitman[A], A] {

  def self = this

  def replaceSamplingModel(newSamplingModel: SamplingModel[A]) = new EwensPitman(newSamplingModel, mass, discount)

  def replaceMass(newMass: Mass) = new EwensPitman(samplingModel, newMass, discount)

  def replaceDiscount(newDiscount: Discount) = new EwensPitman(samplingModel, mass, newDiscount)

  private val logGammaOneMinusDiscount = logGamma(1.0 - discount.value)

  override def logProbability(partition: Partition[A]) = {
    val v1 = Range(0, partition.nClusters).foldLeft(0.0)((sum, i) => sum + log(mass.value + i * discount.value)) - logGamma(mass.value + partition.nItems) + mass.logGammaValue
    val v2 = partition.foldLeft(0.0)((sum, subset) => { sum + logGamma(subset.size - discount.value) - logGammaOneMinusDiscount })
    v1 + v2
  }

  def logCascadingConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    val numerator = if (subset.size == 0) mass.value + partition.nClusters * discount.value
    else subset.size - discount.value
    log(numerator / (mass.value + partition.nItems))
  }

  def logPredictive(partition: Partition[A]): Iterable[(Subset[A],Double)] = {
    val p = partition ++ Iterable(Subset.empty(samplingModel.sample()))
    p.map(subset => {
      val numerator = if (subset.size == 0) mass.value + partition.nClusters * discount.value
      else subset.size - discount.value
      (subset, log(numerator / (mass.value + partition.nItems)))
    })
  }
  
}

object EwensPitman extends TwoParameterNumberOfSubset {

  def apply[A](samplingModel: SamplingModel[A], mass: Mass, discount: Discount) = new EwensPitman(samplingModel, mass, discount)

  def factory[A](samplingModel: SamplingModel[A], massFactory: () => Mass, discountFactory: () => Discount) = () => apply(samplingModel, massFactory(), discountFactory())

}

