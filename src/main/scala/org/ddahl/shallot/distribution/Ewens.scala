package org.ddahl.shallot
package distribution

import org.apache.commons.math3.special.Gamma.logGamma
import org.apache.commons.math3.util.FastMath.log
import parameter._
import parameter.partition._

class Ewens[A] private (val samplingModel: SamplingModel[A], val mass: Mass) extends PartitionModel[A] with HasSamplingModel[Ewens[A], A] with HasMassEscobarWest[Ewens[A], A] {

  def self = this

  def replaceSamplingModel(newSamplingModel: SamplingModel[A]) = new Ewens(newSamplingModel, mass)

  def replaceMass(newMass: Mass) = new Ewens(samplingModel, newMass)

  override def logProbability(partition: Partition[A]) = {
    val v1 = mass.logValue * partition.nClusters + mass.logGammaValue - logGamma(mass.value + partition.nItems)
    val v2 = partition.foldLeft(0.0)((sum, subset) => sum + logGamma(subset.size))
    v1 + v2
  }

  def logCascadingConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    val numerator = if (subset.size == 0) mass.value
    else subset.size
    log(numerator / (mass.value + partition.nItems))
  }

  def logPredictive(partition: Partition[A]): Iterable[(Subset[A],Double)] = {
    val p = partition ++ Iterable(Subset.empty(samplingModel.sample()))
    p.map(subset => {
      val numerator = if (subset.size == 0) mass.value
      else subset.size
      (subset, log(numerator / (mass.value + partition.nItems)))
    })
  }

}

object Ewens extends OneParameterNumberOfSubset {

  def apply[A](samplingModel: SamplingModel[A], mass: Mass) = new Ewens(samplingModel, mass)

  def factory[A](samplingModel: SamplingModel[A], massFactory: () => Mass) = () => apply(samplingModel, massFactory())

}

