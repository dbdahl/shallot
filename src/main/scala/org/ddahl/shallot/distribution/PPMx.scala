package org.ddahl.shallot
package distribution

import org.apache.commons.math3.special.Gamma.logGamma
import org.apache.commons.math3.util.FastMath.log
import parameter._
import parameter.partition._


class PPMx[A] private (val samplingModel: SamplingModel[A], val mass: Mass, val ppmxSimilarity: PPMxSimilarity) extends PartitionModel[A] with HasSamplingModel[PPMx[A], A] with HasMassEscobarWest[PPMx[A], A] {
  
  def self = this

  def replaceSamplingModel(newSamplingModel: SamplingModel[A]) = new PPMx(newSamplingModel, mass, ppmxSimilarity)  
  
  def replaceMass(newMass: Mass) = new PPMx(samplingModel, newMass, ppmxSimilarity)

  override def logProbability(partition: Partition[A]) = {
    val v1 = mass.logValue * partition.nClusters + mass.logGammaValue - logGamma(mass.value + partition.nItems)
    val v2 = partition.foldLeft(0.0)((sum, subset) => sum + ppmxSimilarity.logValue(subset.toArray) + logGamma(subset.size))
    v1 + v2
  }
  
  override val partitionProbabilitiesSumToOne: Boolean = false
  
  def logCascadingConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    val numerator = if (subset.size == 0) mass.value
    else subset.size
    ppmxSimilarity.logValue(i, subset.toArray) + log(numerator / (mass.value + partition.nItems))
  }

}

object PPMx {

  def apply[A](samplingModel: SamplingModel[A], mass: Mass, ppmxSimilarity: PPMxSimilarity) = new PPMx(samplingModel, mass, ppmxSimilarity)

  def factory[A](samplingModel: SamplingModel[A], massFactory: () => Mass, ppmxSimilarity: PPMxSimilarity) = () => apply(samplingModel, massFactory(), ppmxSimilarity)

}
