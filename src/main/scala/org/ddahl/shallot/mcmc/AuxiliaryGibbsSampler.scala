package org.ddahl.shallot
package mcmc

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import parameter._
import parameter.partition._
import distribution._

object AuxiliaryGibbsSampler {

  def apply[A](partition: Partition[A], samplingModel: SamplingModel[A], priorModel: PartitionModel[A], rdg: RDG): (Partition[A], Boolean) = {
    var p = partition
    val nItems = p.nItems
    def logPosterior(i: Int, subset: Subset[A]): Double = {
      priorModel.logFullConditional(i, subset, p) + samplingModel.logDensity(i, subset)
    }
    for (i <- 0 until nItems) {
      val (tmpPartition, formerSubsetForI) = p.removeWithCluster(i)
      p = tmpPartition
      val emptySubset = if (formerSubsetForI.size > 0) Subset.empty(samplingModel.sample())
      else formerSubsetForI
      val elements = p.map(subset => (subset, logPosterior(i, subset)))
      val elementForEmpty = (emptySubset, logPosterior(i, emptySubset))
      val distribution = Distribution(p, i, Some(elementForEmpty) ++ elements, true, true)
      p = distribution.sample(rdg)
    }
    p = p.replace(samplingModel.sample)
    (p, true)
  }

}

