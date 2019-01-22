package org.ddahl.shallot
package distribution

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.exp
import org.apache.commons.math3.util.FastMath.log
import parameter._
import parameter.partition._

trait PartitionModel[A] {

  // Items that must be defined

  val samplingModel: SamplingModel[A]

  def logCascadingConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double

  // Items whose implementation may need to be overridden

  val partitionProbabilitiesSumToOne: Boolean = true // Override if necessary

  def iteratorForForwardSampling(nItems: Int) = Range(0, nItems).iterator // Override if not exchangeable

  def logFullConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = logCascadingConditional(i, subset, partition) // Override if not exchangeable

  def logProbability(partition: Partition[A]): Double = { // Override if there is something more efficient
    var p = Partition.empty[A]()
    var sum = 0.0
    iteratorForForwardSampling(partition.nItems).foreach(i => {
      val sOption = p.find(s => partition.paired(i, s.head))
      val s = if (sOption.isEmpty) Subset.empty(samplingModel.sample()) else sOption.get
      sum += logCascadingConditional(i, s, p)
      p = p.add(i, s)
    })
    sum
  }

  // Final methods

  final def probability(partition: Partition[A]): Double = exp(logProbability(partition: Partition[A]))

  final def fullConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = exp(logFullConditional(i, subset, partition))

  final def probabilitySequential(i: Int, subset: Subset[A], partition: Partition[A]): Double = exp(logCascadingConditional(i, subset, partition))

  final def probabilitySequential(i: Int, partition: Partition[A]): Distribution[A] = {
    val nItems = partition.nItems
    val elements = partition.map(subset => (subset, logCascadingConditional(i, subset, partition)))
    val emptySubset = Subset.empty(samplingModel.sample())
    Distribution(partition, i, Some((emptySubset, logCascadingConditional(i, emptySubset, partition))) ++ elements, true, true)
  }

  final def fullConditional(i: Int, partition: Partition[A]): Distribution[A] = {
    val cache = scala.collection.mutable.HashMap[Subset[A], Double]()
    val nItems = partition.nItems
    val elements = partition.map(subset => (subset, logFullConditional(i, subset, partition)))
    val emptySubset = Subset.empty(samplingModel.sample())
    Distribution(partition, i, Some((emptySubset, logFullConditional(i, emptySubset, partition))) ++ elements, true, true)
  }

  final def forwardSample(nItems: Int, rdg: RDG) = {
    iteratorForForwardSampling(nItems).foldLeft(Partition.empty[A]())((partition, i) => probabilitySequential(i, partition).sample(rdg))
  }

  final def forwardSampler(nItems: Int, rdg: RDG) = () => forwardSample(nItems, rdg)

  final def tabulate(nItems: Int) = {
    Distribution(Partition.enumerate(() => samplingModel.sample(), nItems).map(partition => (partition, logProbability(partition))), true, !partitionProbabilitiesSumToOne)
  }

}

object PartitionModel {

  def forwardSampler[A](nItems: Int, model: PartitionModel[A]): (Int, RDG) => List[Partition[A]] = {
    (nSamples: Int, rdg: RDG) => List.fill(nSamples) { model.forwardSample(nItems, rdg) }
  }

  def forwardSampler[A](nItems: Int, factory: () => PartitionModel[A]): (Int, RDG) => List[Partition[A]] = {
    (nSamples: Int, rdg: RDG) => List.fill(nSamples) { factory().forwardSample(nItems, rdg) }
  }

  def tabulator[A](nItems: Int, factory: () => PartitionModel[A]): (Int) => List[Distribution[A]] = {
    (nSamples: Int) => List.fill(nSamples) { factory().tabulate(nItems) }
  }

}

