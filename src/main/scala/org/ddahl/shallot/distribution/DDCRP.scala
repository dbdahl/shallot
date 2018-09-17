package org.ddahl.shallot
package distribution

import parameter.partition._
import parameter.Mass
import parameter.Distance
import parameter.Similarity
import parameter.SamplingModel
import parameter.decay.DecayFunction
import parameter.decay.ExponentialDecayFactory
import parameter.partition.Partition
import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.{ pow, log, exp }

class DDCRP[A] private (val samplingModel: SamplingModel[A], val mass: Mass, val similarity: Similarity) {

  val nItems = similarity.nItems

  private val denominator = Range(0, nItems).map(i => log(Range(0, nItems).map(j => similarity(i, j)).sum - similarity(i, i) + mass.value)).toArray

  def logProbability(edges: Array[Int]) = {
    Range(0, nItems).map(i => {
      val numerator = if (edges(i) == i) mass.logValue else log(similarity(i, edges(i)))
      numerator - denominator(i)
    }).sum
  }

  def sample(rdg: RDG): Partition[A] = {
    val edges = new Array[Int](nItems)
    for (i <- 0 until nItems) {
      val denom = exp(denominator(i))
      edges(i) = sample(Range(0, nItems).map(j => if (j == i) mass.value / denom else similarity(i, j) / denom), rdg)
    }
    edgesToPartition(edges)
  }

  def toDistribution = {
    val builder = DistributionBuilder[A]()
    enumerateEdges(nItems).foreach(e => {
      builder.process(edgesToPartition(e), exp(logProbability(e)))
    })
    builder.toDistribution
  }

  private def sample(probs: IndexedSeq[Double], rdg: RDG): Int = {
    val unif = rdg.nextUniform(0.0, 1.0, true)
    var cumsum = 0.0
    var index = -1
    probs.find(element => {
      index += 1
      cumsum += element
      cumsum > unif
    })
    index
  }

  private def edgesToPartition(edges: Array[Int]): Partition[A] = {
    val nItems = edges.length
    val labels = Array.fill(nItems) { -1 }
    val backedges = Array.fill(nItems) { -1 }
    var nextLabel = 0
    for (i <- 0 until nItems) {
      for (j <- 0 until nItems) backedges(j) = -1
      var float = i
      while ((labels(float) < 0) && (backedges(edges(float)) < 0)) {
        val oldIndex = float
        float = edges(float)
        backedges(float) = oldIndex
      }
      if (labels(float) < 0) {
        while (float != i) {
          labels(float) = nextLabel
          float = backedges(float)
        }
        labels(i) = nextLabel
        nextLabel += 1
      } else {
        val label = labels(float)
        while (float != i) {
          labels(float) = label
          float = backedges(float)
        }
        labels(i) = label
      }
    }
    Partition((i: Int) => samplingModel.sample(), labels)
  }

  private def enumerateEdges(nItems: Int, state: List[Array[Int]], depth: Int): List[Array[Int]] = {
    if (depth == nItems) state
    else {
      Range(0, nItems).map(i =>
        enumerateEdges(nItems, state.map(x => Array(i, x: _*)).toList, depth + 1)).flatten.toList
    }
  }

  private def enumerateEdges(nItems: Int): List[Array[Int]] = enumerateEdges(nItems, List(Array[Int]()), 0)

}

object DDCRP {

  def apply[A](samplingModel: SamplingModel[A], mass: Mass, distance: Distance, decay: DecayFunction) = {
    val similarity = Similarity(distance, decay)
    new DDCRP(samplingModel, mass, similarity)
  }

  def factory[A](samplingModel: SamplingModel[A], massFactory: () => Mass, attractionFactory: () => Attraction) = () => new DDCRP(samplingModel, massFactory(), attractionFactory().similarity)

}

