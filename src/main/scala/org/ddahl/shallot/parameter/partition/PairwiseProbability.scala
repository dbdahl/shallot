package org.ddahl.shallot
package parameter
package partition

class PairwiseProbability(x: Array[Array[Double]]) extends Matrix {

  val nItems = x.length

  // Just use "lower" triangular part
  def apply(i: Int, k: Int) = if (k < i) x(i)(k) else if (i == k) 1.0 else x(k)(i)

  private def overlap[A](subset1: Subset[A], subset2: Subset[A]): Double = {
    subset1.map(i => subset2.map(j => apply(i, j)).sum).sum / (subset1.size * subset2.size)
  }

  def confidenceComputations[A](partition: Partition[A]): (Array[Int], Array[Double], Array[Array[Double]], Array[Int], Array[Int]) = {
    val matrix = scala.collection.mutable.HashMap[Tuple2[Subset[A], Subset[A]], Double]()
    partition.map(s1 => partition.map(s2 => matrix((s1, s2)) = overlap(s1, s2)))
    val subsetsWithSumOverlap = partition.toList.map(s1 => (s1, partition.map(s2 => matrix(s1, s2)).sum))
    val sortedClusters = subsetsWithSumOverlap.sortWith((t1, t2) => t1._2 < t2._2).map(_._1)
    val map = sortedClusters.zipWithIndex.toMap
    val matrixOutSmall = Array.ofDim[Double](partition.nClusters, partition.nClusters)
    partition.map(c1 => partition.map(c2 => matrixOutSmall(map(c1))(map(c2)) = matrix((c1, c2))))
    val confidence = new Array[Double](nItems)
    for (i <- 0 until partition.nItems) {
      val subset = partition.clusterFor(i)
      confidence(i) = subset.map(j => apply(i, j)).sum / subset.size
    }
    val UNINITIALIZED = -1
    val exemplar = Array.fill(map.size)({ UNINITIALIZED })
    val order = Array.range(0, nItems).sortWith((i, j) => {
      val io = map(partition.clusterFor(i))
      val jo = map(partition.clusterFor(j))
      if (io < jo) {
        true
      } else if (io > jo) {
        false
      } else {
        val iBigger = confidence(i) > confidence(j)
        if (iBigger) {
          if ((exemplar(io) == UNINITIALIZED) || (confidence(exemplar(io)) < confidence(i))) exemplar(io) = i
        } else {
          if ((exemplar(jo) == UNINITIALIZED) || (confidence(exemplar(jo)) < confidence(j))) exemplar(jo) = j
        }
        iBigger
      }
    })
    val labels = Array.range(0, nItems).map(i => map(partition.clusterFor(i)))
    (labels, confidence, matrixOutSmall, order, exemplar)
  }

}

object PairwiseProbability {

  private def tally[A](nItems: Int, partitions: List[Partition[A]]): Array[Array[Int]] = {
    val counts = new Array[Array[Int]](nItems)
    for (i <- 0 until nItems) counts(i) = new Array[Int](i)
    partitions.foreach(partition => {
      partition.foreach(subset => {
        val indices = subset.toArray
        for (i <- 0 until indices.length) {
          val ii = indices(i)
          for (k <- 0 until i) {
            val kk = indices(k)
            if (kk < ii) counts(ii)(kk) += 1
            else counts(kk)(ii) += 1
          }
        }
      })
    })
    counts
  }

  private def reduce(nItems: Int, allCounts: Iterable[Array[Array[Int]]]): Array[Array[Int]] = {
    if (allCounts.isEmpty) {
      val totals = new Array[Array[Int]](nItems)
      for (i <- 0 until nItems) totals(i) = new Array[Int](i)
      totals
    } else {
      val totals = allCounts.head
      allCounts.tail.foreach(counts => {
        for (i <- 0 until nItems) {
          val ti = totals(i)
          val ci = counts(i)
          for (k <- 0 until i) {
            ti(k) += ci(k)
          }
        }
      })
      totals
    }
  }

  def apply[A](partitions: List[Partition[A]], parallel: Boolean) = {
    if (partitions.isEmpty) throw new IllegalArgumentException("List must not be empty")
    val size = partitions.size
    val nItems = partitions.head.nItems
    val totals = if (parallel) {
      val nCores = Runtime.getRuntime.availableProcessors
      val lists = partitions.grouped((size / nCores) + 1).toList.par
      reduce(nItems, lists.map(x => tally(nItems, x)).toList)
    } else {
      tally(nItems, partitions)
    }
    val proportions = new Array[Array[Double]](nItems)
    val sizeAsDouble = size.toDouble
    for (i <- 0 until nItems) {
      proportions(i) = new Array[Double](i)
      val pi = proportions(i)
      val ti = totals(i)
      for (k <- 0 until i) {
        pi(k) += ti(k) / sizeAsDouble
      }
    }
    new PairwiseProbability(proportions)
  }

  def apply(probabilities: Array[Array[Double]]) = {
    val nItems = probabilities.length
    val proportions = new Array[Array[Double]](nItems)
    for (i <- 0 until nItems) {
      proportions(i) = new Array[Double](i)
      val pi = proportions(i)
      val pr = probabilities(i)
      for (k <- 0 until i) {
        val v = pr(k)
        if ((v < 0.0) || (v > 1.0)) throw new IllegalArgumentException("Probabilities must be in [0.0,1.0].")
        pi(k) = pr(k)
      }
    }
    new PairwiseProbability(proportions)
  }

}

