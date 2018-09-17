package org.ddahl.shallot
package parameter
package partition

import scala.util.control.Breaks._

object MinBinder {

  private def allocate[A](i: Int, partition: Partition[A], pp: PairwiseProbability) = {
    val absorbingSubset = partition.minBy(subset => {
      subset.foldLeft(0.0)((sum, k) => sum + (0.5 - pp(i, k)))
    })
    partition.add(i, absorbingSubset)
  }

  private def dropSmallSubsets[A](partition: Partition[A], maxSubsets: Int, pp: PairwiseProbability) = {
    val (toBeKept, toBeAbsorbed) = partition.toList.sortWith(_.size > _.size).splitAt(maxSubsets)
    var p = Partition[A](toBeKept)
    toBeAbsorbed.map(_.iterator).flatten.foreach(i => p = allocate(i, p, pp))
    p
  }

  private def reallocate[A](pp: PairwiseProbability, partition: Partition[A], maxScans: Int) = {
    val nItems = partition.nItems
    var p = partition
    var old: Partition[A] = null
    var scanCount = 0
    while ((scanCount < maxScans) && (p != old)) {
      old = p
      for (i <- 0 until nItems) {
        p = allocate(i, p.remove(i), pp)
      }
      scanCount += 1
    }
    p
  }

  private def search[A](partitions: List[Partition[A]], maxSubsets: Int, maxScans: Int, pp: PairwiseProbability): (Partition[A], Double) = {
    val nItems = pp.nItems
    var ppMin = partitions.head
    var ppSumMin = Double.MaxValue
    partitions.foreach(partition => {
      var p = partition
      if (maxSubsets > 0) p = dropSmallSubsets(p, maxSubsets, pp)
      if (maxScans > 0) p = reallocate(pp, p, maxScans)
      var ppSum = 0.0
      p.foreach(subset => {
        val indices = subset.toArray
        for (i <- 0 until indices.length) {
          val ii = indices(i)
          for (k <- 0 until i) {
            val kk = indices(k)
            ppSum += (0.5 - pp(ii, kk))
          }
        }
      })
      if (ppSum < ppSumMin) {
        ppSumMin = ppSum
        ppMin = p
      }
    })
    (ppMin, ppSumMin)
  }

  def apply[A](pp: PairwiseProbability, partitions: List[Partition[A]], maxSubsets: Int, maxScans: Int, parallel: Boolean): Partition[A] = {
    if (partitions.isEmpty) throw new IllegalArgumentException("List must not be empty")
    if (parallel) {
      val nCores = Runtime.getRuntime.availableProcessors
      val size = partitions.size
      val lists = partitions.grouped((size / nCores) + 1).toList.par
      lists.map(x => search(x, maxSubsets, maxScans, pp)).minBy(_._2)._1
    } else {
      search(partitions, maxSubsets, maxScans, pp)._1
    }
  }

}

