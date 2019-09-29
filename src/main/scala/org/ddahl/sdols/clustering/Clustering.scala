package org.ddahl.sdols
package clustering

import scala.annotation.tailrec

final class Clustering[A](val nItems: Int, val nClusters: Int, protected val x: Set[Cluster[A]]) extends Iterable[Cluster[A]] {

  final private val checks = Clustering.checks

  if (checks) {
    if (nClusters != x.size) throw new RuntimeException("Internal error")
    if (x.foldLeft(0)((sum, cluster) => sum + cluster.size) != nItems) throw new RuntimeException("Internal error")
  }

  def add(cluster: Cluster[A]): Clustering[A] = {
    if (checks) {
      if (cluster.size == 0) throw new IllegalArgumentException("Cluster is empty.")
      if (contains(cluster)) throw new IllegalArgumentException("Clustering already contains this cluster.")
      for ( oldCluster <- x ) if ( ! oldCluster.x.intersect(cluster.x).isEmpty ) throw new IllegalArgumentException("Cluster contains items that are already clustered.")
    }
    new Clustering(nItems + cluster.size, nClusters + 1, x + cluster)
  }

  def add(i: Int, cluster: Cluster[A]): Clustering[A] = {
    if (checks) {
      if (contains(i)) throw new IllegalArgumentException("Clustering already contains " + i + ".")
    }
    val newCluster = cluster.add(i)
    if (contains(cluster)) new Clustering(nItems + 1, nClusters, x - cluster + newCluster)
    else new Clustering(nItems + 1, nClusters + 1, x + newCluster)
  }

  def remove(cluster: Cluster[A]): Clustering[A] = {
    if (checks) {
      if (!contains(cluster)) throw new IllegalArgumentException("Clustering does not contain this cluster.")
    }
    new Clustering(nItems - cluster.size, nClusters - 1, x - cluster)
  }

  def remove(i: Int, cluster: Cluster[A]): Clustering[A] = {
    if (checks) {
      if (!cluster.contains(i)) throw new IllegalArgumentException("Cluster does not contain " + i + ".")
    }
    if (cluster.size == 1) new Clustering(nItems - 1, nClusters - 1, x - cluster)
    else new Clustering(nItems - 1, nClusters, x - cluster + cluster.remove(i))
  }

  def remove(i: Int): Clustering[A] = {
    val s = x.find(_.contains(i))
    if (checks) {
      if (s.isEmpty) throw new IllegalArgumentException("Clustering does not contain " + i + ".")
    }
    val cluster = s.get
    remove(i, cluster)
  }

  def removeWithCluster(i: Int): (Clustering[A], Cluster[A]) = {
    val s = x.find(_.contains(i))
    if (checks) {
      if (s.isEmpty) throw new IllegalArgumentException("Clustering does not contain " + i + ".")
    }
    val cluster = s.get
    val clusterWithoutI = cluster.remove(i)
    if (cluster.size == 1) {
      (new Clustering(nItems - 1, nClusters - 1, x - cluster), clusterWithoutI)
    } else {
      (new Clustering(nItems - 1, nClusters, x - cluster + clusterWithoutI), clusterWithoutI)
    }
  }

  def replace(func: (Cluster[A]) => A): Clustering[A] = {
    Clustering(map(cluster => {
      cluster.replace(func(cluster))
    }))
  }

  def contains(cluster: Cluster[A]): Boolean = x.contains(cluster)

  def contains(i: Int): Boolean = x.exists(_.contains(i))

  def paired(i: Int, k: Int): Boolean = {
    val clusterOption = x.find(y => y.contains(i) || y.contains(k))
    if ( clusterOption.isEmpty ) false
    else {
      val cluster = clusterOption.get
      cluster.contains(i) && cluster.contains(k)
    }
  }

  def clusterFor(i: Int): Cluster[A] = {
    x.find(_.contains(i)).get
  }

  def parameterFor(i: Int): A = {
    clusterFor(i).parameter
  }

  def pairwiseAllocationMatrix: Array[Array[Int]] = {
    val r = Array.ofDim[Int](nItems, nItems)
    x.foreach(s => {
      val indices = s.toList.filter(_ < nItems).sortWith(_ > _)
      var x = indices
      while (!x.isEmpty) {
        val rr = r(x.head)
        rr(x.head) += 1
        x.tail.foreach(rr(_) += 1)
        x = x.tail
      }
    })
    // Symmetrize
    var i = 0
    while ( i < nItems ) {
      val ri = r(i)
      var j = i + 1
      while ( j < nItems ) {
        ri(j) = r(j)(i)
        j += 1
      }
      i += 1
    }
    r
  }

  override def equals(other: Any): Boolean = other match {
    case that: Clustering[A] =>
      if (that.nItems != nItems) false
      else that.x == x
    case _ => false
  }

  override def hashCode: Int = x.hashCode

  def iterator = x.iterator

  override def size = nClusters

  override def toString: String = "{" + map(_.toString).toList.sortWith(_ < _).mkString(",") + "}"

  def toStringTerse: String = "{" + map(_.toStringTerse).toList.sortWith(_ < _).mkString(",") + "}"

  def entropy: Double = {
    val rates = map(_.size.asInstanceOf[Double] / nItems).toList.sortWith(_ > _)
    rates.foldLeft(0.0)((s, p) => { s - (if (p > 0.0) p * math.log(p) else 0.0) })
  }

  def toLabels: Array[Int] = {
    val result = new Array[Int](nItems)
    var label = 0
    toList.sortWith(_.min < _.min).foreach(cluster => {
      cluster.foreach(i => result(i) = label)
      label += 1
    })
    result
  }

  def toLabelsWithParameters: (Array[Int], List[A]) = {
    val result = new Array[Int](nItems)
    var resultParameters = List[A]()
    var label = 0
    toList.sortWith(_.min < _.min).foreach(cluster => {
      cluster.foreach(i => result(i) = label)
      resultParameters = cluster.parameter :: resultParameters
      label += 1
    })
    (result, resultParameters.reverse)
  }

  def write(objOutputStream: java.io.ObjectOutputStream): Unit = {
    objOutputStream.writeInt(nClusters)
    iterator.foreach(_.write(objOutputStream))
  }

}

object Clustering {

  var checks = false

  def empty[A](): Clustering[A] = new Clustering(0, 0, Set[Cluster[A]]())

  def apply[A](sampler: () => A, nItems: Int, allTogether: Boolean): Clustering[A] = {
    if (allTogether) {
      Clustering(Cluster(sampler(), Range(0, nItems)))
    } else {
      Clustering(Range(0, nItems).map(i => Cluster(sampler(), i)): _*)
    }
  }

  @tailrec
  private def makeClustering[A](sampler: Int => A, labelsWithIndex: Iterable[(Int, Int)], list: List[Cluster[A]]): List[Cluster[A]] = {
    val label = labelsWithIndex.head._1
    val (left, right) = labelsWithIndex.partition(_._1 == label)
    val longerList = Cluster(sampler(label), left.map(_._2)) +: list
    if (right.isEmpty) longerList
    else makeClustering(sampler, right, longerList)
  }

  def apply[A](sampler: Int => A, labels: Iterable[Int]): Clustering[A] = {
    if (labels.isEmpty) throw new IllegalArgumentException("Labels may not by empty.")
    apply(makeClustering(sampler, labels.zipWithIndex, List[Cluster[A]]()))
  }

  def apply(labels: Array[Int]): Clustering[Null] = apply((i: Int) => null, labels)

  def apply[A](i: Cluster[A]*): Clustering[A] = apply(i)

  def apply[A](i: Iterable[Cluster[A]]): Clustering[A] = {
    val nItems = i.foldLeft(0)((sum, cluster) => sum + cluster.size)
    new Clustering(nItems, i.size, i.toSet)
  }

  def apply[A](objInputStream: java.io.ObjectInputStream): Clustering[A] = {
    val nClusters = objInputStream.readInt()
    val seq = Seq.fill(nClusters) { Cluster[A](objInputStream) }
    apply(seq)
  }

  def enumerate[A](sampler: () => A, nItems: Int): List[Clustering[A]] = {
    var cache = List[Clustering[A]]()
    def engine(clustering: Clustering[A], nextItem: Int, nItems: Int): Unit = {
      if (nextItem == nItems) {
        cache = clustering +: cache
      } else {
        clustering.foreach(cluster => {
          engine(clustering.add(nextItem, cluster), nextItem + 1, nItems)
        })
        engine(clustering.add(nextItem, Cluster.empty(sampler())), nextItem + 1, nItems)
      }
    }
    engine(empty[A](), 0, nItems)
    cache
  }

}

