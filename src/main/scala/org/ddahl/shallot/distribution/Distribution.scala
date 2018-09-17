package org.ddahl.shallot
package distribution

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.{ exp, pow, round }
import parameter._
import parameter.partition._

trait Distribution[A] {

  def partitions: Iterable[Partition[A]]

  def apply(partition: Partition[A]): Double

  def sample(rdg: RDG): Partition[A]

  def sampler(rdg: RDG) = () => sample(rdg)

  def approximatelyEquals(other: Distribution[A], epsilon: Double): Boolean = {
    val allPartitions = (other.partitions ++ this.partitions).toSet
    val okay = allPartitions.forall(partition => {
      math.abs(this(partition) - other(partition)) <= epsilon
    })
    if (!okay) {
      val sortedPartitions = allPartitions.toList.sortWith(_.toString < _.toString)
      println("<<<<")
      sortedPartitions.foreach(partition => {
        println(("%1.4f %1.4f %s").format(this(partition), other(partition), partition))
      })
      println(">>>>")
    }
    okay
  }

}

private class DistributionWithPartitions[A](private val map: Iterable[(Partition[A], Double)]) extends Distribution[A] {

  def partitions = map.map(_._1)

  def apply(partition: Partition[A]) = map.find(_._1 == partition).getOrElse((partition, 0.0))._2

  def sample(rdg: RDG): Partition[A] = {
    val unif = rdg.nextUniform(0.0, 1.0, true)
    var cumsum = 0.0
    map.find(element => {
      cumsum += element._2
      cumsum > unif
    }).get._1
  }

  override def toString = {
    map.toList.sortWith((x, y) => {
      if (x._1.nClusters < y._1.nClusters) true
      else if (x._1.nClusters > y._1.nClusters) false
      else x._1.toStringTerse < y._1.toStringTerse
    }).map(x => " %1.6f".format(x._2) + " " + x._1.toStringTerse).mkString("\n")
  }

}

private class DistributionWithSubsets[A](private val partition: Partition[A], private val i: Int, private val map: Iterable[(Subset[A], Double)]) extends Distribution[A] {

  def partitions = map.map(x => partition.add(i, x._1))

  def apply(partition: Partition[A]) = map.find(x => partition.add(i, x._1) == this.partition).getOrElse((partition, 0.0))._2

  def sample(rdg: RDG): Partition[A] = {
    val unif = rdg.nextUniform(0.0, 1.0, true)
    var cumsum = 0.0
    val subset = map.find(element => {
      cumsum += element._2
      cumsum > unif
    }).get._1
    partition.add(i, subset)
  }

  override def toString = {
    map.map(x => partition.add(i, x._1).toLabels.mkString(",") + " %1.6f".format(x._2)).mkString("\n")
  }

}

object Distribution {

  def apply[A](elements: Iterable[(Partition[A], Double)], onLogScale: Boolean, normalize: Boolean): Distribution[A] = {
    val x1 = if (onLogScale) {
      elements.map(y => (y._1, exp(y._2)))
    } else elements
    val x2 = if (normalize) {
      val sum = x1.foldLeft(0.0)((sum, y) => sum + y._2)
      x1.map(y => (y._1, y._2 / sum))
    } else x1
    new DistributionWithPartitions(x2)
  }

  def apply[A](partition: Partition[A], i: Int, elements: Iterable[(Subset[A], Double)], onLogScale: Boolean, normalize: Boolean): Distribution[A] = {
    val x1 = if (onLogScale) {
      elements.map(y => (y._1, exp(y._2)))
    } else elements
    val x2 = if (normalize) {
      val sum = x1.foldLeft(0.0)((sum, y) => sum + y._2)
      if (sum == 0.0) {
        val values = elements.map(y => y._2).toArray
        val max = values.max
        val index = values.indexWhere(_ == max)
        val w = new Array[Double](values.length)
        w(index) = 1.0
        elements.map(y => y._1).zip(w)
      } else {
        x1.map(y => (y._1, y._2 / sum))
      }
    } else x1
    new DistributionWithSubsets(partition, i, x2)
  }

}

