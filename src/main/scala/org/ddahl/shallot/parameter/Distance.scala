package org.ddahl.shallot
package parameter

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }

// Immutable because of methods and companion object constructor

class Distance private (private val x: Array[Array[Double]]) extends Matrix {

  implicit val ordering = CrossCompatibility.doubleOrdering

  val nItems = x.length

  // Just use "lower" triangular part
  def apply(i: Int, k: Int) = if (k < i) x(i)(k) else if (i == k) 0.0 else x(k)(i)

  def remove(i: Int): Distance = {
    val subset = x.patch(i,Seq(),1).map(_.patch(i,Seq(),1))
    Distance(subset,false)
  }

  def defaultPermutation: Permutation = {
    Permutation(Range(0,nItems).map(i => {
      (Range(0,nItems).foldLeft(0.0) { (sum,j) => sum + apply(i,j) },i)
    }).sortWith(_._1 > _._1).map(_._2).toArray)
  }
  
  def max = x.tail.map(_.max).max
  def min = x.tail.map(_.min).min

  private def checkTriangleInequality = {
    val d = this
    val triplets = (0 until nItems - 2).map(i => {
      (0 until nItems - 1).map(j => {
        (0 until nItems - 1).map(k => {
          (i, j, k)
        })
      }).flatten
    }).flatten
    triplets.forall(t => {
      val okay1 = d(t._1, t._2) + d(t._2, t._3) >= d(t._1, t._3)
      val okay2 = d(t._1, t._3) + d(t._3, t._2) >= d(t._1, t._2)
      val okay3 = d(t._2, t._1) + d(t._1, t._3) >= d(t._2, t._3)
      if (!okay1) throw new IllegalArgumentException("Problem with d12 + d23 >= d13 for " + t)
      if (!okay2) throw new IllegalArgumentException("Problem with d13 + d23 >= d12 for " + t)
      if (!okay3) throw new IllegalArgumentException("Problem with d12 + d13 >= d23 for " + t)
      okay1 && okay2 && okay3
    })
  }

}

object Distance {

  def apply(x: Array[Array[Double]], checkTriangleInequality: Boolean): Distance = {
    val nItems = x.length
    val y = new Array[Array[Double]](nItems)
    for (i <- 0 until nItems) {
      y(i) = x(i).take(i)
    }
    for (i <- 0 until nItems) {
      if (y(i).length != i) throw new IllegalArgumentException(s"Row ${i + 1} has length ${y(i).length}, but should be ${i}.")
      if (!y(i).forall(_ > 0.0)) throw new IllegalArgumentException("Not all the elements are strictly positive.")
      if (!y(i).forall(_ < Double.PositiveInfinity)) throw new IllegalArgumentException("Not all the elements are strictly less than infinity.")
    }
    val z = new Distance(y)
    if (checkTriangleInequality) z.checkTriangleInequality
    z
  }

  def apply[T](x: Array[T], norm: DistanceFunction[T]): Distance = {
    val nItems = x.length
    val y = new Array[Array[Double]](nItems)
    for (i <- 0 until nItems) {
      y(i) = new Array[Double](i)
      for (j <- 0 until i) {
        y(i)(j) = norm(x(i), x(j))
      }
    }
    apply(y, !norm.satisfiesTriangleInequality)
  }

  def sample(nItems: Int, rdg: RDG): Distance = {
    val nDimensions = 2
    val x = new Array[Array[Double]](nItems)
    for (i <- 0 until nItems) {
      x(i) = new Array[Double](nDimensions)
      for (j <- 0 until nDimensions) {
        x(i)(j) = rdg.nextGamma(1.0, 1.0)
      }
    }
    apply(x, EuclideanDistance)
  }

}
