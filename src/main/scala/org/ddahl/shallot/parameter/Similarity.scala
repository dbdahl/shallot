package org.ddahl.shallot
package parameter

import decay.DecayFunction
import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }

// Immutable because of methods and companion object constructor

class Similarity private (private val x: Array[Array[Double]]) extends Matrix {

  val nItems = x.length

  // Just use "lower" triangular part
  def apply(i: Int, k: Int) = if (k < i) x(i)(k) else if (i == k) 0.0 else x(k)(i)

  // Most similar items are first
  def defaultPermutation = Permutation(x.map(_.sum).zipWithIndex.sortWith(_._1 > _._1).map(_._2))

}

object Similarity {

  def apply(x: Array[Array[Double]]) = {
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
    new Similarity(y)
  }

  def apply(distance: Distance, decay: DecayFunction) = {
    val nItems = distance.nItems
    val y = new Array[Array[Double]](nItems)
    for (i <- 0 until nItems) {
      val z = new Array[Double](i)
      for (j <- 0 until i) {
        z(j) = decay(distance(i, j))
      }
      y(i) = z
    }
    for (i <- 0 until nItems) {
      if (!y(i).forall(_ > 0.0)) throw new IllegalArgumentException("Not all the elements are strictly positive.")
      if (!y(i).forall(_ < Double.PositiveInfinity)) throw new IllegalArgumentException("Not all the elements are strictly less than infinity.")
    }
    new Similarity(y)
  }

}

