package org.ddahl.shallot.parameter

trait DistanceFunction[T] {

  def apply(x: T, y: T): Double
  val satisfiesTriangleInequality: Boolean

}

object SequentialDistance extends DistanceFunction[Int] {

  def apply(x: Int, y: Int): Double = {
    math.abs(x - y)
  }

  val satisfiesTriangleInequality = true

}

object EuclideanDistance extends DistanceFunction[Array[Double]] {

  def apply(x: Array[Double], y: Array[Double]): Double = {
    val n = x.length
    if (y.length != n) throw new IllegalArgumentException("Length mismatch.")
    math.sqrt(Range(0, n).map(i => (x(i) - y(i)) * (x(i) - y(i))).sum)
  }

  val satisfiesTriangleInequality = true

}

object ManhattanDistance extends DistanceFunction[Array[Double]] {

  def apply(x: Array[Double], y: Array[Double]): Double = {
    val n = x.length
    if (y.length != n) throw new IllegalArgumentException("Length mismatch.")
    Range(0, n).map(i => math.abs(x(i) - y(i))).sum
  }

  val satisfiesTriangleInequality = true

}
