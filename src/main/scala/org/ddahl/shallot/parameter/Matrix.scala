package org.ddahl.shallot
package parameter

trait Matrix {

  val nItems: Int
  def apply(i: Int, j: Int): Double
  
  def row(i: Int): Array[Double] = {
    Range(0,nItems).map(apply(i,_)).toArray
  }

  def column(j: Int): Array[Double] = {
    Range(0,nItems).map(apply(_,j)).toArray
  }
  
  def toArray = {
    val x = Array.ofDim[Double](nItems, nItems)
    for (i <- 0 until nItems) {
      for (j <- 0 until nItems) {
        x(i)(j) = apply(i, j)
      }
    }
    x
  }

  override def toString = {
    Range(0, nItems).map(i =>
      Range(0, nItems).map(j =>
        "%1.4f".format(apply(i, j))).mkString(" ")).mkString(System.lineSeparator)
  }

}

