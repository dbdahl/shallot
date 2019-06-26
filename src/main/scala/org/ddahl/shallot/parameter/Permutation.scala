package org.ddahl.shallot
package parameter

import partition._

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }

// Immutable

class Permutation private (private val x: Array[Int]) extends Iterable[Int] {

  val nItems = x.length

  def apply(i: Int) = x(i)

  def iterator = x.iterator

  def nPreceeding(i: Int) = xinv(i)

  def antecedents(i: Int) = x.view.slice(0, xinv(i))

  def thin[A](subset: Subset[A], i: Int): Iterator[Int] = new Iterator[Int]() {
    val nPreceedingI = if (i < 0) 0 else xinv(i)
    private var index = -1
    findNext
    private def findNext = do index += 1 while ((index < nItems) && ((index < nPreceedingI) || !subset.contains(x(index))))
    def hasNext = (index < nItems)
    def next = {
      val result = x(index)
      findNext
      result
    }
  }

  def thin[A](subset: Subset[A]): Iterator[Int] = thin(subset, -1)
  
  private def shuffleSubset(selectedIndices: Array[Int], rdg: RDG): Permutation = {
    val k = selectedIndices.length
    val permutation = rdg.nextPermutation(k, k)
    val y = x.clone
    for (i <- 0 until k) {
      y(selectedIndices(permutation(i))) = x(selectedIndices(i))
    }
    new Permutation(y)
  }

  def shuffleSubset(k: Int, rdg: RDG): Permutation = {
    shuffleSubset(rdg.nextPermutation(nItems, k),rdg)
  }

  def shuffleSubsetExcept(k: Int, rdg: RDG, fixedIndices: Set[Int]): Permutation = {
    shuffleSubset(rdg.nextPermutation(nItems, k).filter(i => !fixedIndices.contains(x(i))),rdg)
  }
  
  // Where i is the index that is not available for shuffling.  This is *not* necessarily the ith item in the permutation.
  def shuffleSubsetExcept(k: Int, rdg: RDG, i: Int): Permutation = {
    shuffleSubsetExcept(k,rdg,Set(i))
  }

  // Where i is the index that is removed.  This is *not* necessarily the ith item in the permutation.
  def remove(i: Int): Permutation = Permutation(x.filter(_ != i).map(j => if ( j >= i ) j-1 else j))

  // Where i is the index that is appended.
  def append(i: Int): Permutation = {
    Permutation(x.map(j => if ( j >= i ) j+1 else j) ++ Array(i))
  }

  override def toString = s"(${x.mkString(",")})"

  private val xinv = new Array[Int](nItems)
  for (k <- 0 until nItems) xinv(x(k)) = k

  def inverse: Permutation = Permutation(xinv)
  
}

object Permutation {

  def apply(x: Array[Int]) = {
    var min = Integer.MAX_VALUE
    var max = Integer.MIN_VALUE
    x.foreach(i => {
      if (i < min) min = i
      if (i > max) max = i
    })
    if (min < 0 || max > x.length - 1) throw new IllegalArgumentException("List is not a permutation.")
    val distinct = x.distinct
    if (x.length != distinct.length) throw new IllegalArgumentException("List is not a permutation.")
    new Permutation(distinct)
  }

  def apply(x: Int*): Permutation = apply(x.toArray)

  def natural(nItems: Int) = apply(Range(0, nItems).toArray)

  def enumerate(nItems: Int): List[Permutation] = Range(0, nItems).toList.permutations.toList.map(x => new Permutation(x.toArray))

  def sample(nItems: Int, rdg: RDG) = new Permutation(rdg.nextPermutation(nItems, nItems))

  def factory(x: Array[Int]) = {
    val permutation = apply(x)
    () => permutation
  }

  def factory(permutation: Permutation) = {
    () => permutation
  }

  def factory(nItems: Int, rdg: RDG) = {
    () => new Permutation(rdg.nextPermutation(nItems, nItems))
  }

}
