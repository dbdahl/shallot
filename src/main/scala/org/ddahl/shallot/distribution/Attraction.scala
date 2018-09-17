package org.ddahl.shallot
package distribution

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.pow
import parameter._
import parameter.decay._
import parameter.partition._

// Immutable
class Attraction private (val distance: Distance, val permutation: Permutation, val decay: DecayFunction, val similarity: Similarity, private val denominator: Array[Double]) extends Matrix {

  val nItems = distance.nItems

  private val cache = {
    val x = new Array[Array[Double]](nItems)
    for (i <- 0 until nItems) {
      val y = new Array[Double](nItems)
      val ii = permutation.nPreceeding(i)
      for (k <- 0 until nItems) {
        y(k) = if (ii <= permutation.nPreceeding(k)) 0.0
        else similarity(i, k) / denominator(i)
      }
      x(i) = y
    }
    x
  }

  def replaceDistance(newDistance: Distance) = {
    if (distance.nItems != nItems) throw new IllegalArgumentException("Inconsistent number of items.")
    val newSimilarity = Similarity(distance, decay)
    val newDenominator = Attraction.computeDenominator(newSimilarity, permutation)
    new Attraction(newDistance, permutation, decay, newSimilarity, newDenominator)
  }

  def replacePermutation(newPermutation: Permutation) = {
    if (newPermutation.nItems != nItems) throw new IllegalArgumentException("Inconsistent number of items.")
    val newDenominator = Attraction.computeDenominator(similarity, newPermutation)
    new Attraction(distance, newPermutation, decay, similarity, newDenominator)
  }

  def replaceDecay(newDecay: DecayFunction) = {
    val similarity = Similarity(distance, newDecay)
    val newDenominator = Attraction.computeDenominator(similarity, permutation)
    new Attraction(distance, permutation, newDecay, similarity, newDenominator)
  }

  def apply(i: Int, k: Int): Double = cache(i)(k)

  def apply[A](i: Int, subset: Subset[A]): Double = {
    subset.foldLeft(0.0)((sum, k) => sum + apply(i, k))
  }
  
  def remove[A](i: Int): (Array[Double],Attraction) = {
    val d = Range(0,nItems).map(distance(i,_)).patch(i,Seq(),1).toArray
    val newDistance = distance.remove(i)
    val newPermutation = permutation.remove(i)
    val newSimilarity = Similarity(newDistance, decay)
    val newDenominator = Attraction.computeDenominator(newSimilarity, newPermutation)
    (d, new Attraction(newDistance, newPermutation, decay, newSimilarity, newDenominator))
  }

}

object Attraction {

  private def computeDenominator(similarity: Similarity, permutation: Permutation): Array[Double] = {
    val nItems = permutation.nItems
    val denominator = new Array[Double](nItems)
    for (i <- 0 until nItems) {
      denominator(i) = permutation.antecedents(i).foldLeft(0.0)((x, k) => x + similarity(i, k)) / permutation.nPreceeding(i)
    }
    denominator
  }

  def constant(nItems: Int) = {
    val x = Array.fill(nItems, nItems)(1.0)
    val distance = Distance(x, false)
    val permutation = Permutation.natural(nItems)
    val decay = ReciprocalDecayFactory(0.0)
    apply(distance, permutation, decay)
  }

  def apply(distance: Distance, permutation: Permutation, decay: DecayFunction) = {
    if (distance.nItems != permutation.nItems) throw new IllegalArgumentException("Inconsistent number of items.")
    val similarity = Similarity(distance, decay)
    val denominator = computeDenominator(similarity, permutation)
    new Attraction(distance, permutation, decay, similarity, denominator)
  }

  def factory(nItems: Int) = {
    val attraction = constant(nItems)
    () => attraction
  }

  def factory(attraction: Attraction) = {
    () => attraction
  }

  def factory(distance: Distance, permutation: Permutation, decay: DecayFunction) = {
    val attraction = Attraction(distance, permutation, decay)
    () => attraction
  }

  def factory(distance: Distance, permutation: Permutation, decayFactory: () => DecayFunction) = {
    () => Attraction(distance, permutation, decayFactory())
  }

  def factory(distance: Distance, permutationFactory: () => Permutation, decay: DecayFunction) = { // Reuse similarity
    if (distance.nItems != permutationFactory().nItems) throw new IllegalArgumentException("Inconsistent number of items.")
    val similarity = Similarity(distance, decay)
    () => {
      val permutation = permutationFactory()
      val denominator = computeDenominator(similarity, permutation)
      new Attraction(distance, permutation, decay, similarity, denominator)
    }
  }

  def factory(distance: Distance, permutationFactory: () => Permutation, decayFactory: () => DecayFunction) = {
    () => Attraction(distance, permutationFactory(), decayFactory())
  }

}

