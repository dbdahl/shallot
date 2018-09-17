package org.ddahl.shallot
package mcmc

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.log
import distribution._
import parameter._
import partition._

object PermutationSampler {

  // Propose a shuffling of k randomly selected items in the permutation of n items
  def update[B <: PartitionModel[A] with HasAttraction[B, A], A](d: HasAttraction[B, A], partition: Partition[A], k: Int, rdg: RDG, fixedIndices: Set[Int]): (B, Boolean) = {
    if (k < 2) throw new IllegalArgumentException("k must be at least 2.")
    val current = d.self
    val proposal = d.replaceAttraction(current.attraction.replacePermutation(current.attraction.permutation.shuffleSubsetExcept(k, rdg, fixedIndices)))
    val logMHRatio = proposal.logProbability(partition) - current.logProbability(partition)
    if (log(rdg.nextUniform(0.0, 1.0, false)) <= logMHRatio) (proposal, true)
    else (current, false)
  }

}

