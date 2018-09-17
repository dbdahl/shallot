package org.ddahl.shallot
package mcmc

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.log
import distribution._
import parameter._
import partition._

object DiscountSampler {

  // Gaussian random walk assuming a Beta(shape1, shape2) prior on discount.
  def gaussianRandomWalk[B <: PartitionModel[A] with HasDiscount[B, A], A](d: HasDiscount[B, A], partition: Partition[A], shape1: Double, shape2: Double, standardDeviation: Double, rdg: RDG): (B, Boolean) = {
    def logPrior(value: Double) = {
      (shape1 - 1) * log(value) + (shape2 - 1) * log(1 - value)
    }
    val current = d.self
    val alpha = rdg.nextGaussian(current.discount.value, standardDeviation)
    if ((alpha <= 0.0) || (alpha >= 1.0)) return (current, false)
    val proposal = d.replaceDiscount(Discount(alpha))
    if ( d.discount.value == 0 ) return (proposal, true)
    val logMHRatio = proposal.logProbability(partition) + logPrior(proposal.discount.value) - current.logProbability(partition) - logPrior(current.discount.value)
    if (log(rdg.nextUniform(0.0, 1.0, false)) <= logMHRatio) (proposal, true)
    else (current, false)
  }

  // Gaussian random walk assuming a Beta(shape1, shape2) prior on discount.
  def spikeAndSlabGaussianRandomWalk[B <: PartitionModel[A] with HasDiscount[B, A], A](d: HasDiscount[B, A], partition: Partition[A], probOf0: Double, shape1: Double, shape2: Double, standardDeviation: Double, rdg: RDG): (B, Boolean) = {
    def logPrior(value: Double) = {
      if ( value != 0.0 ) (shape1 - 1) * log(value) + (shape2 - 1) * log(1 - value)
      else log(probOf0)
    }
    val current = d.self
    if ( current.discount.value == 0.0 ) {
      val alpha = rdg.nextBeta(shape1, shape2)
      val proposal = d.replaceDiscount(Discount(alpha))
      val logMHRatio = proposal.logProbability(partition) - current.logProbability(partition)
      if (log(rdg.nextUniform(0.0, 1.0, false)) <= logMHRatio) (proposal, true)
      else (current, false)
    } else {
      if ( rdg.nextUniform(0.0,1.0) < probOf0 ) {
        val proposal = d.replaceDiscount(Discount(0.0))
        val logMHRatio = proposal.logProbability(partition) - current.logProbability(partition)
        if (log(rdg.nextUniform(0.0, 1.0, false)) <= logMHRatio) (proposal, true)
        else (current, false)
      } else {
        val alpha = rdg.nextGaussian(current.discount.value, standardDeviation)
        if ((alpha < 0.0) || (alpha >= 1.0)) return (current, false)
        val proposal = d.replaceDiscount(Discount(alpha))
        val logMHRatio = proposal.logProbability(partition) - current.logProbability(partition)
        if (log(rdg.nextUniform(0.0, 1.0, false)) <= logMHRatio) (proposal, true)
        else (current, false)
      }
    }
  }

}

