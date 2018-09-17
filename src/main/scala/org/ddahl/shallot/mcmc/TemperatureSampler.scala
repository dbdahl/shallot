package org.ddahl.shallot
package mcmc

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.log
import distribution._
import parameter._
import partition._

object TemperatureSampler {

  // Gaussian random walk assuming a Gamma(shape,rate) prior on temperature.
  def gaussianRandomWalk[B <: PartitionModel[A] with HasAttraction[B, A], A](d: HasAttraction[B, A], partition: Partition[A], shape: Double, rate: Double, standardDeviation: Double, rdg: RDG): (B, Boolean) = {
    def logPrior(value: Double) = {
      (shape - 1) * log(value) - rate * value
    }
    val current = d.self
    val t = rdg.nextGaussian(d.attraction.decay.temperature, standardDeviation)
    if (t < 0.0) return (current, false)
    val proposal = d.replaceAttraction(d.attraction.replaceDecay(d.attraction.decay.update(t)))
    val logMHRatio = proposal.logProbability(partition) + logPrior(proposal.attraction.decay.temperature) - current.logProbability(partition) - logPrior(current.attraction.decay.temperature)
    if (log(rdg.nextUniform(0.0, 1.0, false)) <= logMHRatio) (proposal, true)
    else (current, false)
  }

}

