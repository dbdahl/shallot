package org.ddahl.shallot
package mcmc

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.log
import distribution._
import parameter._
import partition._

object MassSampler {

  // Method from Escobar & West (1995) putting a Gamma(shape,rate) prior on mass.
  // NOTE: THEY GIVE THEIR ALGORITHM IN TERMS OF A RATE EVEN THOUGH THEY CALL IT A "SCALE"!!!
  def escobarWest[B <: PartitionModel[A], A](d: HasMassEscobarWest[B, A], partition: Partition[A], shape: Double, rate: Double, rdg: RDG): (B, Boolean) = {
    val n = partition.nItems
    val k = partition.nClusters
    val eta = rdg.nextBeta(d.mass.value + 1, n)
    val a = shape
    val b = rate
    val bMinusLogEta = b - log(eta)
    val w1 = a + k - 1
    val w2 = n*bMinusLogEta
    val pi = w1 / (w1 + w2)
    val alpha = if (rdg.nextUniform(0.0, 1.0, true) < pi) rdg.nextGamma(a + k, 1.0 / bMinusLogEta)
    else rdg.nextGamma(a + k - 1, 1.0 / bMinusLogEta)
    (d.replaceMass(Mass(alpha)), true)
  }

  // Gaussian random walk assuming a Gamma(shape,rate) prior on mass.
  def gaussianRandomWalk[B <: PartitionModel[A] with HasMass[B, A], A](d: HasMass[B, A], partition: Partition[A], shape: Double, rate: Double, standardDeviation: Double, rdg: RDG): (B, Boolean) = {
    def logPrior(value: Double) = {
      (shape - 1) * log(value) - rate * value
    }
    val current = d.self
    val alpha = rdg.nextGaussian(d.mass.value, standardDeviation)
    if (alpha <= 0.0) return (current, false)
    val proposal = d.replaceMass(Mass(alpha))
    val logMHRatio = proposal.logProbability(partition) + logPrior(proposal.mass.value) - current.logProbability(partition) - logPrior(current.mass.value)
    if (log(rdg.nextUniform(0.0, 1.0, false)) <= logMHRatio) (proposal, true)
    else (current, false)
  }

}

