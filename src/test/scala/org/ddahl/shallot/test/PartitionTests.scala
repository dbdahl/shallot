package org.ddahl.shallot
package test

import org.scalatest.FlatSpec

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import parameter._
import parameter.{ NullSamplingModel => q }
import parameter.decay._
import parameter.partition._
import distribution._
import mcmc._

abstract class EmpericalMatchesTheoretical(model: PartitionModel[NullSamplingModel.tipe]) extends FlatSpec {

  val nItems = 4
  val nSamples = 100000
  val epsilon = 0.005

  val distTheoretical = model.tabulate(nItems)

  "Forward sampling based emperical probabilities" should "be within Monte Carlo error of theoretical probabilities" in {
    val rdg = new RDG()
    val sampler = model.forwardSampler(nItems, rdg)
    val db = DistributionBuilder[NullSamplingModel.tipe]()
    for (r <- 0 until nSamples) db.process(sampler(), 1)
    val distEmperical = db.toDistribution
    assert(distEmperical.approximatelyEquals(distTheoretical, epsilon))
  }

  "MCMC based emperical probabilities" should "be within Monte Carlo error of theoretical probabilities" in {
    val rdg = new RDG()
    var partition = Partition(() => q.sample(), nItems, true)
    val db = DistributionBuilder[NullSamplingModel.tipe]()
    for (r <- 0 until nSamples) {
      partition = AuxiliaryGibbsSampler(partition, q, model, rdg)._1
      db.process(partition, 1)
    }
    val distEmperical = db.toDistribution
    assert(distEmperical.approximatelyEquals(distTheoretical, epsilon))
  }

}

abstract class SpecialCases(mass: Mass) extends FlatSpec {

  val nItems = 4
  val epsilon = 1e-10

  val ewens = Ewens(q, mass).tabulate(nItems)

  "EwensPitman with discount=0" should "be equal to Ewens" in {
    val other = EwensPitman(q, mass, Discount(0.0)).tabulate(nItems)
    assert(ewens.approximatelyEquals(other, epsilon))
  }

  "EwensAttraction with uniform distance" should "be equal to Ewens" in {
    val attraction = Attraction.constant(4)
    val other = EwensAttraction(q, mass, attraction).tabulate(nItems)
    assert(ewens.approximatelyEquals(other, epsilon))
  }

  "EwensAttraction with temperature=0" should "be equal to Ewens" in {
    val distance = Distance(Array(
      Array[Double](),
      Array(1.0),
      Array(2.0, 5.0),
      Array(3.0, 1.0, 4.0)), false)
    val permutation = Permutation(2, 0, 1, 3)
    val attraction = Attraction(distance, permutation, ReciprocalDecayFactory(1e-10))
    val other = EwensAttraction(q, mass, attraction).tabulate(nItems)
    assert(ewens.approximatelyEquals(other, epsilon))
  }

}

object Defaults {
  val smallMass = Mass(0.8)
  val largeMass = Mass(2.6)
  val smallDiscount = Discount(0.1)
  val largeDiscount = Discount(0.5)
  val smallTemperature = ReciprocalDecayFactory(0.5)
  val largeTemperature = ReciprocalDecayFactory(3.0)
  val distance1 = Distance(Array(
    Array[Double](),
    Array[Double](1),
    Array[Double](2, 5),
    Array[Double](3, 1, 4)), false)
  val distance2 = Distance(Array(
    Array[Double](),
    Array[Double](10),
    Array[Double](1, 2),
    Array[Double](1, 6, 2)), false)
  val permutation1 = Permutation(0, 1, 2, 3)
  val permutation2 = Permutation(2, 0, 3, 1)
  val attraction1 = Attraction(distance1, permutation1, smallTemperature)
  val attraction2 = Attraction(distance2, permutation1, smallTemperature)
  val attraction3 = Attraction(distance1, permutation2, smallTemperature)
  val attraction4 = Attraction(distance2, permutation2, smallTemperature)
  val attraction5 = Attraction(distance1, permutation1, largeTemperature)
  val attraction6 = Attraction(distance2, permutation1, largeTemperature)
  val attraction7 = Attraction(distance1, permutation2, largeTemperature)
  val attraction8 = Attraction(distance2, permutation2, largeTemperature)
}

import Defaults._

class EwensSmallMassSpec extends EmpericalMatchesTheoretical(Ewens(q, smallMass))
class EwensLargeMassSpec extends EmpericalMatchesTheoretical(Ewens(q, largeMass))

class EwensPitmanSmallMassSmallDiscountSpec extends EmpericalMatchesTheoretical(EwensPitman(q, smallMass, smallDiscount))
class EwensPitmanSmallMassLargeDiscountSpec extends EmpericalMatchesTheoretical(EwensPitman(q, smallMass, largeDiscount))
class EwensPitmanLargeMassSmallDiscountSpec extends EmpericalMatchesTheoretical(EwensPitman(q, largeMass, smallDiscount))
class EwensPitmanLargeMassLargeDiscountSpec extends EmpericalMatchesTheoretical(EwensPitman(q, largeMass, largeDiscount))

class EwensAttraction11 extends EmpericalMatchesTheoretical(EwensAttraction(q, smallMass, attraction1))
class EwensAttraction21 extends EmpericalMatchesTheoretical(EwensAttraction(q, largeMass, attraction1))
class EwensAttraction12 extends EmpericalMatchesTheoretical(EwensAttraction(q, smallMass, attraction2))
class EwensAttraction22 extends EmpericalMatchesTheoretical(EwensAttraction(q, largeMass, attraction2))
class EwensAttraction13 extends EmpericalMatchesTheoretical(EwensAttraction(q, smallMass, attraction3))
class EwensAttraction23 extends EmpericalMatchesTheoretical(EwensAttraction(q, largeMass, attraction3))
class EwensAttraction14 extends EmpericalMatchesTheoretical(EwensAttraction(q, smallMass, attraction4))
class EwensAttraction24 extends EmpericalMatchesTheoretical(EwensAttraction(q, largeMass, attraction4))
class EwensAttraction15 extends EmpericalMatchesTheoretical(EwensAttraction(q, smallMass, attraction5))
class EwensAttraction25 extends EmpericalMatchesTheoretical(EwensAttraction(q, largeMass, attraction5))
class EwensAttraction16 extends EmpericalMatchesTheoretical(EwensAttraction(q, smallMass, attraction6))
class EwensAttraction26 extends EmpericalMatchesTheoretical(EwensAttraction(q, largeMass, attraction6))
class EwensAttraction17 extends EmpericalMatchesTheoretical(EwensAttraction(q, smallMass, attraction7))
class EwensAttraction27 extends EmpericalMatchesTheoretical(EwensAttraction(q, largeMass, attraction7))
class EwensAttraction18 extends EmpericalMatchesTheoretical(EwensAttraction(q, smallMass, attraction7))
class EwensAttraction28 extends EmpericalMatchesTheoretical(EwensAttraction(q, largeMass, attraction7))

class EwensPitmanAttraction11 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, smallMass, smallDiscount, attraction1))
class EwensPitmanAttraction21 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, largeMass, smallDiscount, attraction1))
class EwensPitmanAttraction12 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, smallMass, smallDiscount, attraction2))
class EwensPitmanAttraction22 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, largeMass, smallDiscount, attraction2))
class EwensPitmanAttraction13 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, smallMass, smallDiscount, attraction3))
class EwensPitmanAttraction23 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, largeMass, smallDiscount, attraction3))
class EwensPitmanAttraction14 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, smallMass, smallDiscount, attraction4))
class EwensPitmanAttraction24 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, largeMass, smallDiscount, attraction4))
class EwensPitmanAttraction15 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, smallMass, smallDiscount, attraction5))
class EwensPitmanAttraction25 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, largeMass, smallDiscount, attraction5))
class EwensPitmanAttraction16 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, smallMass, smallDiscount, attraction6))
class EwensPitmanAttraction26 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, largeMass, smallDiscount, attraction6))
class EwensPitmanAttraction17 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, smallMass, smallDiscount, attraction7))
class EwensPitmanAttraction27 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, largeMass, smallDiscount, attraction7))
class EwensPitmanAttraction18 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, smallMass, smallDiscount, attraction7))
class EwensPitmanAttraction28 extends EmpericalMatchesTheoretical(EwensPitmanAttraction(q, largeMass, smallDiscount, attraction7))

class SpecialCasesSmallMass extends SpecialCases(smallMass)
class SpecialCasesLargeMass extends SpecialCases(largeMass)

class UpdateMass extends FlatSpec {

  type Q = NullSamplingModel.tipe

  val nSamples = 100000
  val epsilon = 0.01

  val rdg = new RDG()

  def makeDistributions[R <: PartitionModel[Q] with HasMass[R, Q]](shape: Double, rate: Double, rWSD: Double, useRandomWalk: Boolean, distribution: R, factory: () => R) = {
    val db2 = DistributionBuilder[Q]()
    var p = Partition(Subset(q.sample(), 0, 1, 2, 3))
    var d = distribution
    for (rep <- 0 until nSamples) {
      p = AuxiliaryGibbsSampler(p, q, d, rdg)._1
      d = if (useRandomWalk) MassSampler.gaussianRandomWalk(d, p, shape, rate, rWSD, rdg)._1
      else MassSampler.escobarWest(d.asInstanceOf[HasMassEscobarWest[R, Q]], p, shape, rate, rdg)._1
      db2.process(p, 1)
    }
    val sampler1 = PartitionModel.forwardSampler(p.nItems, factory)
    val db1 = DistributionBuilder[Q]()
    sampler1(nSamples, rdg).foreach(p => db1.process(p, 1))
    (db1.toDistribution, db2.toDistribution)
  }

  "Equal probabilities subject to Monte Carlo error" should "hold for setting 1" in {
    val attraction = attraction1
    val shape = 2.0
    val rate = 4.0
    val rWSD = 1.0
    for (useRandomWalk <- List(true, false)) {
      val (d1, d2) = makeDistributions(shape, rate, rWSD, useRandomWalk, EwensAttraction(q, Mass(3.0), attraction), EwensAttraction.factory(q, Mass.factory(shape, rate, rdg), Attraction.factory(attraction)))
      assert(d1.approximatelyEquals(d2, epsilon))
    }
  }

  "Equal probabilities subject to Monte Carlo error" should "hold for setting 2" in {
    val attraction = attraction2
    val shape = 3.0
    val rate = 1.0
    val rWSD = 0.5
    for (useRandomWalk <- List(true, false)) {
      val (d1, d2) = makeDistributions(shape, rate, rWSD, useRandomWalk, EwensAttraction(q, Mass(3.0), attraction), EwensAttraction.factory(q, Mass.factory(shape, rate, rdg), Attraction.factory(attraction)))
      assert(d1.approximatelyEquals(d2, epsilon))
    }
  }

}

class UpdateDiscount extends FlatSpec {

  type Q = NullSamplingModel.tipe

  val nSamples = 100000
  val epsilon = 0.01

  val rdg = new RDG()

  def makeDistributions[R <: PartitionModel[Q] with HasDiscount[R, Q]](shape1: Double, shape2: Double, rWSD: Double, distribution: R, factory: () => R) = {
    val db2 = DistributionBuilder[Q]()
    var p = Partition(Subset(q.sample(), 0, 1, 2, 3))
    var d = distribution
    for (rep <- 0 until nSamples) {
      p = AuxiliaryGibbsSampler(p, q, d, rdg)._1
      d = DiscountSampler.gaussianRandomWalk(d, p, shape1, shape2, rWSD, rdg)._1
      db2.process(p, 1)
    }
    val sampler1 = PartitionModel.forwardSampler(p.nItems, factory)
    val db1 = DistributionBuilder[Q]()
    sampler1(nSamples, rdg).foreach(p => db1.process(p, 1))
    (db1.toDistribution, db2.toDistribution)
  }

  "Equal probabilities subject to Monte Carlo error" should "hold for setting 1" in {
    val mass = Mass(1.4)
    val discount = Discount(0.05)
    val shape1 = 2.0
    val shape2 = 4.0
    val rWSD = 1.0
    for (useRandomWalk <- List(true, false)) {
      val (d1, d2) = makeDistributions(shape1, shape2, rWSD, EwensPitman(q, mass, discount), EwensPitman.factory(q, Mass.factory(mass), Discount.factory(shape1, shape2, rdg)))
      assert(d1.approximatelyEquals(d2, epsilon))
    }
  }

  "Equal probabilities subject to Monte Carlo error" should "hold for setting 2" in {
    val mass = Mass(0.8)
    val discount = Discount(0.15)
    val shape1 = 3.0
    val shape2 = 3.0
    val rWSD = 0.5
    for (useRandomWalk <- List(true, false)) {
      val (d1, d2) = makeDistributions(shape1, shape2, rWSD, EwensPitman(q, mass, discount), EwensPitman.factory(q, Mass.factory(mass), Discount.factory(shape1, shape2, rdg)))
      assert(d1.approximatelyEquals(d2, epsilon))
    }
  }

}

class UpdateTemperature extends FlatSpec {

  type Q = NullSamplingModel.tipe

  val nSamples = 100000
  val epsilon = 0.005

  val rdg = new RDG()

  def makeDistributions(shape: Double, rate: Double, rWSD: Double, distribution: EwensAttraction[Q], factory: () => EwensAttraction[Q]) = {
    val db2 = DistributionBuilder[Q]()
    var p = Partition(Subset(q.sample(), 0, 1, 2, 3))
    var d = distribution
    for (rep <- 0 until nSamples) {
      p = AuxiliaryGibbsSampler(p, q, d, rdg)._1
      d = TemperatureSampler.gaussianRandomWalk(d, p, shape, rate, rWSD, rdg)._1
      db2.process(p, 1)
    }
    val sampler1 = PartitionModel.forwardSampler(p.nItems, factory)
    val db1 = DistributionBuilder[Q]()
    sampler1(nSamples, rdg).foreach(p => db1.process(p, 1))
    (db1.toDistribution, db2.toDistribution)
  }

  "Equal probabilities subject to Monte Carlo error" should "hold for setting 1" in {
    val mass = Mass(3.0)
    val nItems = 4
    val distance = Distance.sample(nItems, rdg)
    val permutation = Permutation.sample(nItems, rdg)
    val attraction = Attraction(distance, permutation, ReciprocalDecayFactory(1.0))
    val shape = 2.0
    val rate = 4.0
    val rWSD = 1.0
    val (d1, d2) = makeDistributions(shape, rate, rWSD, EwensAttraction(q, mass, attraction),
      EwensAttraction.factory(q, Mass.factory(mass), Attraction.factory(distance, permutation, ReciprocalDecayFactory.factory(shape, rate, rdg))))
    assert(d1.approximatelyEquals(d2, epsilon))
  }

  "Equal probabilities subject to Monte Carlo error" should "hold for setting 2" in {
    val mass = Mass(0.8)
    val nItems = 4
    val distance = Distance.sample(nItems, rdg)
    val permutation = Permutation.sample(nItems, rdg)
    val attraction = Attraction(distance, permutation, ReciprocalDecayFactory(1.0))
    val shape = 3.0
    val rate = 1.0
    val rWSD = 0.5
    val (d1, d2) = makeDistributions(shape, rate, rWSD, EwensAttraction(q, mass, attraction),
      EwensAttraction.factory(q, Mass.factory(mass), Attraction.factory(distance, permutation, ReciprocalDecayFactory.factory(shape, rate, rdg))))
    assert(d1.approximatelyEquals(d2, epsilon))
  }

}

class UpdatePermutation extends FlatSpec {

  type Q = NullSamplingModel.tipe

  val nSamples = 100000
  val epsilon = 0.005

  val rdg = new RDG()

  def makeDistributions(k: Int, distribution: EwensAttraction[Q], factory: () => EwensAttraction[Q]) = {
    val db2 = DistributionBuilder[Q]()
    var p = Partition(Subset(q.sample(), 0, 1, 2, 3))
    var d = distribution
    for (rep <- 0 until nSamples) {
      p = AuxiliaryGibbsSampler(p, q, d, rdg)._1
      d = PermutationSampler.update(d, p, k, rdg, Set())._1
      db2.process(p, 1)
    }
    val sampler1 = PartitionModel.forwardSampler(p.nItems, factory)
    val db1 = DistributionBuilder[Q]()
    sampler1(nSamples, rdg).foreach(p => db1.process(p, 1))
    (db1.toDistribution, db2.toDistribution)
  }

  "Equal probabilities subject to Monte Carlo error" should "hold for setting 1" in {
    val mass = Mass(3.0)
    val nItems = 4
    val distance = Distance.sample(nItems, rdg)
    val permutation = Permutation.sample(nItems, rdg)
    val decay = ReciprocalDecayFactory(2.5)
    val attraction = Attraction(distance, permutation, decay)
    val k = 3
    val (d1, d2) = makeDistributions(k, EwensAttraction(q, mass, attraction),
      EwensAttraction.factory(q, Mass.factory(mass), Attraction.factory(distance, Permutation.factory(nItems, rdg), decay)))
    assert(d1.approximatelyEquals(d2, epsilon))
  }

}

class JainNeal extends FlatSpec {

  import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
  import org.ddahl.shallot.example._

  val rdg = new RDG()
  val data = Array(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
  val nReps = 20000
  val epsilon1 = 0.05
  val epsilon2 = 0.01

  val samplingModel1 = new NormalNormalModel(data, 0.1, 0.0, 1.0, rdg)
  var partition1 = Partition(() => samplingModel1.sample(), data.length, true)
  var sum1a = 0.0
  var sum1b = 0.0
  for (i <- 0 until nReps) {
    partition1 = AuxiliaryGibbsSampler(partition1, samplingModel1, Ewens(samplingModel1, Mass(1.0)), rdg)._1
    sum1a += partition1.nClusters
    sum1b += partition1.entropy
  }

  val samplingModel2 = new IntegratedNormalNormalModel(data, 0.1, 0.0, 1.0, rdg)
  var partition2 = Partition(() => samplingModel2.sample(), data.length, true)
  var sum2a = 0.0
  var sum2b = 0.0
  for (i <- 0 until nReps) {
    partition2 = AuxiliaryGibbsSampler(partition2, samplingModel2, Ewens(samplingModel2, Mass(1.0)), rdg)._1
    sum2a += partition2.nClusters
    sum2b += partition2.entropy
  }

  "Both conjugate and nonconjugate version of Neal (2000) model" should "given the same summaries" in {
    val mean1a = sum1a / nReps
    val mean1b = sum1b / nReps
    val mean2a = sum2a / nReps
    val mean2b = sum2b / nReps
    if (math.abs(mean1a - mean2a) > epsilon1) {
      println(mean1a + " " + mean2a)
      assert(false)
    }
    if (math.abs(mean1b - mean2b) > epsilon2) {
      println(mean1b + " " + mean2b)
      assert(false)
    }

  }

}

