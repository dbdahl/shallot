package org.ddahl.shallot
package r

import org.ddahl.rscala.RClient
import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import parameter._
import parameter.decay._
import parameter.partition._
import distribution._
import mcmc._

abstract class MCMCSampler {

  val recordMass: Boolean
  val recordDiscount: Boolean
  val recordPermutation: Boolean
  val recordTemperature: Boolean

  val massMonitor = AcceptanceRateMonitor()
  val discountMonitor = AcceptanceRateMonitor()
  val permutationMonitor = AcceptanceRateMonitor()
  val temperatureMonitor = AcceptanceRateMonitor()

  def reset(): Unit = {
    massMonitor.reset()
    discountMonitor.reset()
    permutationMonitor.reset()
    temperatureMonitor.reset()
  }

  def rates: Array[Double] = {
    Array(massMonitor.rate, discountMonitor.rate, permutationMonitor.rate, temperatureMonitor.rate)
  }

  def sample(nIterations: Int, rdg: RDG): AugmentedSample

}

object RInterface {

  def makeMCMCSamplerForR(distribution: String, nItems: Int, distance: Distance,
    massTuple: (Double,Double,Double,Boolean,Double),
    discount: Double, discountShape1: Double, discountShape2: Double, discountFixed: Boolean, discountRWSD: Double,
    permutation: Permutation, permutationFixed: Boolean, permutationGrabSize: Int,
    decayString: String, decayMultiplier: Double,
    temperature: Double, temperatureShape: Double, temperatureRate: Double, temperatureFixed: Boolean, temperatureRWSD: Double,
    samplingModel: RAdapter): MCMCSampler = {
    val (mass, massShape, massRate, massFixed, massRWSD) = massTuple
    val (samplingModelP) = if (samplingModel == null) NullRSamplingModel
    else samplingModel
    var p = Partition(samplingModelP, nItems, true)
    val decay = decayString match {
      case "reciprocal" => ReciprocalDecayFactory(temperature)
      case "exponential" => ExponentialDecayFactory(temperature)
      case "subtraction" => (new SubtractionDecayFactory(decayMultiplier * distance.max)) (temperature)
      case "" => null
    }
    distribution match {
      case "shallot.distribution.ewens" =>
        var priorModel = Ewens(samplingModelP, Mass(mass))
        new MCMCSampler() {
          val recordMass = true
          val recordDiscount = false
          val recordPermutation = false
          val recordTemperature = false
          def sample(nIterations: Int, rdg: RDG) = {
            for (i <- 0 until nIterations) {
              if (!massFixed) priorModel = massMonitor(MassSampler.escobarWest(priorModel, p, massShape, massRate, rdg))
              p = AuxiliaryGibbsSampler(p, samplingModelP, priorModel, rdg)._1
            }
            (p, priorModel.mass.value, discount, permutation, temperature)
          }
        }
      case "shallot.distribution.ewensPitman" =>
        var priorModel = EwensPitman(samplingModelP, Mass(mass), Discount(discount))
        new MCMCSampler() {
          val recordMass = true
          val recordDiscount = true
          val recordPermutation = false
          val recordTemperature = false
          def sample(nIterations: Int, rdg: RDG) = {
            for (i <- 0 until nIterations) {
              if (!massFixed) priorModel = massMonitor(MassSampler.gaussianRandomWalk(priorModel, p, massShape, massRate, massRWSD, rdg))
              if (!discountFixed) priorModel = discountMonitor(DiscountSampler.gaussianRandomWalk(priorModel, p, discountShape1, discountShape2, discountRWSD, rdg))
              p = AuxiliaryGibbsSampler(p, samplingModelP, priorModel, rdg)._1
            }
            (p, priorModel.mass.value, priorModel.discount.value, permutation, temperature)
          }
        }
      case "shallot.distribution.ewensAttraction" =>
        var priorModel = EwensAttraction(samplingModelP, Mass(mass), Attraction(distance, permutation, decay))
        new MCMCSampler() {
          val recordMass = true
          val recordDiscount = false
          val recordPermutation = true
          val recordTemperature = true
          def sample(nIterations: Int, rdg: RDG) = {
            for (i <- 0 until nIterations) {
              if (!massFixed) priorModel = massMonitor(MassSampler.escobarWest(priorModel, p, massShape, massRate, rdg))
              if (!permutationFixed) priorModel = permutationMonitor(PermutationSampler.update(priorModel, p, permutationGrabSize, rdg, Set()))
              if (!temperatureFixed) priorModel = temperatureMonitor(TemperatureSampler.gaussianRandomWalk(priorModel, p, temperatureShape, temperatureRate, temperatureRWSD, rdg))
              p = AuxiliaryGibbsSampler(p, samplingModelP, priorModel, rdg)._1
            }
            (p, priorModel.mass.value, discount, priorModel.attraction.permutation, priorModel.attraction.decay.temperature)
          }
        }
      case "shallot.distribution.ewensPitmanAttraction" =>
        var priorModel = EwensPitmanAttraction(samplingModelP, Mass(mass), Discount(discount), Attraction(distance, permutation, decay))
        new MCMCSampler() {
          val recordMass = true
          val recordDiscount = true
          val recordPermutation = true
          val recordTemperature = true
          def sample(nIterations: Int, rdg: RDG) = {
            for (i <- 0 until nIterations) {
              if (!massFixed) priorModel = massMonitor(MassSampler.gaussianRandomWalk(priorModel, p, massShape, massRate, massRWSD, rdg))
              if (!discountFixed) priorModel = discountMonitor(DiscountSampler.gaussianRandomWalk(priorModel, p, discountShape1, discountShape2, discountRWSD, rdg))
              if (!permutationFixed) priorModel = permutationMonitor(PermutationSampler.update(priorModel, p, permutationGrabSize, rdg, Set()))
              if (!temperatureFixed) priorModel = temperatureMonitor(TemperatureSampler.gaussianRandomWalk(priorModel, p, temperatureShape, temperatureRate, temperatureRWSD, rdg))
              p = AuxiliaryGibbsSampler(p, samplingModelP, priorModel, rdg)._1
            }
            (p, priorModel.mass.value, priorModel.discount.value, priorModel.attraction.permutation, priorModel.attraction.decay.temperature)
          }
        }
      case _ => throw new IllegalArgumentException("Unsupported distribution.")
    }
  }

  def sampleMCMC(nSamples: Int, nIterationsPerSample: Int, rdg: RDG, sampler: MCMCSampler): AugmentedSamples = {
    var list: AugmentedSamples = List[AugmentedSample]()
    for (i <- 0 until nSamples) {
      list = sampler.sample(nIterationsPerSample, rdg) +: list
    }
    list
  }

  def extractNSubsets(list: Samples): Array[Int] = list.map(_.nSubsets).toArray
  def extractEntropy(list: Samples): Array[Double] = list.map(_.entropy).toArray
  def extractSamples(list: AugmentedSamples): Samples = list.map(_._1)
  def extractMasses(list: AugmentedSamples): Array[Double] = list.map(_._2).toArray
  def extractDiscounts(list: AugmentedSamples): Array[Double] = list.map(_._3).toArray
  def extractPermutations(list: AugmentedSamples): List[Permutation] = list.map(_._4)
  def extractTemperatures(list: AugmentedSamples): Array[Double] = list.map(_._5).toArray

  def makeForwardSamplerForR(distribution: String, nItems: Int, distance: Distance,
    mass: Double, massShape: Double, massRate: Double, massFixed: Boolean,
    discount: Double, discountShape1: Double, discountShape2: Double, discountFixed: Boolean,
    permutation: Permutation, permutationFixed: Boolean,
    decayString: String, decayMultiplier: Double,
    temperature: Double, temperatureShape: Double, temperatureRate: Double, temperatureFixed: Boolean,
    ppmxSimilarity: PPMxSimilarity) = (nSamples: Int, rdg: RDG) => {
    val massFactory = if (massFixed) Mass.factory(mass)
    else Mass.factory(massShape, massRate, rdg)
    val discountFactory = if (discountFixed) Discount.factory(discount)
    else Discount.factory(discountShape1, discountShape2, rdg)
    var list = List[Partition[NullRSamplingModel.tipe]]()
    val distributionFactory = distribution match {
      case "shallot.distribution.ewens" =>
        Ewens.factory(NullRSamplingModel, massFactory)
      case "shallot.distribution.ewensPitman" =>
        EwensPitman.factory(NullRSamplingModel, massFactory, discountFactory)
      case "shallot.distribution.ppmx" =>
        PPMx.factory(NullRSamplingModel, massFactory, ppmxSimilarity)
      case "shallot.distribution.ewensAttraction" | "shallot.distribution.ewensPitmanAttraction" | "shallot.distribution.ddcrp" =>
        val attractionFactory = if (temperatureFixed) {
          val decay = decayString match {
            case "reciprocal" => ReciprocalDecayFactory(temperature)
            case "exponential" => ExponentialDecayFactory(temperature)
            case "subtraction" => (new SubtractionDecayFactory(decayMultiplier * distance.max)) (temperature)
          }
          if (permutationFixed) Attraction.factory(distance, permutation, decay)
          else Attraction.factory(distance, Permutation.factory(nItems, rdg), decay)
        } else {
          val decayFactory = decayString match {
            case "reciprocal" => ReciprocalDecayFactory.factory(temperatureShape, temperatureRate, rdg)
            case "exponential" => ExponentialDecayFactory.factory(temperatureShape, temperatureRate, rdg)
            case "subtraction" => (new SubtractionDecayFactory(decayMultiplier * distance.max)).factory(temperatureShape, temperatureRate, rdg)
          }
          if (permutation.nItems != 0) Attraction.factory(distance, permutation, decayFactory)
          else Attraction.factory(distance, Permutation.factory(nItems, rdg), decayFactory)
        }
        if (distribution == "shallot.distribution.ewensAttraction") EwensAttraction.factory(NullRSamplingModel, massFactory, attractionFactory)
        else if (distribution == "shallot.distribution.ewensPitmanAttraction") EwensPitmanAttraction.factory(NullRSamplingModel, massFactory, discountFactory, attractionFactory)
        else if (distribution == "shallot.distribution.ddcrp") {
          val factory = DDCRP.factory(NullRSamplingModel, massFactory, attractionFactory)
          list = List.fill(nSamples) { factory().sample(rdg) }
          EwensAttraction.factory(NullRSamplingModel, massFactory, attractionFactory) // Dummy just to make the return type correct
        } else throw new IllegalArgumentException("Unsupported distribution.")
      case _ => throw new IllegalArgumentException("Unsupported distribution.")
    }
    if (list.isEmpty) {
      val sampler = PartitionModel.forwardSampler(nItems, distributionFactory)
      sampler(nSamples, rdg)
    } else list
  }

  def sampleForward(nSamples: Int, rdg: RDG, sampler: Function2[Int, RDG, Samples], parallel: Boolean): Samples = {
    if (!parallel) sampler(nSamples, rdg)
    else {
      val nCores = Runtime.getRuntime.availableProcessors
      val nSamplesPerCore = (nSamples / nCores) + 1
      val randomGenerator = rdg.getRandomGenerator
      val rdgList = List.fill(nCores) { new RDG(randomGenerator) }
      rdgList.par.map(r => sampler(nSamplesPerCore, r)).toList.flatten
    }
  }

  def mergeSamples(list1: Samples, list2: Samples): Samples = list1 ++ list2

  def mergeAugmentedSamples(list1: AugmentedSamples, list2: AugmentedSamples): AugmentedSamples = list1 ++ list2

  def listToArray(a: List[Double]) = a.toArray

  def matrixToPartitions(partitions: Array[Array[Int]]): Samples = {
    partitions.map(x => Partition(NullRSamplingModel, x)).toList
  }

  def pad(i: Int, length: Int) = {
    if (length == 1) i.toString
    else {
      val width = math.floor(math.log10(length - 1)).toInt + 1
      ("%0" + width + "d").format(i)
    }
  }

  def toLabelsWithParameters(partition: Partition[RPersistentReference]): (Array[Int], Array[RPersistentReference]) = {
    val result1 = new Array[Int](partition.nItems)
    val result2 = new Array[RPersistentReference](partition.nSubsets)
    var label = 0
    partition.toList.sortWith(_.min < _.min).foreach(subset => {
      result2(label) = subset.parameter
      label += 1
      subset.foreach(i => result1(i) = label)
    })
    (result1, result2)
  }

  def partition(x: Array[Int]): Partition[RPersistentReference] = Partition(NullRSamplingModel, x)

  def rotateForConfidencePlot(pp: PairwiseProbability, order: Array[Int]): Array[Array[Double]] = {
    val nItems = pp.nItems
    val xx = Array.ofDim[Double](nItems, nItems)
    for (i <- 0 until nItems) {
      for (j <- 0 until nItems) {
        xx(i)(nItems - j - 1) = pp(order(i) - 1, order(j) - 1)
      }
    }
    xx
  }

}

