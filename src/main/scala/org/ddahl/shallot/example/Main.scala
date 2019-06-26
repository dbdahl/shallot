package org.ddahl.shallot
package example

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import parameter._
import parameter.decay._
import parameter.partition._
import distribution._
import mcmc._

object Main {

  val q = NullSamplingModel

  def test1() = {
    val rdg = new RDG()
    val permutationFactory = Permutation.factory(4, rdg)
    val decay = ReciprocalDecayFactory(2)
    val distance = Distance(Array(
      Array[Double](),
      Array[Double](1),
      Array[Double](2, 2),
      Array[Double](3, 1, 4)), false)
    println(distance)
    for (i <- 0 until 10) {
      val attraction = Attraction(distance, permutationFactory(), decay)
      println(attraction)
    }
    println(Subset(q.sample(), 1, 2, 5))
    println(Partition(Subset(q.sample(), 1, 2, 5), Subset(q.sample(), 3, 0, 3)).add(Subset(q.sample(), 4)).remove(Subset(q.sample(), 4)).remove(Subset(q.sample(), 5, 1, 2)))
    println(Permutation.enumerate(4).mkString("\n"))
    println(Partition.enumerate(() => q.sample(), 4).mkString("\n"))
    println(Partition.enumerate(() => q.sample(), 6).map(_.toLabels.mkString("")).mkString("\n"))
  }

  def test2() = {
    println(Ewens(q, Mass(0.1)).tabulate(4))
    println("---")
    println(EwensPitman(q, Mass(0.1), Discount(0.1)).tabulate(4))
  }

  def test3() = {
    val ewens = Ewens(q, Mass(0.8))
    val partition = Partition(Subset(q.sample(), 1, 2), Subset(q.sample(), 0))
    val rdg = new RDG()
    val distribution = ewens.fullConditional(partition.nItems, partition)
    println(distribution)
    val sampler = distribution.sampler(rdg)
    for (i <- 0 until 10) {
      println(sampler().toLabels.mkString(""))
    }
  }

  def test4() = {
    val rdg = new RDG()
    val distance = Distance(Array(
      Array[Double](),
      Array[Double](1),
      Array[Double](2, 2),
      Array[Double](3, 1, 4)), false)
    val massShape = 2.0
    val massRate = 2.0
    val massRWSD = 2.0
    val temperatureShape = 2.0
    val temperatureRate = 1.0
    val temperatureRWSD = 1.0
    val distribution1 = Ewens(q, Mass(0.8))
    val distribution2 = EwensPitman(q, Mass(0.8), Discount(0.1))
    val distribution3 = EwensAttraction(q, Mass(0.8), Attraction(distance, Permutation(1, 2, 0, 3), ReciprocalDecayFactory(1.0)))
    val sampler1 = distribution1.forwardSampler(4, rdg)
    val sampler2 = distribution2.forwardSampler(4, rdg)
    val sampler3 = PartitionModel.forwardSampler(4, EwensAttraction.factory(q, Mass.factory(massShape, massRate, rdg), Attraction.factory(distance, Permutation(1, 2, 0, 3), ReciprocalDecayFactory(1.0))))
    val viaMCMC = true
    //println(distribution1.tabulate(4))
    //println(distribution2.tabulate(4))
    //println(distribution3.tabulate(4))
    //println("----")
    val nSamples = 1000
    if (!viaMCMC) {
      println(sampler3(nSamples, rdg).map(_.toLabels.mkString("")).mkString("\n"))
    } else {
      var p = Partition(Subset(q.sample(), 0, 1, 2, 3))
      val k = 4
      var d = distribution3
      val monitorMass = AcceptanceRateMonitor()
      val monitorPermutation = AcceptanceRateMonitor()
      val monitorTemperature = AcceptanceRateMonitor()
      for (rep <- 0 until nSamples) {
        p = AuxiliaryGibbsSampler(p, q, d, rdg)._1
        d = monitorMass(MassSampler.gaussianRandomWalk(d, p, massShape, massRate, massRWSD, rdg))
        d = monitorPermutation(PermutationSampler.update(d, p, k, rdg, Set()))
        d = monitorTemperature(TemperatureSampler.gaussianRandomWalk(d, p, temperatureShape, temperatureRate, temperatureRWSD, rdg))
        println("" + d.mass + " " + d.attraction.decay + " " + monitorMass + " " + monitorPermutation + " " + monitorTemperature)
        //println(p.toLabels.mkString(""))
      }
    }
  }

  def test5() = {
    val rdg = new RDG()
    val mass = Mass(1.4)
    val distance = Distance(Array(
      Array[Double](),
      Array[Double](1),
      Array[Double](2, 2),
      Array[Double](3, 1, 4)), false)
    val distanceUniform = Distance(Array(
      Array[Double](),
      Array[Double](1),
      Array[Double](1, 1),
      Array[Double](1, 1, 1)), false)
    val ewensAttr = EwensAttraction(q, mass, Attraction(distance, Permutation(1, 2, 0, 3), ReciprocalDecayFactory(1.0)))
    val sampler = ewensAttr.forwardSampler(4, rdg)
    for (rep <- 0 until 100000) println(sampler().toLabels.mkString(""))
    val ewens = Ewens(q, mass)
    println("---")
    println(ewensAttr.tabulate(4))
  }

  def test6() = {
    val distance = Distance(Array(
      Array[Double](),
      Array[Double](1),
      Array[Double](2, 2),
      Array[Double](3, 1, 4)), false)
    val nItems = distance.nItems
    val rdg = new RDG()
    val massFactory = Mass.factory(1.2, 1.9, rdg)
    val permutationFactory = Permutation.factory(nItems, rdg)
    val decayFactory = ReciprocalDecayFactory.factory(2.5, 2.0, rdg)
    val attractionFactory = Attraction.factory(distance, permutationFactory, decayFactory)
    val distributionFactory = EwensAttraction.factory(q, massFactory, attractionFactory)
    val sampler = PartitionModel.forwardSampler(nItems, distributionFactory)
    val tabulator = PartitionModel.tabulator(nItems, distributionFactory)
    val samples = sampler(100, rdg)
    samples.foreach { x => println(x.toLabels.mkString("")) }
    val distributions = tabulator(100)
    distributions.foreach(println)
  }

  def main(args: Array[String]) = {
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
  }

}

