package org.ddahl.shallot
package example

import distribution.Attraction
import distribution.PartitionModel
import distribution.Ewens
import distribution.EwensPitman
import distribution.EwensPitmanAttraction
import parameter.NullSamplingModel
import parameter.Mass
import parameter.Discount
import parameter.Permutation
import parameter.Similarity
import parameter.partition.Partition
import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }

object CRA2CRP {

  def totalVariationDistance(p: Double, q: Double) = 0.5 * math.abs(p - q)
  def kullbackLeibler(p: Double, q: Double) = p * (if (p == q) 0.0 else math.log(p / q))
  def jensenShannon(p: Double, q: Double) = 0.5 * (kullbackLeibler(p, q) + kullbackLeibler(q, p))

  def compute(partitions: scala.collection.parallel.immutable.ParSeq[Partition[NullSamplingModel.tipe]], d1: PartitionModel[NullSamplingModel.tipe], d2: PartitionModel[NullSamplingModel.tipe], fn: (Double, Double) => Double) = {
    partitions.map(partition => fn(d1.probability(partition), d2.probability(partition))).sum
  }

  /*
  def main9(args: Array[String]): Unit = {
    val nSamples = 1000000
    val nItems = 4
    val mass = Mass.factory(6.0)
    val discount = Discount.factory(0.3)

    val c = EwensPitman.meanNumberOfSubsets(nItems, mass(), discount())
    println(c)
    val aa = EwensPitman.sampleNumberOfSubsets(nItems, mass, discount, nSamples)
    val mean = aa.sum.toDouble / nSamples
    println(mean)
    return

    val start1 = System.nanoTime()
    //val a = Ewens.sampleNumberOfSubsets(nItems, mass, nSamples).sum.toDouble / nSamples
    val a = EwensPitman.sampleNumberOfSubsets(nItems, mass, discount, nSamples).sum.toDouble / nSamples
    val stop1 = System.nanoTime()
    val d1 = EwensPitman.factory(NullSamplingModel, mass, discount)
    //val d1 = Ewens.factory(NullSamplingModel, mass)
    val sampler = PartitionModel.forwardSampler(nItems, d1)
    val rdg = new RDG()
    val start2 = System.nanoTime()
    val b = sampler(nSamples, rdg).map(_.nSubsets).sum.toDouble / nSamples
    val stop2 = System.nanoTime()
    println(a + " " + (stop1 - start1))
    println(b + " " + (stop2 - start2))
    println((stop2 - start2) / (stop1 - start1))
  }

  def main66(args: Array[String]): Unit = {
    val mass = Mass(1.0)
    val discount = Discount(0.05)
    val nItems = 4
    val partitions = Partition.enumerate(NullSamplingModel, nItems).filter(p => p.paired(1, 3))
    val d1 = EwensPitmanAttraction(NullSamplingModel, mass, discount, Attraction.constant(nItems))
    partitions.foreach(p => {
      val x = d1.probability(p)
      println()
      //println(p.toLabels.mkString("") + "" + d1.probability(p))
    })
  }

  def main10(args: Array[String]): Unit = {
    val mass = Mass(0.01)
    val discount = Discount(0.99)
    val nItems = 6
    val partitions = Partition.enumerate(NullSamplingModel, nItems)
    val d1 = EwensPitman(NullSamplingModel, mass, discount)
    for (i <- 1 to nItems) {
      println("" + EwensPitman.probabilityNumberOfSubsets(nItems, i, mass, discount) + "  " + partitions.filter(_.nSubsets == i).map(d1.probability).sum)
    }
    val a = 20
    for (i <- 1 to a) {
      println(EwensPitman.probabilityNumberOfSubsets(a, i, mass, discount))
    }
  }

  def main2(args: Array[String]): Unit = {
    for (nItems <- Seq(4, 5, 6, 7, 8, 9, 10)) {
      val partitions = Partition.enumerate(NullSamplingModel, nItems).par
      for (mass <- Seq(0.01, 0.5, 1.0, 2.0, 5.0).map(Mass(_))) {
        for (discount <- Seq(0.1, 0.2, 0.3, 0.4, 0.99).map(Discount(_))) {
          val d0 = Ewens(NullSamplingModel, mass)
          val d1 = EwensPitman(NullSamplingModel, mass, discount)
          val d2 = EwensPitmanAttraction(NullSamplingModel, mass, discount, Attraction.constant(nItems))
          val kl1 = compute(partitions, d0, d1, kullbackLeibler)
          val kl2 = compute(partitions, d0, d2, kullbackLeibler)
          if (kl1 < kl2) {
            println("nItems = " + nItems)
            println("Mass = " + mass)
            println("Discount = " + discount)
            println(kl1)
            println(kl2)
            println("---")
          }
        }
      }
    }
    println("Done")
  }

  def main3(args: Array[String]): Unit = {
    val nItems = 3
    val partitions = Partition.enumerate(NullSamplingModel, nItems).par
    val mass = Mass(1.0)
    val discount = Discount(0.1)
    val d = EwensPitmanAttraction(NullSamplingModel, mass, discount, Attraction.constant(nItems))
    for (p <- partitions) {
      println(p + ": " + d.probability(p))
    }
  }

  def main4(args: Array[String]): Unit = {
    val nItems = 4
    val mass = Mass(1.0)
    val discount = Discount(0.1)
    val similarity1 = Similarity(Array(
      Array(1, 2, 3),
      Array(2, 1, 6),
      Array(3, 6, 1)))
    val similarity2 = Similarity(Array(
      Array(1, 2, 3, 900),
      Array(2, 1, 6, 8),
      Array(3, 6, 1, 20),
      Array(900, 8, 20, 1)))
    val attraction1 = Attraction(similarity1, Permutation.natural(similarity1.nItems), Temperature(1.0))
    val attraction2 = Attraction(similarity2, Permutation.natural(similarity2.nItems), Temperature(1.0))
    val d1 = EwensPitmanAttraction(NullSamplingModel, mass, discount, attraction1)
    val d2 = EwensPitmanAttraction(NullSamplingModel, mass, discount, attraction2)
    for (p <- Partition.enumerate(NullSamplingModel, similarity1.nItems)) {
      println(p.toLabels.mkString("") + ": %4.4f".format(d1.probability(p)))
    }
    println("---")
    for (p <- Partition.enumerate(NullSamplingModel, similarity2.nItems)) {
      println(p.toLabels.mkString("") + ": %4.4f".format(d2.probability(p)))
    }
  }

  def main5(args: Array[String]): Unit = {
    val nItems = 4
    val mass = Mass(2.0)
    val discount = Discount(0.1)
    val similarity = Similarity(Array(
      Array(1, 1, 1, 1),
      Array(1, 1, 1, 1),
      Array(1, 1, 1, 1),
      Array(1, 1, 1, 1)))
    val permutations = Permutation.enumerate(nItems)
    for (p <- Partition.enumerate(NullSamplingModel, similarity.nItems)) {
      val label = p.toLabels.map(_ + 1).mkString("")
      println("(* " + label + " *)")
      println("p" + label + " = 0")
      permutations.foreach(perm => {
        print("p" + label + " = p" + label + " + ")
        val attraction = Attraction(similarity, perm, Temperature(1.0))
        val d = EwensPitmanAttraction(NullSamplingModel, mass, discount, attraction)
        val prob = d.probability(p)
        println("""""")
      })
      println
    }
  }

  def main7(args: Array[String]): Unit = {
    val nItems = 4
    val mass = Mass(2.0)
    val discount = Discount(0.1)
    val similarity = Similarity(Array(
      Array(1, 1, 1, 1),
      Array(1, 1, 1, 1),
      Array(1, 1, 1, 1),
      Array(1, 1, 1, 1)))
    val permutations = Permutation.enumerate(nItems)
    for (p <- Partition.enumerate(NullSamplingModel, similarity.nItems)) {
      val label = p.toLabels.map(_ + 1).mkString("")
      println("(* " + label + " *)")
      println("p" + label + " = 0")
      permutations.foreach(perm => {
        print("p" + label + " = p" + label + " + ")
        val attraction = Attraction(similarity, perm, Temperature(1.0))
        val d = EwensPitmanAttraction(NullSamplingModel, mass, discount, attraction)
        val prob = d.probability(p)
        println("""""")
      })
      println
    }
  }
  */

}
