package org.ddahl.shallot
package distribution

import org.apache.commons.math3.special.Gamma.logGamma
import org.apache.commons.math3.util.FastMath.log
import parameter._
import parameter.partition._

class EwensPitmanAttraction[A] private (val samplingModel: SamplingModel[A], val mass: Mass, val discount: Discount, val attraction: Attraction) extends PartitionModel[A] with HasSamplingModel[EwensPitmanAttraction[A], A] with HasMass[EwensPitmanAttraction[A], A] with HasDiscount[EwensPitmanAttraction[A], A] with HasAttraction[EwensPitmanAttraction[A], A] {

  def self = this

  def replaceSamplingModel(newSamplingModel: SamplingModel[A]) = new EwensPitmanAttraction(newSamplingModel, mass, discount, attraction)

  def replaceMass(newMass: Mass) = new EwensPitmanAttraction(samplingModel, newMass, discount, attraction)

  def replaceDiscount(newDiscount: Discount) = new EwensPitmanAttraction(samplingModel, mass, newDiscount, attraction)

  def replaceAttraction(newAttraction: Attraction) = new EwensPitmanAttraction(samplingModel, mass, discount, newAttraction)

  override def iteratorForForwardSampling(nItems: Int) = {
    if (nItems != attraction.nItems) throw new IllegalArgumentException("Number of items is inconsistent.")
    attraction.permutation.iterator
  }

  override def logFullConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    val p2 = partition.add(i, subset)
    var p = Partition.empty[A]()
    var sum = 0.0
    val perm = attraction.permutation
    iteratorForForwardSampling(p2.nItems).foreach(k => {
      val sOption = p.find(s => p2.paired(k, s.head))
      val s = if (sOption.isEmpty) Subset.empty(samplingModel.sample()) else sOption.get
      if (perm.nPreceeding(k) >= perm.nPreceeding(i)) sum += logCascadingConditional(k, s, p)
      p = p.add(k, s)
    })
    sum
  }

  def logCascadingConditional(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    val numerator = if (subset.size == 0) {
      mass.value + partition.nClusters * discount.value
    } else {
      attraction(i, subset) * (1.0 - partition.nClusters * discount.value / partition.nItems)
    }
    log(numerator / (mass.value + partition.nItems))
  }
  
  def logPredictive(distanceToHeldOut: Array[Double], partition: Partition[A]): Iterable[(Subset[A],Double)] = {
    val s = distanceToHeldOut.map(attraction.decay(_))
    if (!s.forall(_ > 0.0)) throw new IllegalArgumentException("Not all the elements are strictly positive.")
    if (!s.forall(_ < Double.PositiveInfinity)) throw new IllegalArgumentException("Not all the elements are strictly less than infinity.")
    val p = partition ++ Iterable(Subset.empty(samplingModel.sample()))
    p.map(subset => {
      val numerator = if (subset.size == 0) {
        mass.value + partition.nClusters * discount.value
      } else {
        subset.foldLeft(0.0)((sum, k) => sum + s(k)) * (1.0 - partition.nClusters * discount.value / partition.nItems)
      }
      (subset, log(numerator / (mass.value + partition.nItems)))
    })
  }

  def logCascadingConditionalForSymbolicExpressions(i: Int, subset: Subset[A], partition: Partition[A]): Double = {
    //val sMass = """\[Alpha]"""
    val sMass = """1"""
    //val sDiscount = """\[Delta]"""
    val sDiscount = """0"""
    def similarity(ii: Int, jj: Int) = {
      val (iii, jjj) = if (ii < jj) (ii + 1, jj + 1) else (jj + 1, ii + 1)
      """Subscript[\[Lambda],""" + iii + "," + jjj + "]"
    }
    val numerator = if (subset.size == 0) {
      print("""(%s + %s %s)""".format(sMass, partition.nClusters, sDiscount))
      mass.value + partition.nClusters * discount.value
    } else {
      print("""(%s - %s %s)""".format(partition.nItems, partition.nClusters, sDiscount))
      if (subset.toList.sortWith(_ < _) != partition.flatten.toList.sortWith(_ < _)) {
        print("""((""")
        print(subset.toList.sortWith(_ < _).map(j => similarity(i, j)).mkString(" + "))
        print(""")/(""")
        print(partition.flatten.toList.sortWith(_ < _).map(j => similarity(i, j)).mkString(" + "))
        print("""))""")
      }
      attraction(i, subset) * (1.0 - partition.nClusters * discount.value / partition.nItems)
    }
    log(numerator / (mass.value + partition.nItems))
  }

}

object EwensPitmanAttraction extends TwoParameterNumberOfSubset {

  def apply[A](samplingModel: SamplingModel[A], mass: Mass, discount: Discount, attraction: Attraction) = new EwensPitmanAttraction(samplingModel, mass, discount, attraction)

  def factory[A](samplingModel: SamplingModel[A], massFactory: () => Mass, discountFactory: () => Discount, attractionFactory: () => Attraction) = () => apply(samplingModel, massFactory(), discountFactory(), attractionFactory())

}

