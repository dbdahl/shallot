package org.ddahl.shallot
package example

import parameter._
import parameter.decay._
import parameter.partition._
import distribution._
import mcmc._
import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }

object Neal2000 {

  def main(args: Array[String]) = {
    val rdg = new RDG()
    val data = Array(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
    val nItems = data.length
    val samplingModel = new NormalNormalModel(data, 0.1, 0.0, 1.0, rdg)
    val massShape = 2.0
    val massRate = 4.0
    val massRWSD = 1.0
    val nReps = 10000
    var partition = Partition(() => samplingModel.sample(), data.length, true)
    val what = if (args.length > 0) args(0)
    else "ea"
    what match {
      case "neal2000" =>
        var priorModel = Ewens(samplingModel, Mass(1.0))
        for (i <- 0 until nReps) {
          partition = AuxiliaryGibbsSampler(partition, samplingModel, priorModel, rdg)._1
          println(partition.nClusters + " " + partition.entropy + " # " + partition)
        }
      case "neal2000integrated" =>
        val samplingModel2 = new IntegratedNormalNormalModel(data, 0.1, 0.0, 1.0, rdg)
        var partition2 = Partition(() => samplingModel2.sample(), data.length, true)
        var priorModel2 = Ewens(samplingModel2, Mass(1.0))
        for (i <- 0 until nReps) {
          partition2 = AuxiliaryGibbsSampler(partition2, samplingModel2, priorModel2, rdg)._1
          println(partition2.nClusters + " " + partition2.entropy + " # " + partition2)
        }
      case "ea" => {
        var priorModel = EwensAttraction(samplingModel, Mass(1.0), Attraction(Distance.sample(nItems, rdg), Permutation.sample(nItems, rdg), ReciprocalDecayFactory(1.0)))
        val temperatureShape = 2.0
        val temperatureRate = 4.0
        val temperatureRWSD = 1.0
        val grabSize = 4
        val monitor = AcceptanceRateMonitor()
        for (i <- 0 until nReps) {
          partition = AuxiliaryGibbsSampler(partition, samplingModel, priorModel, rdg)._1
          priorModel = MassSampler.escobarWest(priorModel, partition, massShape, massRate, rdg)._1
          priorModel = TemperatureSampler.gaussianRandomWalk(priorModel, partition, temperatureShape, temperatureRate, temperatureRWSD, rdg)._1
          priorModel = monitor(PermutationSampler.update(priorModel, partition, grabSize, rdg, Set()))
          println(partition.nClusters + " " + partition.entropy + " " + priorModel.mass + " " + priorModel.attraction.decay + " " + monitor + " # " + partition)
        }
      }
    }
  }

}

