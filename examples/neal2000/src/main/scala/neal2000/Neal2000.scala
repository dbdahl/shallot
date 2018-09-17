package neal2000

import org.ddahl.shallot.parameter.Mass
import org.ddahl.shallot.parameter.partition.Partition
import org.ddahl.shallot.distribution.Ewens
import org.ddahl.shallot.mcmc.AuxiliaryGibbsSampler
import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }

object Main {

  def main(args: Array[String]) = {
    val rdg = new RDG()
    val data = Array(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
    val samplingModel = new NormalNormalModel(data, 0.1, 0.0, 1.0, rdg)
    val nReps = 10000
    val priorModel = Ewens(samplingModel, Mass(1.0))
    var partition = Partition(samplingModel, data.length, true)
    for (i <- 0 until nReps) {
      partition = AuxiliaryGibbsSampler(partition, samplingModel, priorModel, rdg)._1
      println(partition.nSubsets + " " + partition.entropy + " # " + partition)
    }
  }

}

