package org.ddahl.shallot
package example

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import parameter._
import partition._
import r._
import distribution._
import mcmc._

object Neal2000InR {

  def doit[A](nItems: Int, model: SamplingModel[A]) = {
    val rdg = new RDG()
    val ewens = Ewens(model, Mass(1.0))
    var p = Partition(model, nItems, true)
    val nReps = 2000
    for (i <- 0 until nReps) {
      p = AuxiliaryGibbsSampler(p, model, ewens, rdg)._1
      println(p.nSubsets + " " + p.entropy + " # " + p)
    }
  }

  def main(model: RAdapter): Unit = {
    doit(10, model)
  }

}

