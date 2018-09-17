package org.ddahl.shallot
package mcmc

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.exp
import parameter.partition._

object PosteriorPrediction {

  def apply[A](partition: Partition[A], logPosteriorPredictive: (Partition[A]) => Iterable[(Subset[A],Double)], rdg: RDG): Subset[A] = {
    val x0 = logPosteriorPredictive(partition)
    val x1 = x0.map(y => (y._1, exp(y._2)))
    val sum = x1.foldLeft(0.0)((sum, y) => sum + y._2)
    val x2 = if (sum == 0.0) {
      val values = x1.map(y => y._2).toArray
      val max = values.max
      val index = values.indexWhere(_ == max)
      val w = new Array[Double](values.length)
      w(index) = 1.0
      x1.map(y => y._1).zip(w)
    } else {
      x1.map(y => (y._1, y._2 / sum))
    }
    val unif = rdg.nextUniform(0.0, 1.0, true)
    var cumsum = 0.0
    val selectedSubset = x2.find(element => {
      cumsum += element._2
      cumsum > unif
    }).get._1
    selectedSubset
  }

}
