package org.ddahl.shallot
package parameter

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.special.Gamma.logGamma
import org.apache.commons.math3.util.FastMath.log
import parameter._
import parameter.partition._

// Immutable

class Mass(value: Double) extends PositiveRealNumber(value) {

  val logValue = log(value)
  val logGammaValue = logGamma(value)

}

object Mass {

  def apply(value: Double) = new Mass(value)

  def factory(mass: Mass) = {
    () => mass
  }

  def factory(value: Double) = {
    val mass = Mass(value)
    () => mass
  }

  def factory(shape: Double, rate: Double, rdg: RDG) = {
    val scale = 1.0 / rate
    () => Mass(rdg.nextGamma(shape, scale))
  }

}

