package org.ddahl.shallot
package parameter

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }

// Immutable

class Discount(value: Double) extends UnitIntervalRealNumber(value)

object Discount {

  def apply(value: Double) = new Discount(value)

  def factory(discount: Discount) = {
    () => discount
  }

  def factory(value: Double) = {
    val discount = Discount(value)
    () => discount
  }

  def factory(shape1: Double, shape2: Double, rdg: RDG) = {
    () => Discount(rdg.nextBeta(shape1, shape2))
  }

}

