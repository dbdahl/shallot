package org.ddahl.shallot
package parameter
package decay

import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
import org.apache.commons.math3.util.FastMath.{ pow, log, exp }

trait DecayFunction {

  def apply(x: Double): Double
  val temperature: Double
  def update(temperature: Double): DecayFunction

  if (temperature < 0.0) throw new IllegalArgumentException("Temperature must be nonnegative.")

}

abstract class DecayFunctionFactory[+B] {

  def apply(temperature: Double): B

  def factory(temperature: Double): () => B = {
    val decay = apply(temperature)
    () => decay
  }

  def factory(shape: Double, rate: Double, rdg: RDG): () => B = {
    val scale = 1.0 / rate
    () => apply(rdg.nextGamma(shape, scale))
  }

}

class ReciprocalDecay(val temperature: Double) extends DecayFunction {

  def apply(x: Double) = 1.0 / pow(x, temperature)
  def update(temp: Double) = ReciprocalDecayFactory(temp)

}

object ReciprocalDecayFactory extends DecayFunctionFactory[ReciprocalDecay] {

  def apply(temperature: Double) = new ReciprocalDecay(temperature)

}

class ExponentialDecay(val temperature: Double) extends DecayFunction {

  def apply(x: Double) = exp(-temperature * x)
  def update(temp: Double) = ExponentialDecayFactory(temp)

}

object ExponentialDecayFactory extends DecayFunctionFactory[ExponentialDecay] {

  def apply(temperature: Double) = new ExponentialDecay(temperature)

}

class SubtractionDecay(val temperature: Double, val maxValue: Double) extends DecayFunction {

  def apply(x: Double) = pow(maxValue - x, temperature)

  def update(temp: Double) = (new SubtractionDecayFactory(maxValue)) (temp)

}

class SubtractionDecayFactory(val maxValue: Double) extends DecayFunctionFactory[SubtractionDecay] {

  def apply(temperature: Double) = new SubtractionDecay(temperature, maxValue)

}

object SubtractionDecayFactory {

  def apply(maxValue: Double) = new SubtractionDecayFactory(maxValue)

}
