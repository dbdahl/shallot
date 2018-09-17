package org.ddahl.shallot
package distribution

import parameter._
import parameter.partition._

trait HasSamplingModel[+B, A] {

  val samplingModel: SamplingModel[A]
  def self: B
  def replaceSamplingModel(newSamplingModel: SamplingModel[A]): B

}

trait HasMass[+B, A] {

  val mass: Mass
  def self: B
  def replaceMass(newMass: Mass): B

}

trait HasMassEscobarWest[+B, A] extends HasMass[B, A]

trait HasDiscount[+B, A] {

  val discount: Discount
  def self: B
  def replaceDiscount(newDiscount: Discount): B

}

trait HasAttraction[+B, A] {

  val attraction: Attraction
  def self: B
  def replaceAttraction(newAttraction: Attraction): B

}

