package org.ddahl.shallot.parameter

import partition._

trait SamplingModel[A] {

  type tipe = A

  def sample(): A

  def sample(subset: Subset[A]): A

  def copy(x: A): A = x

  def logDensity(i: Int, subset: Subset[A]): Double

}

class GeneralNullSamplingModel[A] extends SamplingModel[A] {

  def sample(): A = null.asInstanceOf[A]

  def sample(subset: Subset[A]) = null.asInstanceOf[A]

  def logDensity(i: Int, subset: Subset[A]): Double = 0.0

}

object NullSamplingModel extends GeneralNullSamplingModel[Null]

