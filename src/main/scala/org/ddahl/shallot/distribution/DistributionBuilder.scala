package org.ddahl.shallot
package distribution

import parameter._
import parameter.partition._

class DistributionBuilder[A] {

  private var map = new scala.collection.mutable.HashMap[Partition[A], Double]() {
    override def default(key: Partition[A]) = 0.0
  }

  def process(partition: Partition[A], increment: Double) = {
    map(partition) += increment
  }

  def toDistribution = Distribution(map, false, true)

}

object DistributionBuilder {

  def apply[A]() = new DistributionBuilder[A]()

}

