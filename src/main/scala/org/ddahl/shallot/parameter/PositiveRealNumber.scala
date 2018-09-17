package org.ddahl.shallot
package parameter

// Immutable

class PositiveRealNumber protected (val value: Double) {

  if (value <= 0.0) throw new IllegalArgumentException("Value must be greater than 0.")
  override def toString = value.toString

}

