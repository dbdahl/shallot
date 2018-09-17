package org.ddahl.shallot
package parameter

// Immutable

class UnitIntervalRealNumber protected (val value: Double) {

  if ((value < 0.0) || (value >= 1.0)) throw new IllegalArgumentException("Value must be in [0,1).")
  override def toString = value.toString

}

