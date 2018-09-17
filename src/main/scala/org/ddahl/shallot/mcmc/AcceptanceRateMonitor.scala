package org.ddahl.shallot.mcmc

class AcceptanceRateMonitor {

  private var _nAttempts = 0L
  private var _nAcceptances = 0L

  def nAcceptances = _nAcceptances
  def nAttempts = _nAttempts
  def rate = _nAcceptances / _nAttempts.toDouble

  def apply[B](x: (B, Boolean)): B = {
    _nAttempts += 1
    if (x._2) _nAcceptances += 1
    x._1
  }

  def reset() = {
    _nAttempts = 0
    _nAcceptances = 0
  }

  override def toString = rate.toString

}

object AcceptanceRateMonitor {

  def apply() = new AcceptanceRateMonitor()

}

