package org.ddahl.shallot
package r

import org.ddahl.rscala.{RClient, RPersistentReference}
import parameter._
import parameter.partition._

class RAdapter(engine: RClient, val usesPredictiveDensity: Boolean, logDensityFunction: RObject, sampleFunction: RObject) extends SamplingModel[RPersistentReference] {

  def sample = if (usesPredictiveDensity) null else engine.evalR(s"""$sampleFunction(c(),NULL)""")

  def sample(subset: Subset[RPersistentReference]) = if (usesPredictiveDensity) null
  else engine.evalR(s"""$sampleFunction(c(${subset.map(_ + 1).mkString(",")}),${subset.parameter})""")

  def copy(x: RPersistentReference) = x

  def logDensity(i: Int, subset: Subset[RPersistentReference]) = {
    val expr = if (usesPredictiveDensity) s"""${logDensityFunction}(${i + 1},c(${subset.map(_ + 1).mkString(",")}))"""
    else s"""${logDensityFunction}(${i + 1},${subset.parameter})"""
    val ld = engine.evalD0(expr)
    if ( java.lang.Double.isNaN(ld) || ( ld == Double.NegativeInfinity ) || ( ld == Double.PositiveInfinity ) ) {
      val sop = if ( usesPredictiveDensity ) "subset {"+subset.map(_ + 1).mkString(",")+"}"
      else "parameter: " + engine.evalS1(s"capture.output(print(${subset.parameter}))").mkString("\n") 
      throw new RuntimeException(s"logDensity for item ${i+1} given ${sop} returned ${ld}.")
    }
    ld
  }

}

object RAdapter {

  def apply(engine: RClient, usesPredictiveDensity: Boolean, logDensityFunction: RObject, sampleFunction: RObject) = new RAdapter(engine, usesPredictiveDensity, logDensityFunction, sampleFunction)

}

