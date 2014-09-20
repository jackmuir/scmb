/*
The MIT License (MIT)

Copyright (c) 2014 Jack Broderick Muir

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
package scmb.cmbhmc

import breeze.linalg.{DenseVector, sum}

import scmb.shdatasets._

import scmb.hmc._

import math.{log, Pi}

class cmbHMC(ll: Int, ee: Double, maxl: Int, dataSets: List[Map[String, String]]) extends HMC(ll, ee) {
  val cmbSets = for (dataSet <- dataSets) yield {
    dataSet("type") match {
      case "residuals" => new ResTT(maxl, dataSet)
      case "PcP-P" => new PcPmP(maxl, dataSet)
      case "P4KP-PcP" => new P4KPmPcP(maxl, dataSet)
    }
  }

  def calcU(q: Position, sigs: HierarchicalList): Double = {
    val pSums = for (cmbSet <- cmbSets) yield cmbSet.pCalcU(cmbSet.gMatrix, q, cmbSet.residuals, sigs(cmbSets.indexOf(cmbSet)))
    pSums.reduceLeft(_ + _)
    }

  def calcGradU(q: Position, sigs: HierarchicalList): DenseVector[Double] = {
    val pGrads = for (cmbSet <- cmbSets) yield cmbSet.pCalcGradU(cmbSet.gMatrix, q, cmbSet.residuals, sigs(cmbSets.indexOf(cmbSet)))
    pGrads.reduceLeft(_ + _)
  }

  def mLogLike(q: Position, sigs: HierarchicalList): Double = {
    calcU(q, sigs) +
    (for (cmbSet <- cmbSets) yield cmbSet.gMatrix.rows.toDouble * log(sigs(cmbSets.indexOf(cmbSet)))).reduceLeft(_ + _) +
    (for (cmbSet <- cmbSets) yield cmbSet.gMatrix.rows.toDouble).reduceLeft(_ + _) * log(2.0 * Pi) / 2.0
  }

  def mLogLikeGrad(q: Position, sigs: HierarchicalList): DenseVector[Double] = {
    val k = q.length + sigs.length
    val mllGrad = DenseVector.zeros[Double](k)
    mllGrad(0 to sigs.length - 1) := DenseVector.tabulate(sigs.length){i => cmbSets(i).gMatrix.rows.toDouble / sigs(i)}
    mllGrad(sigs.length to k - 1) := calcGradU(q, sigs)
    mllGrad
  }

  def gibbsUpdate(q: Position): HierarchicalList = {
    val noParameters = for (cmbSet <- cmbSets) yield cmbSet.residuals.length
    val sqMisfits = for (cmbSet <- cmbSets) yield cmbSet.sqMisfit(cmbSet.gMatrix, q, cmbSet.residuals)
    for (sqm <- sqMisfits) yield stdDevSample(sqm, noParameters(sqMisfits.indexOf(sqm)))
  }
}
