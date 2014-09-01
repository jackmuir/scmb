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

import breeze.linalg._

import scmb.shdatasets._

import scmb.hmc._


class cmbHMC(ll: Int, ee: Double, maxl: Int, dataSets: List[Map[String, String]]) extends HMC(ll, ee) {
  val cmbSets = for (dataSet <- dataSets) yield {
    match dataSet("type") {
      case "residuals" => new ResTT(maxl, dataSet)
      case "PcP-P" => new PcPmP(maxl, dataSet)
      case "P4KP-PcP" => new P4KPmPcP(maxl, dataSet)
    }
  }

  def calcU(q: DenseVector[Double], sigs: List[Double]): Double = {
    val pSums = for (cmbSet <- cmbSets) yield cmbSet.pCalcU(cmbSet.gm, q, cmbSet.residuals, sigs(cmbSets.indexOf(cmbSet)))
    sum(pSums)
    }

  def calcUgrad(q: DenseVector[Double], sigs: List[Double]): DenseVector[Double] = {
    val pGrads = for (cmbSet <- cmbSets) yield cmbSet.pCalcUgrad(cmbSet.gm, q, cmbSet.residuals, sigs(cmbSets.indexOf(cmbSet)))
    pGrads.reduceLeft((a, b) => a + b)
  }

  def gibbsUpdate(q: DenseVector[Double]): List[Double] = {
    val ns = for (cmbSet <- cmbSets) yield cmbSet.residuals.length
    val sqMisfits = for (cmbSet <- cmbSets) yield cmbSet.sqMisfit(cmbSet.gm, q, cmbSet.residuals)
    for (sqm <- sqMisfits) yield stdDevSample(sqm, ns(sqMisfits.indexOf(sqm)))
  }
}

/*
class resSHA(ll: Int, ee: Double, maxl: Int, dataSets: List[Map[String, String]]) extends HMC(ll, ee) {
  val resSets = for (dataSet <- dataSets) yield new ResTT(maxl, dataSet)

  def calcU(q: DenseVector[Double], sigs: List[Double]): Double = {
    val pSums = for (resSet <- resSets) yield resSet.pCalcU(resSet.gm, q, resSet.residuals, sigs(resSets.indexOf(resSet)))
    sum(pSums)
    }

  def calcUgrad(q: DenseVector[Double], sigs: List[Double]): DenseVector[Double] = {
    val pGrads = for (resSet <- resSets) yield resSet.pCalcUgrad(resSet.gm, q, resSet.residuals, sigs(resSets.indexOf(resSet)))
    pGrads.reduceLeft((a, b) => a + b)
  }

  def gibbsUpdate(q: DenseVector[Double]): List[Double] = {
    val ns = for (resSet <- resSets) yield resSet.residuals.length
    val sqMisfits = for (resSet <- resSets) yield resSet.sqMisfit(resSet.gm, q, resSet.residuals)
    for (sqm <- sqMisfits) yield stdDevSample(sqm, ns(sqMisfits.indexOf(sqm)))
  }
}
*/
