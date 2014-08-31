package scmb.hmc

import breeze.linalg._

import breeze.stats.distributions.Uniform

import breeze.stats.distributions.Gaussian

import breeze.stats.distributions.Gamma

import math.log, math.sqrt

import scmb.shdatasets._

abstract class HMC(ll: Int, ee: Double) {
  val uDev = new Uniform(0.0,1.0)
  val nDev = new Gaussian(0.0,1.0)
  val l = ll
  val ep = ee

  def metCheck(ho: Double, hn: Double): Boolean = log(uDev.sample()) < ho - hn

  def calcU(q: DenseVector[Double]): Double

  def calcUgrad(q: DenseVector[Double]): DenseVector[Double]

  def gibbsUpdate(q: DenseVector[Double]): List[Double]

  def hmcUpdate(qo: DenseVector[Double], sigs: List[Double], l: Int, ep: Double, acc: Int): (DenseVector[Double], Int) = {
    val uep = 0.9 * ep + 0.2 * uDev.sample() // perturb epsilon pm 10% to avoid cyclic orbits
    val po = DenseVector(nDev.sample(qo.length).toArray)
    val ho = calcU(qo, sigs) + (po dot po)

    def leapFrog(qo: DenseVector[Double], po: DenseVector[Double], sigs: List[Double], l: Int, ep: Double): (DenseVector[Double], DenseVector[Double]) = {
      if (l == 1) {
        val qn = qo + (po :* ep)
        val pn = qo - (calcUgrad(qo, sigs) :* ep) / 2.0
        (qn, pn)
      }
      else {
        val qn = qo + (po :* ep)
        val pn = qo - (calcUgrad(qo, sigs) :* ep)
        leapFrog(qn, pn, sigs, l - 1, ep)
      }
    }

    val (qn, pn) = leapFrog(qo, po - (calcUgrad(qo, sigs) :* uep) / 2.0, l, uep)
    val hn = calcU(qn, sigs) + (pn dot pn)
    if (metCheck(ho,hn)) (qn, acc + 1) else (qo, acc)
  }

  def stdDevSample(sqmis: Double, n: Int) = sqrt(
    sqmis / (2.0 * Gamma((n.toDouble + 3) / 2.0, 1.0).sample())
    )

  //Pass initial values to rl
  def hmcRunner(l: Int, ep: Double, qlen: Int, report: Int, iter: Int, acc: Int
    rl: List[(DenseVector[Double], List[Double], Double)]):
    List[(DenseVector[Double], List[Double], Double)] = {
    if (iter == 1) {
      val (qn, accn) = hmcUpdate(rl(0)(0), rl(0)(1), l, ep, acc)
      val sigsn = gibbsUpdate(qn)
      (qn, sign, accn) :: r
    }
    else {
      val (qn, accn) = hmcUpdate(rl(0)(0), rl(0)(1), l, ep, acc)
      val sigsn = gibbsUpdate(qn)
      hmcRunner(l, ep * min(sigsn) * min(sigsn), qlen, report, iter - 1, accn, (qn, sign, accn) :: rl)
    }
  }
}

class resSHA(maxl: Int, dataSets: List[Map[String, String]], ll: Int, ee: Double) extends HMC(ll, ee) {
  val resSets = for (dataSet <- dataSets) yield new ResTT(maxl, dataSet)

  def calcU(q: DenseVector[Double], sigs: List[Double]): Double = {
    val pSums = for (resSet <- resSets) yield resSet.pCalcU(resSet.gm, q, resSet.residuals, sigs(resSets.indexOf(resSet)))
    sum(pSums)
    }

  def calcUgrad(q: DenseVector[Double], sigs: List[Double]): DenseVector[Double] = {
    val pGrads = for (resSet <- resSets) yield resSet.pCalcUgrad(resSet.gm, q, resSet.residuals, sigs(resSets.indexOf(resSet)))
    pGrads.reduceLeft((a, b) => a + b)
  }
}


//class sCMB(datasets: List[Map[String, String]]) extends HMC {}
