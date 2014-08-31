package scmb.hmc

import breeze.linalg._

import breeze.stats.distributions.Uniform

import breeze.stats.distributions.Gaussian

import breeze.stats.distributions.Gamma

import math.log, math.sqrt

import scala.annotation.tailrec

import scmb.shdatasets._

abstract class HMC(ul: Int, ue: Double) {
  val uDev = new Uniform(0.0,1.0)
  val nDev = new Gaussian(0.0,1.0)
  val l = ul
  val ep = ue

  private def metCheck(ho: Double, hn: Double): Boolean = log(uDev.sample()) < ho - hn

  def calcU(q: DenseVector[Double], sigs: List[Double]): Double

  def calcUgrad(q: DenseVector[Double], sigs: List[Double]): DenseVector[Double]

  def gibbsUpdate(q: DenseVector[Double]): List[Double]

  final def hmcUpdate(qo: DenseVector[Double], sigs: List[Double], acc: Int): (DenseVector[Double], Int) = {
    val uep = ep * min(sigs) * min(sigs) * (0.9 + 0.2 * uDev.sample()) // perturb epsilon pm 10% to avoid cyclic orbits
    val ul = (l.toDouble * (0.9 + 0.2 * uDev.sample())).toInt
    val po = DenseVector(nDev.sample(qo.length).toArray)
    val ho = calcU(qo, sigs) + (po dot po)
    @tailrec
    def leapFrog(qo: DenseVector[Double], po: DenseVector[Double], ll: Int): (DenseVector[Double], DenseVector[Double]) = {
      if (ll == 1) {
        val qn = qo + (po :* uep)
        val pn = po - (calcUgrad(qn, sigs) :* uep) / 2.0
        (qn, pn)
      }
      else {
        val qn = qo + (po :* uep)
        val pn = po - (calcUgrad(qn, sigs) :* uep)
        leapFrog(qn, pn, ll - 1)
      }
    }

    val (qn, pn) = leapFrog(qo, po - (calcUgrad(qo, sigs) :* uep) / 2.0, ul)
    val hn = calcU(qn, sigs) + (pn dot pn)
    if (metCheck(ho,hn)) (qn, acc + 1) else (qo, acc)
  }

  def stdDevSample(sqmis: Double, n: Int) = sqrt(
    sqmis / (2.0 * Gamma((n.toDouble + 3) / 2.0, 1.0).sample())
    )

  //Pass initial values to rl
  @tailrec
  final def hmcRunner(qlen: Int, maxIt: Int, report: Int, iter: Int, acc: Int,
    rl: List[(DenseVector[Double], List[Double])]):
    List[(DenseVector[Double], List[Double])] = {
    if (iter == maxIt) {
      val (qn, accn) = hmcUpdate(rl(0)._1, rl(0)._2, acc)
      val sigsn = gibbsUpdate(qn)
      (qn, sigsn) :: rl
    }
    else {
      val (qn, accn) = hmcUpdate(rl(0)._1, rl(0)._2, acc)
      val sigsn = gibbsUpdate(qn)
      val currU = calcU(qn, sigsn)
      if (iter % report == 0) {
        print("\nIteration "); print(iter); print(" of "); println(maxIt)
        print("U: "); println(currU)
        print("Acceptance: "); println(acc.toDouble / iter.toDouble)
        print("Sigmas: "); println(min(sigsn))
      }
      hmcRunner(qlen, maxIt, report, iter + 1, accn, (qn, sigsn) :: rl)
    }
  }
}

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


//class sCMB(datasets: List[Map[String, String]]) extends HMC {}
