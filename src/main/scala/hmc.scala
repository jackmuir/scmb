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
package scmb.hmc

import breeze.linalg._

import breeze.stats.distributions.Uniform

import breeze.stats.distributions.Gaussian

import breeze.stats.distributions.Gamma

import math.log, math.sqrt

import scala.annotation.tailrec

abstract class HMC(ul: Int, ue: Double) {
  val uDev = new Uniform(0.0,1.0)
  val nDev = new Gaussian(0.0,1.0)
  val l = ul
  val ep = ue

  final def metCheck(ho: Double, hn: Double): Boolean = log(uDev.sample()) < ho - hn

  def calcU(q: DenseVector[Double], sigs: List[Double]): Double

  def calcUgrad(q: DenseVector[Double], sigs: List[Double]): DenseVector[Double]

  def gibbsUpdate(q: DenseVector[Double]): List[Double]

  final def hmcUpdate(qo: DenseVector[Double], sigs: List[Double], acc: Int): (DenseVector[Double], Int) = {
    val uep = ep * min(sigs) * min(sigs) * (0.9 + 0.2 * uDev.sample()) // perturb epsilon pm 10% to avoid cyclic orbits
    val ul = (l.toDouble * (0.9 + 0.2 * uDev.sample())).toInt
    val po = DenseVector(nDev.sample(qo.length).toArray)
    val ho = calcU(qo, sigs) + (po dot po) / 2.0
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
    val hn = calcU(qn, sigs) + (pn dot pn) / 2.0
    if (metCheck(ho,hn)) (qn, acc + 1) else (qo, acc)
  }

  final def stdDevSample(sqmis: Double, n: Int) = sqrt(
    sqmis / (2.0 * Gamma((n.toDouble + 3) / 2.0, 1.0).sample())
    )

  //Pass initial values to rl
  @tailrec
  final def hmcRunner(maxIt: Int, report: Int, iter: Int, acc: Int,
    rl: List[(DenseVector[Double], List[Double])]):
    List[(DenseVector[Double], List[Double])] = {
    if (iter == maxIt) {
      val (qn, accn) = hmcUpdate(rl(0)._1, rl(0)._2, acc)
      val gibbsVar = gibbsUpdate(qn)
      (qn, gibbsVar) :: rl
    }
    else {
      val (qn, accn) = hmcUpdate(rl(0)._1, rl(0)._2, acc)
      val gibbsVar = gibbsUpdate(qn)
      val currU = calcU(qn, gibbsVar)
      if (iter % report == 0) {
        print("\nIteration "); print(iter); print(" of "); println(maxIt)
        print("U: "); println(currU)
        print("Acceptance: "); println(acc.toDouble / iter.toDouble)
      }
      hmcRunner(maxIt, report, iter + 1, accn, (qn, gibbsVar) :: rl)
    }
  }
}
