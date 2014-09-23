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

import breeze.linalg.{DenseVector, min}

import breeze.stats.distributions.{Uniform, Gaussian, Gamma}

import math.{log, sqrt}

import scala.annotation.tailrec

/*
This algorithm, with the described alterations (perturbing epsilon and L,
performing Gibbs updates of heirarchical variables) may be found in R. Neal
`` MCMC using Hamiltonian Dynamics'' http://arxiv.org/abs/1206.1901
*/

abstract class HMC(baseSteps: Int, baseEpsilon: Double) {
  type Position = DenseVector[Double]
  type Momentum = DenseVector[Double]
  type HierarchicalList = List[Double]
  val uniformDev = new Uniform(0.0,1.0)
  val normalDev = new Gaussian(0.0,1.0)

  final def metCheck(ho: Double, hn: Double): Boolean = log(uniformDev.sample()) < ho - hn

  def calcU(q: Position, heirarchicalVar: HierarchicalList): Double

  def calcGradU(q: Position, hierarchicalVar: HierarchicalList): DenseVector[Double]

  def misfits(q: Position): List[DenseVector[Double]]

  def varianceReduction(q: Position): List[Double]

  def mLogLike(q: Position, heirarchicalVar: HierarchicalList): Double

  def mLogLikeGrad(q: Position, heirarchicalVar: HierarchicalList): DenseVector[Double]

  final def totalEnergy(q: Position, p: Momentum, hierarchicalVar: HierarchicalList) = calcU(q, hierarchicalVar) + (p dot p) / 2.0

  def gibbsUpdate(q: Position): HierarchicalList

  final def hmcUpdate(q: Position, hierarchicalVar: HierarchicalList): (Position, Int) = {
    val minHeirarchical = min(min(hierarchicalVar), 1.0) // making maximum multiplier 1 stops the integration from going wildly unstable when in low probability areas
    val epsilon = baseEpsilon * minHeirarchical * minHeirarchical * (0.9 + 0.2 * uniformDev.sample()) // perturb epsilon pm 10% to avoid cyclic orbits
    val steps = (baseSteps.toDouble * (0.9 + 0.2 * uniformDev.sample())).toInt
    val p = DenseVector(normalDev.sample(q.length).toArray)
    val ho = totalEnergy(q, p, hierarchicalVar)
    @tailrec
    def leapFrog(qo: Position, po: Momentum, stepsToGo: Int): (Position, Momentum) = stepsToGo match {
      case 1 => {
        val qn = qo + (po :* epsilon)
        val pn = po - (calcGradU(qn, hierarchicalVar) :* epsilon) / 2.0
        (qn, pn)
      }
      case _ => {
        val qn = qo + (po :* epsilon)
        val pn = po - (calcGradU(qn, hierarchicalVar) :* epsilon)
        leapFrog(qn, pn, stepsToGo - 1)
      }
    }

    val (qn, pn) = leapFrog(q, p - (calcGradU(q, hierarchicalVar) :* epsilon) / 2.0, steps)
    val hn = totalEnergy(qn, pn, hierarchicalVar)
    if (metCheck(ho, hn)) (qn, 1) else (q, 0)
  }

  final def stdDevSample(sqMis: Double, freeParameters: Int) = sqrt(
    sqMis / (2.0 * Gamma((freeParameters.toDouble + 3.0) / 2.0, 1.0).sample())
    )

  //Pass initial values to rl
  @tailrec
  final def hmcRunner(maxIt: Int, report: Int, iter: Int, acc: Int,
    rl: List[(Position, HierarchicalList, Double, List[Double])]):
    List[(Position, HierarchicalList, Double, List[Double])] = iter match {
    case `maxIt` => {
      val (q, acceptedNew) = hmcUpdate(rl(0)._1, rl(0)._2)
      val hierarchicalVar = gibbsUpdate(q)
      val mLL = mLogLike(q, hierarchicalVar)
      val vR = varianceReduction(q)
      (q, hierarchicalVar, mLL, vR) :: rl
    }
    case _ => {
      val (q, acceptedNew) = hmcUpdate(rl(0)._1, rl(0)._2)
      val hierarchicalVar = gibbsUpdate(q)
      val currU = calcU(q, hierarchicalVar)
      val mLL = mLogLike(q, hierarchicalVar)
      val vR = varianceReduction(q)
      if (iter % report == 0) {
        println(s"\nIteration $iter of $maxIt")
        println(s"U: $currU")
        println(s"Acceptance: ${100.0 * (acc + acceptedNew).toDouble / iter.toDouble}")
      }
      hmcRunner(maxIt, report, iter + 1, acc + acceptedNew, (q, hierarchicalVar, mLL, vR) :: rl)
    }
  }
}
