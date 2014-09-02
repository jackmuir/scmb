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

import breeze.linalg.DenseVector

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

  def calcU(q: Position, heirarchicalVar: List[Double]): Double

  def calcGradU(q: Position, hierarchicalVar: HierarchicalList): DenseVector[Double]

  final def totalEnergy(q: Position, p: Momentum, hierarchicalVar: HierarchicalList) = calcU(q, hierarchicalVar) + (p dot p) / 2.0

  def gibbsUpdate(q: Position): HeirarchicalList

  final def hmcUpdate(q: Position, hierarchicalVar: HierarchicalList): (Position, Int) = {
    val minHeirarchical = min(hierarchicalVar)
    val epsilon = baseEpsilon * minHeirarchical * minHeirarchical * (0.9 + 0.2 * uniformDev.sample()) // perturb epsilon pm 10% to avoid cyclic orbits
    val steps = (baseSteps.toDouble * (0.9 + 0.2 * uDev.sample())).toInt
    val p = DenseVector(normalDev.sample(q.length).toArray)
    val ho = totalEnergy(q, hierarchicalVar, p)
    @tailrec
    def leapFrog(qo: Position, po: Momentum, stepsToGo: Int): (Position, Momentum) = match stepsToGo {
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

    val (qn, pn) = leapFrog(qo, po - (calcGradU(qo, hierarchicalVar) :* epsilon) / 2.0, steps)
    val hn = totalEnergy(qn, hierarchicalVar, pn)
    if (metCheck(ho, hn)) (qn, 1) else (qo, 0)
  }

  final def stdDevSample(sqMis: Double, freeParameters: Int) = sqrt(
    sqMis / (2.0 * Gamma((freeParameters.toDouble + 3.0) / 2.0, 1.0).sample())
    )

  //Pass initial values to rl
  @tailrec
  final def hmcRunner(maxIt: Int, report: Int, iter: Int, acc: Int,
    rl: List[(Position, HierarchicalList)]):
    List[(Position, HierarchicalList)] = iter match {
    case maxit => {
      val (qn, acceptedNew) = hmcUpdate(rl(0)._1, rl(0)._2)
      val hierarchicalVar = gibbsUpdate(qn)
      (qn, hierarchicalVar) :: rl
    }
    case _ => {
      val (qn, acceptedNew) = hmcUpdate(rl(0)._1, rl(0)._2)
      val hierarchicalVar = gibbsUpdate(qn)
      val currU = calcU(qn, hierarchicalVar)
      if (iter % report == 0) {
        println("\nIteration $iter of $maxIt")
        println("U: $currU")
        println("Acceptance: ${(acc + acceptedNew).toDouble / iter.toDouble}")
      }
      hmcRunner(maxIt, report, iter + 1, acc + acceptedNew, (qn, hierarchicalVar) :: rl)
    }
  }
}
