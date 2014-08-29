import breeze.linalg._

import breeze.stats.distributions.Uniform

import breeze.stats.distributions.Gaussian

class HMC {

  val uDev = new Uniform(0.0,1.0)

  val nDev = new Gaussian(0.0,1.0)

  def metCheck(ho: Double, hn: Double): Boolean = math.log(uDev.sample()) < ho - hn

  // For now calcU and calcUgrad return 0's, not sure what a good default is
  def calcU(q: DenseVector[Double]): Double = 0

  def calcUgrad(q: DenseVector[Double]): DenseVector[Double] = DenseVector.zeros[Double](q.length)

  def hmcUpdate(qo: DenseVector[Double], l: Int, ep: Double, acc: Int): (DenseVector[Double], Int) = {
    val uep = 0.9 * ep + 0.2 * uDev.sample() // perturb epsilon pm 10% to avoid cyclic orbits
    val po = DenseVector(nDev.sample(qo.length).toArray)
    val ho = calcU(qo) + (po dot po)

    def leapFrog(qo: DenseVector[Double], po: DenseVector[Double], l: Int, ep: Double): (DenseVector[Double], DenseVector[Double]) = {
      if (l == 1) {
        val qn = qo + (po :* ep)
        val pn = qo - (calcUgrad(qo) :* ep) / 2.0
        (qn, pn)
      }
      else {
        val qn = qo + (po :* ep)
        val pn = qo - (calcUgrad(qo) :* ep)
        leapFrog(qn, pn, l - 1, ep)
      }
    }

    val (qn, pn) = leapFrog(qo, po - (calcUgrad(qo) :* uep) / 2.0, l, uep)
    val hn = calcU(qn) + (pn dot pn)
    if (metCheck(ho,hn)) (qn, acc + 1) else (qo, acc)
  }
}
