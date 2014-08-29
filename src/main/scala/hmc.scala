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
    val po = nDev.sample(qo.Length)
    val ho = calcU(q0) + po dot po
    val (qn, pn) = leapFrog(qo, po - uep * calcUgrad(qo) / 2.0, l, uep)
    val hn = calcU(qn) + pn dot pn

    if (metCheck(ho,hn)) (qn, acc + 1) else (qo, acc)

    def leapFrog(qo: DenseVector[Double], po: DenseVector[Double], l: Int, ep: Dobule) {
      if (l == 1) {
        val qn = qo + ep * po
        val pn = qo - ep * calcUgrad(qo) / 2.0
        (qn, pn)
      }
      else {
        val qn = qo + ep * po
        val pn = qo - ep * calcUgrad(qo)
        leapFrog(qn, pn, l - 1, ep)
      }
    }
}

  val hmcu = (hmcupdate _).tupled
}
