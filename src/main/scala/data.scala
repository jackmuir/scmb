import breeze.linalg._

import breeze.stats.distributions.Gamma

import breeze.constants.Pi

import math._

class SHDataObject(maxl: Int) {

  val ml = maxl

  val hlist = harmlister(ml, Nil)

  def ylm(l: Int, m: Int, th: Double, ph: Double): Double = {
    val pi = 3.1415926536
    //compute (m)!!/sqrt((m)!) (supply 2* m for xfact fortran code)
    def xfact(m: Int): Double = {
      if (m <= 1) 1.0
      else {
        if (m % 2 == 1) m.toDouble / sqrt(m.toDouble) * xfact(m - 1)
        else 1.0 / sqrt(m.toDouble) * xfact(m - 1)
      }
    }
    //this is a very un-scala function.... + no exception handling
    def modpLgndr(l: Int, m: Int, x: Double): Double = {
      val dl = l.toDouble
      val dm = m.toDouble
      val norm = sqrt(2.0 * dl + 1.0) / sqrt(4.0 * pi)
      var pmm = norm
      if (m != 0) pmm = (pow(-1, m)).toDouble * pmm * xfact(2 * m) * pow((1.0-x * x), (dm / 2.0))
      if (l == m) pmm
      else {
        var pmmp1 = x * pmm * sqrt(2.0 * m + 1.0)
        if (l == m + 1) pmmp1
        else {
          var pll = 0.0
          var dll = 0.0
          for (ll <- m + 2 to l) {
            dll = ll.toDouble
            pll = (x * (2.0 * dll - 1.0) * pmmp1 - sqrt(pow((dll - 1.0), 2.0) - dm * dm) * pmm) / sqrt(pow(dll, 2.0) - pow(dm, 2.0))
            pmm = pmmp1
            pmmp1 = pll
          }
          pll
        }
      }
    }

    if (m > 0) modpLgndr(l, m, cos(th)) * cos(m * ph) * sqrt(2.0)
    else if (m < 0 ) modpLgndr(l, -m, cos(th)) * sin(m * ph) * sqrt(2.0) else modpLgndr(l, m, cos(th))

  }

  def harmLister(l: Int, hl: List[List[Int]]): List[List[Int]] = {
    if (l == 0) List(0, 0) :: hl
    else {
      val sl = (for(m <- -l to l) yield List(l,m)).toList
      harmlister(l - 1, sl ::: hl)
    }
  }

  def misfit(g: DenseMatrix[Double], q: DenseVector[Double], m: DenseVector[Double]): DenseVector[Double] = g * q - m

  def sqMisfit(g: DenseMatrix[Double], q: DenseVector[Double], m: DenseVector[Double]): Double = {
    val mis = misfit(g,q,m)
    mis dot mis
  }

  def pcalcU(g: DenseMatrix[Double], q: DenseVector[Double], m: DenseVector[Double], sig: Double): Double = {
    val phi = sqMisfit(g,q,m)
    phi / (2.0 * sig)
  }

  def pcalcUgrad(g: DenseMatrix[Double], q: DenseVector[Double], m: DenseVector[Double], sig: Double): DenseVector[Double] = {
    val mis = misfit(g,q,m)
    g.t * mis / sig
  }

  def stdDevSample(sqmis: Double, n: Int) = sqrt(
    sqmis / (2.0 * Gamma((n.toDouble + 3) / 2.0, 1.0).sample())
    )
}
