import breeze.linalg._

import breeze.stats.distributions.Gamma

class DataObject(maxl: Int) {

  val ml = maxl

  val hlist = harmlister(ml, Nil)

  def harmlister(l: Int, hl: List[List[Int]]): List[List[Int]] = {
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

  def stdDevSample(sqmis: Double, n: Int) = math.sqrt(
    sqmis / (2.0 * Gamma((n.toDouble + 3) / 2.0, 1.0).sample())
    )
}
