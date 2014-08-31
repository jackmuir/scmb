package scmb.shdatasets

import breeze.linalg._

import breeze.stats.distributions.Gamma

import math.pow, math.abs, math.sqrt, math.Pi, math.cos, math.sin

import scala.io.Source

abstract class SHData(maxl: Int) {
  val ml = maxl
  val hList = harmLister(ml, Nil)

  def ylm(l: Int, m: Int, ph: Double, th: Double): Double = {
    //compute (m)!!/sqrt((m)!) (supply 2* m for xfact fortran code)
    def xfact(m: Int): Double = {
      if (m <= 1) 1.0
      else {
        if (m % 2 == 1) m.toDouble / sqrt(m.toDouble) * xfact(m - 1)
        else 1.0 / sqrt(m.toDouble) * xfact(m - 1)
      }
    }
    // notice lat lon swapped relative to ph th to conform with geographical coordinate standards

    def modpLgndr(l: Int, m: Int, x: Double): Double = {
      assert(0 <= m && m <= l && abs(x) <= 1.0)
      val dl = l.toDouble
      val dm = m.toDouble
      val norm = sqrt(2.0 * dl + 1.0) / sqrt(4.0 * Pi)
      val pmm = if (m == 0) norm else pow(-1, m).toDouble * norm * xfact(2 * m) * pow((1.0-x * x), (dm / 2.0))
      if (l == m) pmm
      else {
        val pmmp1 = x * pmm * sqrt(2.0 * m + 1.0)
        if (l == m + 1) pmmp1
        else {
          def mplacc(ll: Int, acc1: Double, acc2: Double): Double = {
            val dll = ll.toDouble
            val pll = (x * (2.0 * dll - 1.0) * acc2 - sqrt(pow((dll - 1.0), 2.0) - dm * dm) * acc1) / sqrt(pow(dll, 2.0) - pow(dm, 2.0))
            if (ll == l) pll
            else mplacc(ll + 1, acc2, pll)
          }
          mplacc(m + 2, pmm, pmmp1)
        }
      }
    }

    if (m > 0) modpLgndr(l, m, cos(th)) * cos(m * ph) * sqrt(2.0)
    else if (m < 0 ) pow(-1,m).toDouble * modpLgndr(l, -m, cos(th)) * sin(-m * ph) * sqrt(2.0) else modpLgndr(l, m, cos(th))
  }

  def llcylm(l: Int, m: Int, lat: Double, lon: Double): Double = ylm(l, m, Pi * lon / 180.0, -Pi * (lat - 90) / 180.0)

  def harmLister(l: Int, hl: List[List[Int]]): List[List[Int]] = {
    if (l == 0) List(0, 0) :: hl
    else {
      val sl = (for(m <- -l to l) yield List(l,m)).toList
      harmLister(l - 1, sl ::: hl)
    }
  }

  def misfit(g: DenseMatrix[Double], q: DenseVector[Double], d: DenseVector[Double]): DenseVector[Double] = g * q - d

  def sqMisfit(g: DenseMatrix[Double], q: DenseVector[Double], d: DenseVector[Double]): Double = {
    val mis = misfit(g,q,d)
    mis dot mis
  }

  def pCalcU(g: DenseMatrix[Double], q: DenseVector[Double], d: DenseVector[Double], sig: Double): Double = {
    val phi = sqMisfit(g,q,d)
    phi / (2.0 * sig * sig)
  }

  def pCalcUgrad(g: DenseMatrix[Double], q: DenseVector[Double], d: DenseVector[Double], sig: Double): DenseVector[Double] = {
    val mis = misfit(g,q,d)
    (g.t * mis) / (sig * sig)
  }

  def stringListToDouble(strList: List[String]): List[Double] = for (str <- strList) yield str.toDouble

  def oTimesLoader(filename: String): List[List[Double]] = {
    //otimes files have file data on the first line that we don't want; hence the drop
    val oTimesFile = Source.fromFile(filename).getLines().toList.drop(1)
    for(line <- oTimesFile) yield stringListToDouble(line.split(" ").toList)
  }

}

class ResTT(maxl: Int, dataM: Map[String,String]) extends SHData(maxl) {
  assert(dataM.keySet.exists(_ == "midpoints"), "Data map must contain a midpoints file")
  val midPoints = midpointsLoader(dataM("midpoints"))
  val residuals = DenseVector(midPoints.transpose.apply(5).toArray)
  val gm = gMatrix(midPoints,hList)

  def midpointsLoader(filename: String): List[List[Double]] = {
    val midpointsFile = Source.fromFile(filename).getLines().toList
    for(line <- midpointsFile) yield stringListToDouble(line.split(" ").toList)
  }

  def gMatrix(mp: List[List[Double]], hl: List[List[Int]]): DenseMatrix[Double] = {
    DenseMatrix.tabulate(mp.length, hl.length){case (i, j) => llcylm(hl(j)(0), hl(j)(1), mp(i)(6), mp(i)(7))}
  }
}
