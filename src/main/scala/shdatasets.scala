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
package scmb.shdatasets

import breeze.linalg.{DenseVector, DenseMatrix}

import breeze.stats.distributions.Gamma

import math.{pow, abs, sqrt, Pi, cos, sin}

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

  def pCalcGradU(g: DenseMatrix[Double], q: DenseVector[Double], d: DenseVector[Double], sig: Double): DenseVector[Double] = {
    val mis = misfit(g,q,d)
    (g.t * mis) / (sig * sig)
  }

  def stringListToDouble(strList: List[String]): List[Double] = for (str <- strList if str.isEmpty == false) yield str.toDouble

  def splitAtLengthOne[T](listWithLengthOnes: List[List[T]]): (List[List[List[T]]], List[T]) = {
    def splitAcc[T](result: List[List[List[T]]], rejects: List[T], remainder: List[List[T]]): (List[List[List[T]]], List[T]) = {
      remainder.length match {
        case 0 => (result, rejects)
        case _ => {
          val (rem, resPlus) = remainder.splitAt(remainder.lastIndexWhere(_.length == 1))
          splitAcc(resPlus.tail :: result, resPlus.head.head :: rejects, rem)
        }
      }
    }
    splitAcc(Nil, Nil, listWithLengthOnes)
  }

  def fortranInConvert(filename: String): List[String] = Source.fromFile(filename).getLines().toList.drop(1)

  def fileToListListDouble(infile: List[String]): List[List[Double]] = for(line <- infile) yield stringListToDouble(line.split(" ").toList)

  def oTimesLoader(filename: String): List[Double] = {
    //otimes files have file data on the first line that we don't want; hence the drop
    val oTimesFile = fortranInConvert(filename)
    val splitOTimesFile = fileToListListDouble(oTimesFile)
    for (line <- splitOTimesFile) yield line(1) + line(2)
  }

  def pathLoader(filename: String): (List[List[List[Double]]], List[Double]) = {
    //raypath files have file data on the first line that we don't want; hence the drop
    val pathFile = fortranInConvert(filename)
    val splitPathFile = fileToListListDouble(pathFile)
    splitAtLengthOne(splitPathFile)
  }

  def combineInPairs[T](listToPair: List[T]): List[List[T]] = {
    assert(listToPair.length % 2 == 0, "List must be of even length")
    def pairUp[T](results: List[List[T]], remainder: List[T]): List[List[T]] = {
      remainder.length match {
        case 0 => results
        case _ => pairUp(remainder.take(2) :: results, remainder.drop(2))
      }
    }
    pairUp(Nil, listToPair).reverse
  }

  def getPathsAndReflections(dataM: Map[String, String]): (List[List[List[Double]]], List[List[List[List[Double]]]], List[List[Double]]) = {
    val rayPaths = pathLoader(dataM("raypaths"))._1
    val (topoRefs, rayParams) = pathLoader(dataM("topoin"))
    (rayPaths, combineInPairs(topoRefs), combineInPairs(rayParams))
  }

  def checkDataContains(dataM: Map[String,String]) {
    assert(dataM.keySet.exists(_ == "otimes"), "Data map must contain an otimes file")
    assert(dataM.keySet.exists(_ == "raypaths"), "Data map must contain a raypath file")
    assert(dataM.keySet.exists(_ == "topoin"), "Data map must contain a topoin file")
  }
}

class ResTT(maxl: Int, dataM: Map[String,String]) extends SHData(maxl) {
  assert(dataM.keySet.exists(_ == "midpoints"), "Data map must contain a midpoints file")
  val midPoints = midpointsLoader(dataM("midpoints"))
  val residuals = DenseVector(midPoints.transpose.apply(5).toArray)
  val gm = gMatrix(midPoints,hList)

  def midpointsLoader(filename: String): List[List[Double]] = {
    val midpointsFile = fortranInConvert(filename)
    fileToListListDouble(midpointsFile)
  }

  def gMatrix(mp: List[List[Double]], hl: List[List[Int]]): DenseMatrix[Double] = {
    DenseMatrix.tabulate(mp.length, hl.length){case (i, j) => llcylm(hl(j)(0), hl(j)(1), mp(i)(6), mp(i)(7))}
  }
}


class PcPmP(maxl: Int, dataM: Map[String,String]) extends SHData(maxl) {
  checkDataContains(dataM)
  val otimes = oTimesLoader(dataM("otimes"))
  val (rayPaths, topoRefs, rayParams) = getPathsAndReflections(dataM)
}

class P4KPmPcP(maxl: Int, dataM: Map[String,String]) extends SHData(maxl) {
  checkDataContains(dataM)
  val otimes = oTimesLoader(dataM("otimes"))
  val (rayPaths, topoRefs, rayParams) = getPathsAndReflections(dataM)
}
