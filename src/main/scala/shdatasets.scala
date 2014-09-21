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
  val gMatrix: DenseMatrix[Double]
  val residuals: DenseVector[Double]
  val rCore = 3479.5
  val vPlus = 13.6601
  val vMinus = 8.0

  def topLeg(rayParam: Double): Double = - sqrt(rCore * rCore / (vPlus * vPlus) - rayParam * rayParam) / rCore

  def bottomLeg(rayParam: Double): Double = sqrt(rCore * rCore / (vMinus * vMinus) - rayParam * rayParam) / rCore

  def reflectTop(rayParam: Double): Double = 2.0 * topLeg(rayParam)

  def reflectBottom(rayParam: Double): Double = 2.0 * bottomLeg(rayParam)

  def transmitThrough(rayParam: Double): Double = topLeg(rayParam) + bottomLeg(rayParam)

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
    else if (m < 0 ) modpLgndr(l, -m, cos(th)) * sin(-m * ph) * sqrt(2.0) else modpLgndr(l, m, cos(th))
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

  def oTimesLoader(filename: String): Array[Double] = {
    //otimes files have file data on the first line that we don't want; hence the drop
    val oTimesFile = fortranInConvert(filename)
    val splitOTimesFile = fileToListListDouble(oTimesFile)
    //(for (line <- splitOTimesFile) yield line(1) + line(2)).toArray
    // lets try just the deviation...
    (for (line <- splitOTimesFile) yield line(1)).toArray
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

  def getPathsAndReflections(dataM: Map[String, String]): (List[List[List[Double]]], List[List[List[Double]]], List[Double]) = {
    val rayPaths = pathLoader(dataM("raypaths"))._1
    val (topoRefs, rayParams) = pathLoader(dataM("topoin"))
    (rayPaths, topoRefs, rayParams)
  }

  def checkDataContains(dataM: Map[String,String]) {
    assert(dataM.keySet.exists(_ == "otimes"), "Data map must contain an otimes file")
    assert(dataM.keySet.exists(_ == "raypaths"), "Data map must contain a raypath file")
    assert(dataM.keySet.exists(_ == "topoin"), "Data map must contain a topoin file")
  }

  def gTomoMatrix(rayPaths: List[List[List[Double]]]): DenseMatrix[Double] = {
    def gTomoMatElement(selectedRayPath: Int, harmonic: Int): Double = {
      val List(l, m) = hList(harmonic)
      val rayPath = rayPaths(selectedRayPath)
      val pathSubSums = for (pathElement <- rayPath) yield pathElement(3) * llcylm(l, m, pathElement(1), pathElement(2))
      (0.0 /: pathSubSums) (_ + _)
    }
    DenseMatrix.tabulate(rayPaths.length, hList.length){case (i, j) => gTomoMatElement(i, j)}
  }

}

class ResTT(maxl: Int, dataM: Map[String,String]) extends SHData(maxl) {
  assert(dataM.keySet.exists(_ == "midpoints"), "Data map must contain a midpoints file")
  val midPoints = midpointsLoader(dataM("midpoints"))
  val residuals = DenseVector(midPoints.transpose.apply(5).toArray)
  val gMatrix = DenseMatrix.tabulate(midPoints.length, hList.length){
    case (i, j) => llcylm(hList(j)(0), hList(j)(1), midPoints(i)(6), midPoints(i)(7))
    }

  def midpointsLoader(filename: String): List[List[Double]] = {
    val midpointsFile = fortranInConvert(filename)
    fileToListListDouble(midpointsFile)
  }
}


class PcPmP(maxl: Int, dataM: Map[String,String]) extends SHData(maxl) {
  checkDataContains(dataM)
  val residuals = DenseVector(oTimesLoader(dataM("otimes")))
  val (rayPaths, topoRefs, rayParams) = getPathsAndReflections(dataM)
  val gTopoMatrix = {
    def gTopoMatElement(selectedTopo: Int, harmonic: Int): Double = {
      val List(l, m) = hList(harmonic)
      val topoLat = topoRefs(selectedTopo)(0)(3)
      val topoLon = topoRefs(selectedTopo)(0)(4)
      val rayParam = rayParams(selectedTopo)
      llcylm(l, m, topoLat, topoLon) * reflectTop(rayParam)
    }
    DenseMatrix.tabulate(topoRefs.length, hList.length){case (i, j) => gTopoMatElement(i, j)}
  }
////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!
  val gMatrix = DenseMatrix.horzcat(gTomoMatrix(rayPaths), gTopoMatrix :* 1000.0)
// times 1000 to approximately equalize elements; have to multiply resulting parameters by 1000 to get true
// radial values
}

class P4KPmPcP(maxl: Int, dataM: Map[String,String]) extends SHData(maxl) {
  checkDataContains(dataM)
  val residuals = DenseVector(oTimesLoader(dataM("otimes")))
  val (rayPaths, topoRefs, rayParams) = getPathsAndReflections(dataM)
  val pairedTopoRefs = combineInPairs(topoRefs)
  val pairedRayParams = combineInPairs(rayParams)
  val gTopoMatrix = {
    def gTopoMatElement(selectedTopo: Int, harmonic: Int): Double = {
      val List(l, m) = hList(harmonic)
      val p4kptopoLat = for (topo <- pairedTopoRefs(selectedTopo)(0)) yield topo(3)
      val p4kptopoLon = for (topo <- pairedTopoRefs(selectedTopo)(0)) yield topo(4)
      val p4kprayParam = pairedRayParams(selectedTopo)(0)
      val pcptopoLat = pairedTopoRefs(selectedTopo)(1)(0)(3)
      val pcptopoLon = pairedTopoRefs(selectedTopo)(1)(0)(4)
      val pcprayParam = pairedRayParams(selectedTopo)(1)
      // Ugly, but is there really a better way to do it?
      val len = p4kptopoLat.length
      val n = (len - 3) / 2
      // This is the diffracted case; the transmission coefficients are at 2 different locations
      // and so are split into top and bottom legs
      //ACTUALLY THIS CAPTURES BOTH DIFFRACTED AND NON DIFFRACTED CASES
      llcylm(l, m, p4kptopoLat(0), p4kptopoLon(0)) * topLeg(p4kprayParam) +
      llcylm(l, m, p4kptopoLat(n - 1), p4kptopoLon(n - 1)) * bottomLeg(p4kprayParam) +
      llcylm(l, m, p4kptopoLat(n), p4kptopoLon(n)) * reflectBottom(p4kprayParam) +
      llcylm(l, m, p4kptopoLat(n + 1), p4kptopoLon(n + 1)) * reflectBottom(p4kprayParam) +
      llcylm(l, m, p4kptopoLat(n + 2), p4kptopoLon(n + 2)) * reflectBottom(p4kprayParam) +
      llcylm(l, m, p4kptopoLat(n + 3), p4kptopoLon(n + 3)) * bottomLeg(p4kprayParam) +
      llcylm(l, m, p4kptopoLat(len - 1), p4kptopoLon(len - 1)) * topLeg(p4kprayParam) -
      llcylm(l, m, pcptopoLat, pcptopoLon) * reflectTop(pcprayParam)
      }
    }
    DenseMatrix.tabulate(pairedTopoRefs.length, hList.length){case (i, j) => gTopoMatElement(i, j)}
  }
  ////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!
  val gMatrix = DenseMatrix.horzcat(gTomoMatrix(rayPaths), gTopoMatrix :* 1000.0)
  // times 1000 to approximately equalize elements; have to multiply resulting parameters by 1000 to get true
  // radial values
}


class PKPabmPKPbc(maxl: Int, dataM: Map[String,String]) extends SHData(maxl) {
  checkDataContains(dataM)
  val residuals = DenseVector(oTimesLoader(dataM("otimes")))
  val (rayPaths, topoRefs, rayParams) = getPathsAndReflections(dataM)
  val pairedTopoRefs = combineInPairs(topoRefs)
  val pairedRayParams = combineInPairs(rayParams)
  val gTopoMatrix = {
    def gTopoMatElement(selectedTopo: Int, harmonic: Int): Double = {
      val List(l, m) = hList(harmonic)
      val pkpabtopoLat = for (topo <- pairedTopoRefs(selectedTopo)(0)) yield topo(3)
      val pkpabtopoLon = for (topo <- pairedTopoRefs(selectedTopo)(0)) yield topo(4)
      val pkpabrayParam = pairedRayParams(selectedTopo)(0)
      val pkpbctopoLat = for (topo <- pairedTopoRefs(selectedTopo)(1)) yield topo(3)
      val pkpbctopoLon = for (topo <- pairedTopoRefs(selectedTopo)(1)) yield topo(4)
      val pcprayParam = pairedRayParams(selectedTopo)(1)
      llcylm(l, m, pkpabtopoLat(0), pkpabtopoLon(0)) * transmitThrough(pkpabprayParam) +
      llcylm(l, m, pkpabtopoLat(1), pkpabtopoLon(1)) * transmitThrough(pkpabprayParam) -
      llcylm(l, m, pkpbctopoLat(0), pkpbctopoLon(0)) * transmitThrough(pkpabprayParam) -
      llcylm(l, m, pkpbctopoLat(1), pkpbctopoLon(1)) * transmitThrough(pkpabprayParam)
      }
    }
    DenseMatrix.tabulate(pairedTopoRefs.length, hList.length){case (i, j) => gTopoMatElement(i, j)}
  }
  ////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!
  val gMatrix = DenseMatrix.horzcat(gTomoMatrix(rayPaths), gTopoMatrix :* 1000.0)
  // times 1000 to approximately equalize elements; have to multiply resulting parameters by 1000 to get true
  // radial values
}
