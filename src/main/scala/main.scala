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

/*
HMC / Gibbs MCMC based hiearchical Bayesian inversion for the lowermost mantle.

CALL SIGNATURE: (after loading sbt in the scmb directory)
run maxdeg L epsilon maxIterations burnIterations reportIterations dataFiles

INPUTS:

maxdeg[Integer] -  maximum degree of spherical harmonic expansion

L[Integer] - base number of Leapfrog steps
          (tuning parameter; if autocorrelation of output is high, increase. High L makes the program slower)

epsilon[Double] - base length of Leapfrog steps;
          (tuning parameter; if the acceptance rate is very low, decrease.
           Optimal acceptance is 65%; 20% - 90% acceptable. Lower epsilon makes the output lower quality)

maxIterations[Integer] - total number of HMC/Gibbs updates to perform. Increase until you get good statistics

burnIterations[Integer] - the code removes the first burnIteration iterations automatically to get rid of chain burn in

reportIterations[Integer] - print some statistics as to how the program is going every reportIteration iterations

dataFiles[.xml] - xml file with the paths to the data files. An example is the file datasets.xml

OUTPUTS

tomo_parameters[.dat] - file containing spherical harmonic coefficients of (slowness/velocity?) perturbations;
                        one set of coefficients per line
topo_parameters[.dat] - file containing spherical harmonic coefficients of CMB radius perturbations;
                        one set of coefficients per line
noise_parameters[.dat] - file containing the ``sigma'' noise parameters for each dataset. One set of sigmas per line

All of these files should be the same length (ie. length corresponds to the number of iterations)

GOTCHAS

The code performs better when the topography and tomography coefficients end up with similar magnitudes.
Therefore, it is necessary to multiply the topography matrix with a large constant factor
This means that both the topography and tomography coefficients will be small and similar in size.
Consequently, when printing the topography results, we need to multiply the results by this factor to get
the 'true' coefficients. Every time this happens, it is marked by a line
////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!
preceeding the multiplication. You will need to change 1 line in main.scala, and 1-3 lines in shdatasets.scala

*/

import scmb.cmbhmc.cmbHMC

import breeze.linalg.{DenseVector, min}

import java.io.{File, PrintWriter}

import xml.{XML, Node, Elem}

import breeze.optimize.{DiffFunction, LBFGS}


object sCMB {
  def main(args: Array[String]) {
    val maxdeg = args(0).toInt
    val qlen = (maxdeg + 1) * (maxdeg + 1)
    val l  = args(1).toInt
    val ep = args(2).toDouble
    val maxIt = args(3).toInt
    val burn = args(4).toInt
    val report = args(5).toInt
    def xmlToListMaps(xmldata: Elem): List[Map[String, String]] = {
      def xmlNodeToMap(node: Node): Map[String, String] = {
        Map("type" -> (node \ "type").text, "otimes" -> (node \ "otimes").text,
            "raypaths" -> (node \ "raypaths").text, "topoin" -> (node \ "topoin").text)
      }
      (for (set <- xmldata \ "set") yield xmlNodeToMap(set)).toList
    }
    val inxml = XML.loadFile(args(6))
    val datamaps = xmlToListMaps(inxml)
    val initialGuess =  List((DenseVector.zeros[Double](qlen * 2), (for (set <- datamaps) yield 1.0).toList, Double.PositiveInfinity, (for (set <- datamaps) yield 0.0).toList))
    val runner = new cmbHMC(l, ep, maxdeg, datamaps)
    val results = runner.hmcRunner(maxIt, report, 1, 0, initialGuess)

    def printParametersAndNoise(f1: File, f2: File, f3: File, f4: File) {
      val p1 = new PrintWriter(f1)
      val p2 = new PrintWriter(f2)
      val p3 = new PrintWriter(f3)
      val p4 = new PrintWriter(f4)
      for (res <- results.dropRight(burn)) {
        for (tomo <- res._1.apply(0 to qlen - 1)) p1.print(s"$tomo ")
        ////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!
        for (topo <- res._1.apply(qlen to 2 * qlen - 1)) p2.print(s"${topo * 1000.0} ") // see shdatasets for 1000
        for (noise <- res._2) p3.print(s"$noise ")
        for (vr <- res._4) p4.print(s"$vr ")
        p1.print("\n"); p2.print("\n"); p3.print("\n"); p4.print("\n")
      }
      p1.close(); p2.close(); p3.close(); p4.close()
    }

    def printAvMisfits() {
      val p = for (cmbSet <- runner.cmbSets) yield new PrintWriter(new File(s"av_misfit_residuals_l_${maxdeg}_set_${runner.cmbSets.indexOf(cmbSet)}.dat"))
      val allMisfits = for (res <- results.dropRight(burn)) yield runner.misfits(res._1)
      val arrangedMisfits = allMisfits.transpose
      val noResults = results.dropRight(burn).length.toDouble
      val avResiduals = for (resList <- arrangedMisfits) yield resList.reduceLeft(_ + _) / noResults
      for (printer <- p) {
        for (residual <- avResiduals(p.indexOf(printer))) {
          printer.print(s"$residual \n")
        }
        printer.close()
      }
    }

    printParametersAndNoise(new File(s"tomo_parameters_l_$maxdeg.dat"),
                            new File(s"topo_parameters_l_$maxdeg.dat"),
                            new File(s"noise_parameters_l_$maxdeg.dat"),
                            new File(s"variance_reduction_l_$maxdeg.dat"))
    printAvMisfits()
  }
}

object aicCalc {
  def main(args: Array[String]) {
    val l = args(0).toInt
    def xmlToListMaps(xmldata: Elem): List[Map[String, String]] = {
      def xmlNodeToMap(node: Node): Map[String, String] = {
        Map("type" -> (node \ "type").text, "otimes" -> (node \ "otimes").text,
            "raypaths" -> (node \ "raypaths").text, "topoin" -> (node \ "topoin").text)
      }
      (for (set <- xmldata \ "set") yield xmlNodeToMap(set)).toList
    }
    val inxml = XML.loadFile(args(1))
    val datamaps = xmlToListMaps(inxml)
    val qlen = (l + 1) * (l + 1)
    val initialGuess =  List((DenseVector.zeros[Double](qlen * 2), (for (set <- datamaps) yield 1.0).toList, Double.PositiveInfinity, 0.0, (for (set <- datamaps) yield 0.0).toList))
    val runner = new cmbHMC(1, 1.0, l, datamaps) // just dummy l, epsilon: we don't use them
    val mll = new DiffFunction[DenseVector[Double]] {
      def calculate(x: DenseVector[Double]) = {
        val sigs = x(0 to runner.cmbSets.length - 1).toArray.toList
        val q = x(runner.cmbSets.length to -1)
        (runner.mLogLike(q, sigs), runner.mLogLikeGrad(q, sigs))
      }
    }
    val lbfgs = new LBFGS[DenseVector[Double]](maxIter=40000, m=7)
    val k = runner.cmbSets(0).gMatrix.cols + runner.cmbSets.length
    val n = (for (set <- runner.cmbSets) yield set.gMatrix.rows).reduceLeft(_ + _)
    val initial = DenseVector.zeros[Double](k)
    initial(0 to runner.cmbSets.length - 1) := DenseVector(for (sigma <- args.drop(2)) yield sigma.toDouble)
    val minimized = lbfgs.minimize(mll, initial)
    val minMLogLike = mll(minimized)
    println(s"Sigma at min = ${minimized(0 to runner.cmbSets.length -1)}")
    println(s"- Max Log Likelihood = $minMLogLike")
    val aicc = 2.0 * minMLogLike + 2.0 * k.toDouble +
              (2.0 * k.toDouble * (k.toDouble + 1.0)) / (n.toDouble - k.toDouble - 1.0)
    println(s"AICc Estimate = $aicc")
    }
}
