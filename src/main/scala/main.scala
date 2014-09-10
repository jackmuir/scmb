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
    val initialGuess =  List((DenseVector.zeros[Double](qlen * 2), (for (set <- datamaps) yield 1.0).toList, Double.PositiveInfinity))
    val runner = new cmbHMC(l, ep, maxdeg, datamaps)
    val results = runner.hmcRunner(maxIt, report, 1, 0, initialGuess)

    def printParametersAndNoise(f1: File, f2: File, f3: File) {
      val p1 = new PrintWriter(f1)
      val p2 = new PrintWriter(f2)
      val p3 = new PrintWriter(f3)
      for (res <- results.dropRight(burn)) {
        for (tomo <- res._1.apply(0 to qlen - 1)) p1.print(s"$tomo ")
        ////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!////////!!!!!!!!
        for (topo <- res._1.apply(qlen to 2 * qlen - 1)) p2.print(s"${topo * 1000.0} ") // see shdatasets for 1000
        for (noise <- res._2) p3.print(s"$noise ")
        p1.print("\n"); p2.print("\n"); p3.print("\n")
      }
      p1.close(); p2.close(); p3.close()
    }

    printParametersAndNoise(new File(s"tomo_parameters_l_$maxdeg.dat"),
                            new File(s"topo_parameters_l_$maxdeg.dat"),
                            new File(s"noise_parameters_l_$maxdeg.dat"))
    val numData = (for (set <- runner.cmbSets) yield set.gMatrix.rows).reduceLeft(_ + _)
    val minmLL = min(for (res <- results.dropRight(burn)) yield res._3)
    val k = runner.cmbSets(0).gMatrix.cols.toDouble//qlen.toDouble
    val aicc = 2.0 * minmLL + 2.0 * k + (2.0 * k * (k + 1.0)) / (numData.toDouble - k - 1.0)
    println(s"Max Log Likelihood = ${-minmLL}")
    println(s"AICC = $aicc")
  }
}

/*object ttSHA {
  def main(args: Array[String]) {
    val maxl = args(0).toInt
    val l  = args(1).toInt
    val ep = args(2).toDouble
    val maxIt = args(3).toInt
    val burn = args(4).toInt
    val report = args(5).toInt
    val testfile = "/Users/jackmuir/Dropbox/Honours/Data/CMB/midpoints_H_M_C_PcP.final"
    val datamap = List(Map("midpoints"->testfile))
    val sha = new resSHA(l, ep, maxl, datamap)
    val results = sha.hmcRunner(81, maxIt, report, 1, 0, List((DenseVector.rand[Double](81),List(1.0))))

    //Inspiration: Rex Kerr - http://stackoverflow.com/questions/4604237/how-to-write-to-a-file-in-scala
    def printParametersAndNoise(f1: java.io.File, f2: java.io.File) {
      val p1 = new java.io.PrintWriter(f1)
      val p2 = new java.io.PrintWriter(f2)
      for (res <- results.dropRight(burn)) {
        for (param <- res._1) {p1.print(param); p1.print(" ")}
        for (noise <- res._2) {p2.print(noise); p2.print(" ")}
        p1.print("\n"); p2.print("\n")
      }
      p1.close(); p2.close()
    }

    printParametersAndNoise(new File("sha_parameters.dat"), new File("noise_parameters.dat"))


  }
}

val maxdeg = 5
val qlen = (maxdeg + 1) * (maxdeg + 1)
val l  = 1
val ep = 1.0
val maxIt = 1
val burn = 0
val report = 1
def xmlToListMaps(xmldata: Elem): List[Map[String, String]] = {
  def xmlNodeToMap(node: Node): Map[String, String] = {
    Map("type" -> (node \ "type").text, "otimes" -> (node \ "otimes").text,
        "raypaths" -> (node \ "raypaths").text, "topoin" -> (node \ "topoin").text)
  }
  (for (set <- xmldata \ "set") yield xmlNodeToMap(set)).toList
}
val inxml = XML.loadFile("datasets.xml")
val datamaps = xmlToListMaps(inxml)
val initialGuess =  List((DenseVector.zeros[Double](qlen * 2), (for (set <- datamaps) yield 1.0).toList))
val runner = new cmbHMC(l, ep, maxdeg, datamaps)

*/
