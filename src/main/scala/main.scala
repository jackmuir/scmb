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

import scmb.cmbhmc.cmbHMC

import breeze.linalg.DenseVector

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
    val initialGuess =  List((DenseVector.zeros[Double](qlen * 2), (for (set <- datamaps) yield 1.0).toList))
    val runner = new cmbHMC(l, ep, maxdeg, datamaps)
    val results = runner.hmcRunner(maxIt, report, 1, 0, initialGuess)

    def printParametersAndNoise(f1: File, f2: File, f3: File) {
      val p1 = new PrintWriter(f1)
      val p2 = new PrintWriter(f2)
      val p3 = new PrintWriter(f3)
      for (res <- results.dropRight(burn)) {
        for (tomo <- res._1.apply(0 to qlen - 1)) p1.print(s"$tomo ")
        for (topo <- res._1.apply(qlen to 2 * qlen - 1)) p2.print(s"$topo ")
        for (noise <- res._2) p3.print(s"$noise ")
        p1.print("\n"); p2.print("\n"); p3.print("\n")
      }
      p1.close(); p2.close(); p3.close()
    }

    printParametersAndNoise(new File("tomo_parameters.dat"),
                            new File("topo_parameters.dat"),
                            new File("noise_parameters.dat"))


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
}*/
