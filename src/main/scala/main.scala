import scmb.hmc._

import breeze.linalg._

object ttSHA {
  def main(args: Array[String]) {
    val l  = args(0).toInt
    val ep = args(1).toDouble
    val maxIt = args(2).toInt
    val report = args(3).toInt
    val testfile = "/Users/jackmuir/Dropbox/Honours/Data/CMB/midpoints_H_M_C_PcP.final"
    val sha = new resSHA(l, ep, 8, List(Map("midpoints"->testfile)))
    sha.hmcRunner(81, maxIt, report, 1, 0, List((DenseVector.rand[Double](81),List(1.0))))
  }
}
