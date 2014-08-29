import breeze.linalg._

class HMC {

  // For now calcU and calcUgrad return 0's, not sure what a good default is
  def calcU(q: DenseVector[Double]): Double = 0

  def calcUgrad(q: DenseVector[Double]): DenseVector[Double] = DenseVector.zeros[Double](q.length)

  def hmcupdate(q: DenseVector[Double], acc: Int): (DenseVector[Double], Int) = (q, acc)

  val hmcu = (hmcupdate _).tupled
}
