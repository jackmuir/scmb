import breeze.linalg._

import breeze.stats.distributions.Gamma

def stdDevSample(phi: Double, n: Int) = math.sqrt(
  phi / (2.0 * Gamma((n.toDouble + 3) / 2.0, 1.0).sample())
  )

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
