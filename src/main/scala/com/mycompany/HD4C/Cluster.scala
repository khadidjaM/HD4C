package fr.inria.zenith.hd4c

import scala.collection.mutable.ArrayBuffer
import Jama.Matrix
import umontreal.ssj.rng.GenF2w32
import umontreal.ssj.stochprocess.OrnsteinUhlenbeckProcess

abstract class AbstractCluster(var id : Int, var size : Double, var y_ : Array[Double], var phi : Array[Double]) extends Serializable{
  
  val dim = y_.length
  
  def empty() = {size == 0}
  
  def calculateMean()={}
  
  def updateId(newId : Int) = {
    
      id = newId
      
  }
 
}

class Cluster(id : Int, s : Double, y : Array[Double], p : Array[Double], var beta : Double, var newCluster : Boolean) extends AbstractCluster(id, s, y , p){
 
  def isNew() = newCluster
    
  def addPoint() = {
      
    size += 1
      
  }
    
  def removePoint(pId : Int) = {
      
    size -= 1
      
  }
  
  def calculateMean(dataCluster : ArrayBuffer[Point])={
    
    if (empty)
      
      phi = Array.fill(dim)(0.0)
    
    else{
      
      val m = new Array[Double](dim)

      for(j <- 0 to dim - 1){

        var s : Double = 0

        for(d <- dataCluster)
            s = s + d.vector(j)

        m(j) = s / dataCluster.length

      }

      y_ = m
    }
    
  }
    
}

class GlobalCluster(var subIds : ArrayBuffer[(Int, Int)], id : Int, size : Double, y_ : Array[Double], p : Array[Double], var mean : Array[Double], var SigmaC : Array[Array[Double]])extends AbstractCluster(id, size, y_ , p){
  
  def updateY_(means : ArrayBuffer[_<:AbstractCluster])={
    
    var sum = Array.fill(dim)(0.0)

    var sumSize : Double = 0
    for (mean <- means){
      for(i <- 0 to dim - 1)
        sum(i) += mean.size * mean.y_(i)
      sumSize += mean.size
    }
    
    for(i <- 0 to dim - 1)
      y_(i) = sum(i) / sumSize
    
  }
  
  def updatePhi(mu : Array[Double], Sigma0 : Array[Array[Double]], Sigma : Array[Array[Double]], sigma : Double, beta : Double, tStar : Array[Double]) = {
   
  //val T = postSigma(Sigma0, Sigma)
  
  val S = postMean(/* /* mu,  */ */Sigma0, Sigma)
  
  /*  println("post mean : ")
  for(i <- 0 to dim - 1)
    print(S(i) + " ")
  println
    
  println("post sigma0 : ")
  for(i <- 0 to dim - 1){
    for(j <- 0 to dim - 1)
      print(T(i)(j) + " ")
  println
  }  */
   
 /*  val N = new MultivariateNormalDistribution(S, T)
  
    phi = N.sample */
   
    //val N = new NormalDistribution(0, Math.pow(sigma, 2)/(2 * beta * nc))
    
    //val x0 = N.sample
    
    val x0 = S(0)

    val OU = new OrnsteinUhlenbeckProcess(x0, beta, 0, sigma/Math.sqrt(size), (new GenF2w32))

    OU.setObservationTimes(tStar, dim-1)
      
    phi = OU.generatePath
    
    phi(0) = S(0)
      
    for(i <- 1 to dim - 1)
      phi(i) += S(i)
      
  }
  
  def postMean(/* mu : Array[Double],  */Sigma0 : Array[Array[Double]], Sigma : Array[Array[Double]]) = {
    
    val y = Array.ofDim[Double](dim, 1)
    for(i <- 0 to dim - 1)
      y(i)(0) = y_(i)
      
    val m = Array.ofDim[Double](dim, 1)
    for(i <- 0 to dim - 1)
      /* m(i)(0) = mu(i) */
    m(i)(0) = mean(i)
    
    var sigma = Array.ofDim[Double](dim, dim)
    
    for (i <- 0 to dim - 1)
      for (j <- 0 to dim - 1)
        sigma(i)(j) = Sigma(i)(j) / size
    
    val s = matrixSum(m, matrixProduct(transposeMatrix(Sigma0), matrixProduct(inverseMatrix(matrixSum(sigma, Sigma0)), matrixSub(y, m))))  
    
    var meanPost = new Array[Double](dim)
    
    for (i <- 0 to dim - 1)
      mean(i) = s(i)(0)
    
    for (i <- 0 to dim - 1)
      meanPost(i) = s(i)(0)
    
    meanPost
    
  }
    
  def postSigma(Sigma0 : Array[Array[Double]], Sigma : Array[Array[Double]]) = {
    
    var sigma = Array.ofDim[Double](dim, dim)
    
    for (i <- 0 to dim - 1)
      for (j <- 0 to dim - 1)
        sigma(i)(j) = Sigma(i)(j) / size
    
    val s = matrixSub(Sigma0, matrixProduct(transposeMatrix(Sigma0), matrixProduct(inverseMatrix(matrixSum(Sigma0, sigma)), Sigma0)))
    
     for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1)
        if(s(i)(j) < Math.pow(10, -5))
          s(i)(j) = 0
     
    s
    
  }
  
  def matrixProduct(A : Array[Array[Double]], B : Array[Array[Double]]) : Array[Array[Double]] = {
    
    val m = A.length
    val n = B.length
    val p = B(0).length
    var C = Array.ofDim[Double](m, p)
    
    for(i <- 0 to m -1)
      for(j <- 0 to p - 1)
        {
          var som : Double = 0
          
          for(k <- 0 to n - 1)
            som += A(i)(k) * B(k)(j)
          
          C(i)(j) = som
        }
    C
  }

  def matrixSum(A : Array[Array[Double]], B : Array[Array[Double]]) : Array[Array[Double]] = {
    
    val n = A.length
    val m = A(0).length
    var C = Array.ofDim[Double](n, m)

    for(i <- 0 to n - 1)
      for(j <- 0 to m - 1)
        C(i)(j) = A(i)(j) + B(i)(j)
        
    C 
  }
  
  def matrixSub(A : Array[Array[Double]], B : Array[Array[Double]]) : Array[Array[Double]] = {
    
    val n = A.length
    val m = A(0).length
    var C = Array.ofDim[Double](n, m)

    for(i <- 0 to n - 1)
      for(j <- 0 to m - 1)
        C(i)(j) = A(i)(j) - B(i)(j)
        
    C 
  }
  
  def inverseMatrix (A : Array[Array[Double]]) = {
    
    val M = new Matrix(A)
    
    M.inverse.getArray
    
  }
  
  def transposeMatrix (A : Array[Array[Double]]) = {
    
    val M = new Matrix(A)
    
    M.transpose.getArray
    
  }
  
  
}

class ResumeCluster(val workerId : Int, id : Int, size : Double, y_ : Array[Double], phi : Array[Double]) extends AbstractCluster(id, size, y_ , phi ){}
