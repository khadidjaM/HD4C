package fr.inria.zenith.hd4c

import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator
import org.apache.commons.math3.analysis.interpolation.DividedDifferenceInterpolator
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator
import org.apache.commons.math3.distribution.BetaDistribution
import org.apache.commons.math3.distribution.GammaDistribution
import scala.collection.mutable.ArrayBuffer

class Master(var sigma : Double, var beta : Double, var sigma0 : Double, var beta0 : Double, n : Int, var gamma : Double, parameters : ArrayBuffer[(Array[Double], Array[Array[Double]])], dim : Int, tStar : Array[Double]) {
  
  val H : Int = 1
  
  var c : Map[Int, Int] = Map()
  
  var clusterVector = new ArrayBuffer[GlobalCluster]
  
  var individuals = new ArrayBuffer[GlobalCluster]
  
  val m = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      m(i) = 0
  
  var Sigma = Array.ofDim[Double](dim, dim)
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1){
        val s = (Math.pow(sigma, 2) / (2 * beta)) * Math.exp(- beta * Math.abs(tStar(i)-tStar(j))) 
        if (s < Math.pow(10, -3) )
          Sigma(i)(j) = 0
        else
          Sigma(i)(j) = s
      }
      
  // Variance between centers
  
  var Sigma0 = Array.ofDim[Double](dim, dim)
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1){
        val s = (Math.pow(sigma0, 2) / (2 * beta0)) * Math.exp(- beta0 * Math.pow(tStar(i)-tStar(j), 2))
        if (s < Math.pow(10, -3) )
          Sigma0(i)(j) = 0
        else
          Sigma0(i)(j) = s
      }
 

  def createNewCluster(y : GlobalCluster) = {
    
    val id = clusterVector.size
    
    /*
    val m = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      m(i) = y.mean(i) 
    
    val yi = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      yi(i) = y.y_(i)
    
    val f = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      f(i) = y.phi(i)
    
    val s = Array.ofDim[Double](dim, dim)
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1)
        s(i)(j) = y.SigmaC(i)(j)
    
    val siz = y.size
    
    */

    var superCluster = new GlobalCluster(new ArrayBuffer[(Int, Int)], id, y.size, y.y_, y.phi, y.mean, y.SigmaC)
    
    superCluster.updatePhi(m, Sigma0, Sigma, sigma, beta, tStar)
    
    clusterVector += superCluster
    c += y.id -> id
    
  }
  
  def getCluster(id : Int) = {
    
    val index = clusterVector.indexWhere(cluster => (cluster.id == id))
    clusterVector(index)
    
  }
  
  def addCluster(index : Int, c : GlobalCluster) = {
    
    clusterVector.insert(index, c)
    
  }
  
  def removeCluster(id : Int) = {
    
   // val index = clusterVector.indexWhere(cluster => (cluster.id == id))
    clusterVector.remove(id)
    
  }
  
  def updateClusterId(id : Int, newId : Int) = {
    
    val index = clusterVector.indexWhere(cluster => (cluster.id == id))
    clusterVector(index).updateId(newId)
    
  }
  
  def removePointFromCluster(pId : Int, clusterId : Int)={
    
    val index2 = clusterVector.indexWhere(cluster => (cluster.id == clusterId))
    clusterVector(index2).size -= 1
 
  }
  
  def getGamma() = gamma
  
  
  def initialize() = {
    
    clusterVector = new ArrayBuffer[GlobalCluster]
    
   /*  println("nb individuals : " + individuals.length)
    println */
    
    for(y <- individuals){
      createNewCluster(y)
    }
    
  }
  
  def Label(pointId : Int) = {
    
    val clusterId = c(pointId)
    
    //Update Clusters
    
    for(cluster <- clusterVector)
      
      if (cluster.id > clusterId)
        
        updateClusterId(cluster.id, cluster.id - 1)
    
    //Update C
    
    for(j <- (c - pointId).keys)
      
      if (c(j) > c(pointId)){
        
        val temp = c(j) - 1
        c += (j -> temp)
        
      }
    
  }

  def gibbsSampling(s : Double, b : Double, s0 : Double, b0 : Double)={
    
    //println("Gibbs ...")

    beta = b
    beta0 = b0
    sigma = s
    sigma0 = s0
    
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1){
        val s = (Math.pow(sigma, 2) / (2 * beta)) * Math.exp(- beta * Math.abs(tStar(i)-tStar(j))) 
        if (s < Math.pow(10, -3) )
          Sigma(i)(j) = 0
        else
          Sigma(i)(j) = s
      }
      
    // Variance between centers

    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1){
        val s = (Math.pow(sigma0, 2) / (2 * beta0)) * Math.exp(- beta0 * Math.pow(tStar(i)-tStar(j), 2))
        if (s < Math.pow(10, -3) )
          Sigma0(i)(j) = 0
        else
          Sigma0(i)(j) = s
      }
    
   /*  println
    println("Sigma0 : ")
    for(i <- 0 to dim - 1){
      for(j <- 0 to dim - 1)
        print(Sigma0(i)(j) + " ")
    println
    }
    
    println
    println("Sigma : ")
    for(i <- 0 to dim - 1){
      for(j <- 0 to dim - 1)
        print(Sigma(i)(j) + " ")
    println
    }
    
    println */
   
    var indice = 0

    for(iter <- 0 to 2){
      
      //println("Master iteration : " + iter)
    
      // for each data point
      for(y <- individuals){
          
          val oldC = c(y.id)
          
        var retire = false
          
          if(getCluster(oldC).size == y.size){
            
            //Delete it
            removeCluster(oldC)
            Label(y.id)
            retire = true
          }
         else 
        removePointFromCluster(y.id, c(y.id))
          
          val K = clusterVector.size
          // innovation          
          addAuxiliaryParameters(H, y)

        c += (y.id -> max(calculateProbaWithLog(y)))        
         
        if(c(y.id) >= K){
     
          var cluster = getCluster(c(y.id))
          cluster.id = K
          cluster.size = y.size
          c += (y.id -> K)
          addCluster(K, cluster)
     
        }
        else 
          getCluster(c(y.id)).size += y.size
        
        if(! retire){
          val effective = y.size
          val index2 = clusterVector.indexWhere(cluster => (cluster.id == oldC))
          clusterVector(index2).size -= effective - 1
         }
        
        for(j <- 0 to H - 1){
          
          val index = clusterVector.length - 1
          clusterVector.remove(index)
          
        }
   
      }
      
   //   println("update means")
     
      for(cluster <- clusterVector)    { 
        
        var dataOfCluster = new ArrayBuffer[GlobalCluster]
        // data of cluster
        for(y <- individuals)
          if (c(y.id) == cluster.id)
            dataOfCluster += y
 
       cluster.updateY_(dataOfCluster)
        

      }
      
      /* println("nb Clusters : " + clusterVector.length)
      println */
      
      // update phis
      
      for(cluster <- clusterVector)    { 
        
        cluster.updatePhi(m, Sigma0, Sigma, sigma, beta, tStar)
        
      }
      
      
      gamma = gammaInference(gamma, 1, 0.5, clusterVector.size, n)
      //println("Master Gamma : " + gamma)
    /*    println("Itération " + iter + " finie, nb clusters : " + clusterVector.length) */
      
      
  }
    
    var count = 0
    var count2 = 0
    
    for(cluster <- clusterVector)    { 
      
        // data of cluster
        for(y <- individuals.filter({case(p) => c(p.id) == cluster.id})){
          for(sub <- y.subIds){
            cluster.subIds += sub
            count += 1
          }
        }
        
      count2 += cluster.subIds.size

      }
    
    clusterVector
    
  }
  
def addAuxiliaryParameters(H : Int, y : GlobalCluster) = {

    val m = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      m(i) = 0
    
    val a = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      a(i) = 0
    
    for (j <- 0 to H-1){

    val id = clusterVector.length
    
    val superCluster = new GlobalCluster(new ArrayBuffer[(Int, Int)], id, gamma, a, y.phi, m, Sigma0)
    
    clusterVector += superCluster
    
    }
    
  }
  
    def calculateProbaWithLog(y : GlobalCluster)={
    
    var proba : Map[Int, Double] = Map()
    for(cluster <- clusterVector){

        proba += (cluster.id -> (Math.log(cluster.size) + logLikelihood(y.y_, cluster.phi, sigma, beta)- Math.log(n - y.size + gamma)))

    }
  
  val max = proba.values.max
   
    for ((cluster, p) <- proba){
      proba += cluster -> Math.exp(p-max)
  
    }
    
    val sum = proba.foldLeft(0.0){ case (a, (k, v)) => a + v } 
    
    for ((cluster, p) <- proba){
      proba += (cluster -> p/sum)
  
    }

    proba
  }
  
  /* def normalLikelihood(data : Array[Double], mean : Array[Double], sigma1 : Array[Array[Double]]) : Double ={
    var vraisemblance : Double = 1
    val N = new MultivariateNormalDistribution(mean, sigma1)
    N.density(data)
 
  } */

  def max(proba : Map[Int, Double])={
   
    proba.maxBy(_._2)._1
    
  }
  
  def updateExistingClasses(ExistingK : Int, cent : ArrayBuffer[ResumeCluster]) = {
    
    var centers = new ArrayBuffer[ResumeCluster]
    for(c <- cent)
      centers += c
   
    individuals.clear  
    
    var id = 0
    
    if(ExistingK != 0)
      for(i <- 0 to ExistingK-1){

        var mergedClusters = new ArrayBuffer[ResumeCluster]
        var subIds = new ArrayBuffer[(Int, Int)]

        for(y <- centers){

          if(y.id == i){

            mergedClusters += y
            subIds += y.workerId -> i

          }

        }

        if(!mergedClusters.isEmpty){

          for(y <- mergedClusters)
            centers -= y

    
          val a = new Array[Double](dim)
          for(i <- 0 to dim - 1)
            a(i) = 0

          val sum = mergedClusters.foldLeft(0.0){ case (a, b) => a + b.size} 
          var superCluster = new GlobalCluster(subIds, id, sum, a, mergedClusters(0).phi, parameters(mergedClusters(0).id)._1, parameters(mergedClusters(0).id)._2)
          superCluster.updateY_(mergedClusters)
          individuals += superCluster

          id += 1
        }

      }
    
    for(y <- centers){

      var subIds = new ArrayBuffer[(Int, Int)]
      subIds += y.workerId -> y.id
      val superCluster = new GlobalCluster(subIds, id, y.size, y.y_ , y.phi, m, Sigma0)
      individuals += superCluster
      
      id += 1
      
    }
    
  }
 
    def gammaInference(gamma : Double, a : Double, b : Double, k : Int, n : Int) = {
    val Beta = new BetaDistribution(gamma + 1, n)
    val eta = Beta.sample
    val pi = (a + k - 1)/(a + k - 1 + n * (b - Math.log(eta)))  
    val rand = scala.util.Random
    val u = rand.nextDouble
   
    if (u <= pi){
      val Gamma = new GammaDistribution(a + k, b - Math.log(eta))
      Gamma.sample
    }else {
      val Gamma = new GammaDistribution(a + k - 1, b - Math.log(eta))
      Gamma.sample
    }    
  }
  
  def scalarProduct(f : Array[Double], g : Array[Double], sigma : Double, beta : Double)={
    
    val tDiff = tStar(1) - tStar(0)
    
    //val spline = new SplineInterpolator
    
    val spline = new LoessInterpolator
    
    val fg = f.zip(g).map{ case (a, b) => a * b }
    
    val fgSpline = spline.interpolate(tStar, fg)
    
    val fSpline = spline.interpolate(tStar, f)
    
    val gSpline = spline.interpolate(tStar, g)
    
    val fDeriv = fSpline.derivative
    
    val gDeriv = gSpline.derivative
    
    val fgDeriv = fgSpline.derivative
    
    val integr = new TrapezoidIntegrator
    
    val A = tDiff * integr.integrate(org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT, fgSpline, tStar(0), tStar(tStar.length-1))

    val B = tDiff * integr.integrate(org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT, fgDeriv, tStar(0), tStar(tStar.length-1))
    
    (1 / Math.pow(sigma, 2)) * (B + Math.pow(beta, 2) * A) + beta / Math.pow(sigma, 2) * (f(0) * g(0) + f(tStar.length-1) * g(tStar.length-1))
    
  }
  
  def logLikelihood(data : Array[Double], mean : Array[Double], sigma : Double, beta : Double) : Double ={
    
    scalarProduct(data, mean, sigma, beta) - scalarProduct(mean, mean, sigma, beta) / 2
    
  }
}
