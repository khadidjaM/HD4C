package fr.inria.zenith.hd4c

import scala.collection.mutable.ArrayBuffer
import org.apache.commons.math3.distribution.MultivariateNormalDistribution
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator
import org.apache.commons.math3.analysis.interpolation.DividedDifferenceInterpolator
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator
import org.apache.commons.math3.distribution.BetaDistribution
import org.apache.commons.math3.distribution.GammaDistribution
import umontreal.ssj.rng.GenF2w32
import umontreal.ssj.stochprocess.OrnsteinUhlenbeckProcess
import org.apache.commons.math3.analysis.polynomials.PolynomialsUtils


class Worker(workerId : Int, val data : ArrayBuffer[Point], gamma : Double, var sigma : Double, var beta : Double, var sigma0 : Double, var beta0 : Double, var betaU : Double, tStar : Array[Double]) extends Serializable {
  
  
  
  /* if(workerId < 4)
      println("create worker : " + workerId) */
     
  //var Phis = new ArrayBuffer[Array[Double]]
  
// dimension
  val dim = data(0).vector.length
  
  val nbPoint : Int = data.size
  
  println("worker " + workerId + " : " + data(0).id + " - " + data(nbPoint-1).id)
  
  // assignements
  var c : Map[Int, Int] = Map()
  
  var clusterVector = new ArrayBuffer[Cluster]
  
  val m = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      m(i) = 0
    
  // Variance in each cluster
  
  
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
      
  var alpha : Double = 1
  
  val H : Int = 3
        
  def getAlpha = alpha
  def getC = c
  //def getPhis = Phis
  def getCluster(id : Int) = {
    
    val index = clusterVector.indexWhere(cluster => (cluster.id == id))
    clusterVector(index)
    
  }
  
  def addCluster(index : Int, c : Cluster) = {
    
    clusterVector.insert(index, c)
    
  }
  
/*   def removeCluster(id : Int) = {
    
    val index = clusterVector.indexWhere(cluster => (cluster.id == id))
    clusterVector.remove(index)
    
  } */
 
  def removePointFromCluster(clusterId : Int)={
    
    val index = clusterVector.indexWhere(cluster => (cluster.id == clusterId))
    
    clusterVector(index).size -= 1
    
  }
  
  def updateClusterId(id : Int, newId : Int) = {
    
    val index = clusterVector.indexWhere(cluster => (cluster.id == id))
    clusterVector(index).updateId(newId)
    
  }
  
/*   def Label(pointId : Int) = {
    
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
    
  } */
  
   def startWithClustering(globalClusters : ArrayBuffer[((Array[Double], Double), ArrayBuffer[(Int, Int)])], u : Double) = {
   
    betaU = u
     
    clusterVector.clear
    var newC : Map[Int, Int] = Map()
    
     var id = 0
    
    for (cluster <- globalClusters){
      
    val a = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      a(i) = 0
      
      val clust = new Cluster(id, 0, a, cluster._1._1, cluster._1._2, false)
      
      
      for(subC <- cluster._2.filter(_._1 == workerId))
        for(ci <- c.filter({case(k,v) => v == subC._2})){
          newC += ci._1 -> id
          clust.addPoint()
        }
      
      clusterVector += clust
      id+=1
      
    }
    
    c = newC
    
    //println("worker " + workerId + " c size : " + c.size + ", data size : " + data.size)
        
    this
        
  }
  
  def initializeWithClustering(initial : ArrayBuffer[(Array[Double], Double)], u : Double, s : Double, b :Double) = {
    
   /*  if(workerId < 4)
      println("initialize worker : " + workerId) */
     
  /*   println("initialize worker " + workerId + " ...") */
    
    beta = b
    beta0 = b
    sigma = s
    sigma0 = s
    
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
    
    betaU = u
    clusterVector.clear
    c = Map()
    
    var id = 0
    
    for (tuple <- initial){

    val a = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      a(i) = 0
      
      val clust = new Cluster(id, 0, a, tuple._1, tuple._2, false)
      clusterVector += clust
      id += 1
      
    }
    
    for(y <- data){
      
      //println("initialize data " + y.id + " ...")
      
      val clusterId : Int = max(calculateProbaWithLog(y))
      //val clusterId : Int = multinomiale(calculateProbaWithLog(y))
        
        c += (y.id -> clusterId)
        
        getCluster(clusterId).addPoint()
     
    }
    
    alpha = alphaInference(alpha, 1, 0.5, clusterVector.size, nbPoint)
    //println("worker " + workerId + " alpha : " + alpha)
    
   /*  println("worker " + workerId + " initilized")
    
    print("worker " + workerId + " nb clusters : " + clusterVector.size + " : ")
        for(c <- clusterVector)
          print(c.size + " ")
        println */
        
    this

  }
  
  def gibbsSampling(s : Double, b : Double, s0 : Double, b0 : Double)={
    
    /* println("worker " + workerId + " : Gibbs Sampling ...") */
    
    //Phis = new ArrayBuffer[Array[Double]]
    
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
    /*
    println
    println("Sigma0 : ")
    for(i <- 0 to dim - 1){
      for(j <- 0 to dim - 1)
        print(Sigma0(i)(j) + " ")
    println
    println
    
   
    }
    */
      
      
    for(iter <- 0 to 9){
      
    /*   println("Worker " + workerId + " itération : " + iter + ", nb clusters : " + clusterVector.size) */
            
      // for each data point
      for(y <- data){
        
        
        // remove it from its cluster
        removePointFromCluster(c(y.id))
        
        // the number of distinct cj for j != i  
        val K = clusterVector.size
          
        val yCluster = getCluster(c(y.id))
        // innovation (new clusters)
        //println("Worker " + workerId + " add auxiliary parameters ...")
        //
        /* if(yCluster.isNew && yCluster.size == 0 && clusterVector.last.id != yCluster.id){
            
            //Delete it
            removeCluster(yCluster.id)
            Label(y.id)
            addAuxiliaryParameters(H)
          }
          else
            // Cluster empty and it is the last
            if (yCluster.isNew && yCluster.size == 0)
              addAuxiliaryParameters(H-1)
          // Cluster not empty
          if(yCluster.size != 0 || (!yCluster.isNew && yCluster.size == 0 ))
            addAuxiliaryParameters(H) */
        //addAuxiliaryParameters(y.vector, H)
        addAuxiliaryParameters(H, y.vector(0))
        //println("Worker " + workerId + " add auxiliary parameters finished")
        // assign y 
        c += (y.id -> max(calculateProbaWithLog(y)))
        //c += (y.id -> multinomiale(calculateProbaWithLog(y)))
        
        if(c(y.id) >= K){
          val Beta = new BetaDistribution(1, gamma)
          val b = Beta.sample
          var cluster = getCluster(c(y.id))
          cluster.id = K
          cluster.size = 1
          cluster.beta = b * betaU
          c += (y.id -> K)
          addCluster(K, cluster)
          
          betaU = (1-b) * betaU
          
        }
        else
                   
          getCluster(c(y.id)).addPoint()
        
        for(j <- 0 to H - 1){
          
          val index = clusterVector.length - 1
          clusterVector.remove(index)
          
        }
        
        /* if(workerId == 0)
          println("Worker " + workerId + " data : " + y.id + ", nb clusters : " + clusterVector.size) */
        
      }
        
       /*  print("worker " + workerId + " nb clusters : " + clusterVector.size + " : ")
        for(c <- clusterVector)
          print(c.size + " ")
        println */
        alpha = alphaInference(alpha, 1, 0.5, clusterVector.size, nbPoint)
       /*  println("worker " + workerId + " alpha : " + alpha)
        
      print("worker " + workerId + " nb clusters : " + clusterVector.size + " : ")
        for(c <- clusterVector)
          print(c.size + " ")
        println */
     
    }
     
    //calculate means of clusters
    for(cluster <- clusterVector)    { 
        
      var dataOfCluster = new ArrayBuffer[Point]
      // data of cluster
      for(d <- data)
        if (c(d.id) == cluster.id)
          dataOfCluster += d
        
      cluster.calculateMean(dataOfCluster)

    }
    
    var result = new ArrayBuffer[ResumeCluster]
    
    for(cluster <- clusterVector){
      
      if(cluster.size != 0){
      
      val sub = new ResumeCluster(workerId, cluster.id, cluster.size, cluster.y_ , cluster.phi)
      result += sub
      }
      
    }
    
    println("worker " + workerId + " : Gibbs Sampling finished") 
    
    result

  }
  
  /* def addAuxiliaryParameters(H : Int) = {
    
    
    
    
    /* val N2 = new MultivariateNormalDistribution(m, Sigma0) */
    
    for (j <- 0 to H-1){

      //val auxiliaryPhi = N2.sample()
      
      val N = new NormalDistribution(0, Math.pow(sigma, 2)/(2 * beta))
    
      val x0 = N.sample

      val OU = new OrnsteinUhlenbeckProcess(x0, beta, 0, sigma, (new GenF2w32))

      OU.setObservationTimes(tStar, dim-1)
      
      val auxiliaryPhi = OU.generatePath
      
      Phis += auxiliaryPhi

      val id = clusterVector.length

      val a = new Array[Double](dim)
      for(i <- 0 to dim - 1)
        a(i) = 0

      val cluster = new Cluster(id, 0, a, auxiliaryPhi, betaU/H, true)

      clusterVector += cluster
    
    }
    
  } */
  
  /* def addAuxiliaryParameters(y : Array[Double], H : Int) = {
    
    
    
    
    /* val N2 = new MultivariateNormalDistribution(m, Sigma0) */
    
    for (j <- 0 to H-1){

      //val auxiliaryPhi = N2.sample()
      
      val N = new NormalDistribution(0, Math.pow(sigma, 2)/(2 * beta))
    
      val x0 = N.sample

      val OU = new OrnsteinUhlenbeckProcess(x0, beta, 0, sigma, (new GenF2w32))

      OU.setObservationTimes(tStar, dim-1)
      
      var auxiliaryPhi = OU.generatePath
      
      for(i <- 0 to dim - 1)
        auxiliaryPhi(i) += y(i)
      
      Phis += auxiliaryPhi

      val id = clusterVector.length

      val a = new Array[Double](dim)
      for(i <- 0 to dim - 1)
        a(i) = 0

      val cluster = new Cluster(id, 0, a, auxiliaryPhi, betaU/H, true)

      clusterVector += cluster
    
    }
    
  } */
  
  def addAuxiliaryParameters(H : Int, x0 : Double) = {
    
    for (j <- 0 to H-1){
      
      //val N = new NormalDistribution(0, Math.pow(sigma, 2)/(2 * beta))
    
      //val x0 = N.sample
      
      val OU = new OrnsteinUhlenbeckProcess(x0, beta, 0, sigma, (new GenF2w32))

      OU.setObservationTimes(tStar, dim-1)
      
      var auxiliaryPhi = OU.generatePath
      
      /* val r = scala.util.Random
      
      val P = PolynomialsUtils.createLegendrePolynomial(r.nextInt(10)) 
      //val P = PolynomialsUtils.createLegendrePolynomial(3)
      
      for(i <- 0 to dim - 1)
        auxiliaryPhi(i) += P.value(i/100) */
      
      val P1 = PolynomialsUtils.createLegendrePolynomial(1) 
      val P2 = PolynomialsUtils.createLegendrePolynomial(2) 
      val P3 = PolynomialsUtils.createLegendrePolynomial(3) 
      //val P4 = PolynomialsUtils.createLegendrePolynomial(4) 
      //val P5 = PolynomialsUtils.createLegendrePolynomial(5) 
      //val P = PolynomialsUtils.createLegendrePolynomial(3)
      
      val a = dirichletSample(Array(1, 1, 1))
      
      val step : Double = 2 / (dim - 1)
      
      var x : Double = -1
      
      for(i <- 0 to dim - 1){
        
        auxiliaryPhi(i) += a(0) * P1.value(x) + a(1) * P2.value(x) + a(2) * P3.value(x)
        x += step
        
      }
    /*   if(Phis.length < 1000)
        Phis += auxiliaryPhi */

      val id = clusterVector.length

      val y = new Array[Double](dim)
      for(i <- 0 to dim - 1)
        y(i) = 0

      val cluster = new Cluster(id, 0, y, auxiliaryPhi, betaU/H, true)

      clusterVector += cluster
    
    }
    
  }
  
  def calculateProbaWithLog(y : Point)={
    
    var proba : Map[Int, Double] = Map()
    for(cluster <- clusterVector){

        proba += (cluster.id -> (Math.log(cluster.size + alpha * cluster.beta) + logLikelihood(y.vector, cluster.phi, sigma, beta)- Math.log(nbPoint - 1 + alpha)))

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
  
  /* def calculateProbaWithLogI(y : Point)={
    
    var proba : Map[Int, Double] = Map()
    for(cluster <- clusterVector){

        proba += (cluster.id -> (Math.log(cluster.size + alpha * cluster.beta) + logLikelihoodI(y.vector, cluster.phi, sigma, beta)- Math.log(nbPoint - 1 + alpha)))

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
  } */
  
  /* def normalLikelihood(data : Array[Double], mean : Array[Double], sigma1 : Array[Array[Double]]) : Double ={
    var vraisemblance : Double = 1
    val N = new MultivariateNormalDistribution(mean, sigma1)
    N.density(data)
  
  } */
  
  def max(proba : Map[Int, Double])={
  
    val max = proba.maxBy(_._2)._1
    
    max
  }
 
  def alphaInference(alpha : Double, a : Double, b : Double, k : Int, n : Int) = {
    val Beta = new BetaDistribution(alpha + 1, n)
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
  
  def contingencyTable(nbGlobalClusters : Int, nbRealClusters : Int) = {
     
    var mat = Array.fill[Int](nbGlobalClusters, nbRealClusters)(0)        

    for(y <- data)
      mat(c(y.id))(y.realCluster - 1) += 1
    
    mat  
        
  }

  def dataCluster(clusterId : Int) = {
    
    var dataOfCluster = new ArrayBuffer[Point]
    
    for(d <- data)
      if (c(d.id) == clusterId)
        dataOfCluster += d
      
    dataOfCluster
    
  }
  
  
 
 def localRSS() = {
   
    var RSS = 0.0
   
    for(clust <- clusterVector){
      
      val dataOfCluster = dataCluster(clust.id)
      
      for(y <- dataOfCluster)
        for(i <- 0 to dim - 1)
          RSS += Math.pow(y.vector(i) - clust.phi(i), 2) / (sigma/10)
      
    }
    
    RSS
    
 }
 
  
 
  def scalarProduct(f : Array[Double], g : Array[Double], sigma : Double, beta : Double)={
    
   /*  if(workerId == 0){
    println
    print(" worker " + workerId + " f : ")
    for(i <- 0 to dim - 1)
      print(f(i) + " ")
     println
    print(" worker " + workerId + " g : ")
    for(i <- 0 to dim - 1)
      print(g(i) + " ")
    println
    } */
    
    val tDiff = tStar(1) - tStar(0)
    
    //val spline = new SplineInterpolator
    val spline = new LoessInterpolator

    
    val fg = f.zip(g).map{ case (a, b) => a * b }
    /* if(workerId == 0){
    println
     print(" worker " + workerId + " fg : ")
    for(i <- 0 to dim - 1)
      print(fg(i) + " ")
    println
    } */
    
    val fgSpline = spline.interpolate(tStar, fg)
    
   /*  if(workerId == 0)
      println(" worker " + workerId + " spline calculated ") */
    
  //  val fSpline = spline.interpolate(tStar, f)
    
  //  val gSpline = spline.interpolate(tStar, g)
    
   // val fDeriv = fSpline.derivative
  //  
   // val gDeriv = gSpline.derivative
    
    val fgDeriv = fgSpline.derivative
    
    //println("worker " + workerId + " : calcul integrale ...")
    
    val integr = new TrapezoidIntegrator
    
    val  A = tDiff * integr.integrate(org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT, fgSpline, fgSpline.getKnots()(0), fgSpline.getKnots().last)

    val B = tDiff * integr.integrate(org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT, fgDeriv, fgSpline.getKnots()(0), fgSpline.getKnots().last)
    
    //println("worker " + workerId + " : integrale calculé")
    
    (1 / Math.pow(sigma, 2)) * (B + Math.pow(beta, 2) * A) + (beta / Math.pow(sigma, 2)) * (f(0) * g(0) + f(tStar.length-1) * g(tStar.length-1))
    
  }
  
  /* def scalarProductI(f : Array[Double], g : Array[Double], sigma : Double, beta : Double)={
    
    val tDiff = tStar(1) - tStar(0)
    
    val spline = new SplineInterpolator
    
    val fg = f.zip(g).map{ case (a, b) => a * b }
    
    val fgSpline = spline.interpolate(tStar, fg)
    
    val fSpline = spline.interpolate(tStar, f)
    
    val gSpline = spline.interpolate(tStar, g)
    
    val fDeriv = fSpline.derivative
    
    val gDeriv = gSpline.derivative
    
    val fgDeriv = fgSpline.derivative
    
    val integr = new TrapezoidIntegrator
    
    /* val  A = tDiff * integr.integrate(org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT, fgSpline, tStar(0), tStar(tStar.length-1))

    val B = tDiff * integr.integrate(org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT, fgDeriv, tStar(0), tStar(tStar.length-1)) */
   
    val  A = tDiff * integr.integrate(org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT, fgSpline, fgSpline.getKnots()(0), fgSpline.getKnots().last)

    val B = tDiff * integr.integrate(org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT, fgDeriv, fgSpline.getKnots()(0), fgSpline.getKnots().last)
    
    (1 / Math.pow(sigma, 2)) * (B + Math.pow(beta, 2) * A) + (beta / Math.pow(sigma, 2)) * (f(0) * g(0) + f(tStar.length-1) * g(tStar.length-1))
    
  } */
  
  def logLikelihood(data : Array[Double], mean : Array[Double], sigma : Double, beta : Double) : Double ={
    
    scalarProduct(data, mean, sigma, beta) - scalarProduct(mean, mean, sigma, beta) / 2
    
  }
  
 /*  def logLikelihoodI(data : Array[Double], mean : Array[Double], sigma : Double, beta : Double) : Double ={
    
    scalarProductI(data, mean, sigma, beta) - scalarProduct(mean, mean, sigma, beta) / 2
    
  } */
  
  def parameters()={
    
    //if(workerId == 0)
      /* println("worker : " + workerId + " Parameters : ") */
    
    var B = 0.0
    var D = 0.0
    
    //if(workerId == 0)
   /*    println("clusterVectorSize : " + clusterVector.length) */
            
    for(z <- clusterVector){
      
      val size = z.size.toInt
      
      if(size > 0){
      
     // if(workerId == 0)
    /*   println("z Size : " + size) */
      
      var x0 = new Array[Double](size)
      var x1 = new Array[Double](size)
      
      val dataZ = dataCluster(z.id)
      
      for(i <- 0 to size - 1){
        
        x0(i) = dataZ(i).vector(0)
        x1(i) = dataZ(i).vector(1)
        
      }
    
      val moy0 = x0.sum / size
      val moy1 = x1.sum / size
      
      //if(workerId == 0)
/*       println("moy0 : " + moy0 + " moy1 : " + moy1) */
      
      for(i <- 0 to size - 1){
        
        x0(i) -= moy0
        x1(i) -= moy1
        
      }
      
      
      /* var sum0 = 0.0
      for(i <- 0 to size -1)
        sum0 += x0(i) - moy0
        //sum0 += Math.abs(x0(i) - moy0)
      
      var sum1 = 0.0
      for(i <- 0 to size -1)
        sum1 += x1(i) - moy1
        //sum1 += Math.abs(x1(i) - moy1)
      
      var sum0sqr = 0.0
      for(i <- 0 to size -1)
        sum0sqr += Math.pow(x0(i) - moy0, 2)
      
      var sum1sqr = 0.0
      for(i <- 0 to size -1)
        sum1sqr += Math.pow(x1(i) - moy1, 2) */
       
      for(i <- 0 to size - 1){
      
        B += x0(i) * x1(i)
        D += Math.pow(x0(i), 2) + Math.pow(x1(i), 2)
      
      }
      
    }
      
    }
    
    //if(workerId == 0)
   /*    println("B : " + B + " D : " + D) */
    
    (B, D, clusterVector.length)
    
  } 
  
  /* def multinomiale(proba : Map[Int, Double])={
   // println("multinomiale")
   
  var k : Int = 0
  
    /*
    var max : Int = 0
    for(k <- 0 to proba.length - 1)
      if (proba(k) > proba(max))
        max = k
    */
   
   // val max = proba.maxBy(_._2)._1
    
    proba.foldLeft(0.0){ case (a, (k, v)) => a + v } 
    val r = scala.util.Random
    val v = r.nextDouble
    println("v : " + v)
    var s : Double = 0
    do{
      s += proba(k)
      println("s : " + s)
      k += 1
    }while(k <= proba.size && s <= v)
    
   //println("max : " + max)
    k - 1
    
 //max
  } */
  
   def dirichletSample(alpha : Array[Double]) = {
    
    var x = new Array[Double](alpha.length)
    var y = new Array[Double](alpha.length)
    
    for(i <- 0 to alpha.length - 1){
      
      val gamma = new GammaDistribution(alpha(i), 1)
      y(i) = gamma.sample
      
    }
    
    val sum = y.foldLeft(0.0){ case (a, b) => a + b } 
    
    for(i <- 0 to alpha.length - 1)
      x(i) = y(i)/sum
    
    x
      
  }
 
  def firstParameters()={
     
    /* println("worker " + workerId + " : Estimate First Parameters ...") */
    
    var B = 0.0
    var D = 0.0
            
    //for(d <- data){
      
      var x0 = new Array[Double](nbPoint)
      var x1 = new Array[Double](nbPoint)
      
      for(i <- 0 to nbPoint - 1){
        
        x0(i) = data(i).vector(0)
        x1(i) = data(i).vector(1)
        
      }
    
      val moy0 = x0.sum / nbPoint
      val moy1 = x1.sum / nbPoint
      
      for(i <- 0 to nbPoint - 1){
        
        x0(i) -= moy0
        x1(i) -= moy1
        
      }
       
      for(i <- 0 to nbPoint - 1){
      
        B += x0(i) * x1(i)
        D += Math.pow(x0(i), 2) + Math.pow(x1(i), 2)
      
      }
      
    //}
    
   /*  println("worker " + workerId + " : First Parameters estimated") */
    
    (B, D, 1)
    
  }
  
}
