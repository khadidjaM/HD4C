/*
 * Copyright 2019 Khadidja Meguelati <khadidja.meguelati@inria.fr>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package fr.inria.zenith.hd4c

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf 
import scala.collection.mutable.ArrayBuffer
//import java.io.File
//import java.io.PrintWriter
import org.apache.commons.math3.distribution.GammaDistribution
import org.apache.spark.mllib.clustering.KMeans
import org.apache.spark.mllib.linalg.Vectors


object App { 
  
  def convertIteratorToArrayBuffer(data : Iterator[Point]) = {
    
    var dataBuffer = new ArrayBuffer[Point]
    var index = 0
    
    while (data.hasNext){
      
      dataBuffer.insert(index, data.next)
      index += 1
      
    }
      
    dataBuffer
    
  }
  
  def arrayToArrayBuffer(A : ArrayBuffer[ResumeCluster])= {
    
    var B = new ArrayBuffer[ResumeCluster]
    
    for (i <- 0 to A.length - 1)
      if(A(i).size > 0){
        B += A(i)
      
      }
    
    B
     
  }
 
  def toPoint(d : Array[Double], i : Int, r : Int)={
    
    val p = new Point(i, d, r)
    
    p
  }
  
  //simuler un vecteur selon une loi de Dirichlet 
  
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
 
  def lineData(line : String, dim : Int, withReal : Int) = {
    
    val d = new Array[Double](dim)
    val temp = line.split(" ")
    
    for (j <- 0 to dim - 1)
      d(j) = temp(j).toDouble
    if(withReal == 1 )
      (d, temp(dim).toDouble.toInt)
    else
      (d, 0)
    
  }
 
  def unionBuffers( A: Array[ArrayBuffer[ResumeCluster]]) = {
      
    var B = new ArrayBuffer[ResumeCluster]
      
    for(a <- A)
      B ++= a
        
    B
      
   }
  
  def main(args: Array[String]):Unit={
    
    if(args.length < 7) {
            System.err.println( "Use: HD4C <dimensions> <number of workers> <number of distributions> <target to data file> <number of clusters for Kmeans> <number of real clusters> <1 if the real number of clusters is known>")
            System.exit(1)
    }
    
      
    val conf = new SparkConf().setAppName("ParallelDP")
    val sc = new SparkContext(conf) 
    
    val dim = args(0).toInt
    val nbWorkers = args(1).toInt
    val nbDist = args(2).toInt
    val data = sc.textFile(args(3), nbWorkers).persist
    
    var numClusters = args(4).toInt
    val numRealClusters = args(5).toInt
    val withReal = args(6).toInt
    
    //val mult = args(7).toDouble
    val mult = 10.0
    
    /* val writer = new PrintWriter(new File(args(3) + "C"))
    val writerPhis = new PrintWriter(new File(args(3) + "Phis"))
    val writerResult = new PrintWriter(new File(args(3) + "Result")) */
    
    var tStar = new Array[Double](dim)
    
    for(i <- 0 to dim -1)
      tStar(i) = i
    
    val n = data.count.toInt

    val parsedData = data.map(s => Vectors.dense(s.split(' ').take(dim).map(_.toDouble)))
    
    var gamma : Double = 1
    // intialize by Kmeans
    val model = KMeans.train(parsedData, numClusters, 10)
    
    val centers = model.clusterCenters
    
    val clusterInd = model.predict(parsedData)

    val clusterSizes = clusterInd.countByValue.values
     
    parsedData.unpersist()
     
    val dataRDD = data.zipWithIndex().map{ case (line, i) => {
          val p = lineData(line, dim, withReal) 
          toPoint(p._1, i.toInt, p._2)}}.persist
   
    val m = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      m(i) = 0

    // create workers
    
/*     println
    println(" *********** create workers ************* ")
    println */
    
    var workers = dataRDD.mapPartitionsWithIndex{(index, dataPartition) => {Iterator(new Worker(index, convertIteratorToArrayBuffer(dataPartition), gamma, 2.5, 15, 2.5, 15, 1, tStar))}}.persist
    
    var time : Long = 0  
    var beginTime = System.currentTimeMillis()
    
    // Estimate parameters
    var param = workers.mapPartitions{ w => {Iterator(w.next.firstParameters)}}.collect
    
    var distParam = DistributionsParameters(param, n)
    
    var beta = distParam._1
    var beta0 = distParam._1
    var sigma = distParam._2
    var sigma0 = distParam._2
    
    var Sigma0 = Array.ofDim[Double](dim, dim)
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1){
        val s = (Math.pow(sigma0, 2) / (2 * beta0)) * Math.exp(- beta0 * Math.pow(tStar(i)-tStar(j), 2))
        if (s < Math.pow(10, -3) )
          Sigma0(i)(j) = 0
        else
          Sigma0(i)(j) = s
      }
  
    var Sigma = Array.ofDim[Double](dim, dim)
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1){
        val s = (Math.pow(sigma, 2) / (2 * beta)) * Math.exp(- beta * Math.abs(tStar(i)-tStar(j))) 
        if (s < Math.pow(10, -3) )
          Sigma(i)(j) = 0
        else
          Sigma(i)(j) = s
      }
 
    var effectives = new Array[Double](numClusters + 1) 
    var  i = 0
    for (n <- clusterSizes){
      effectives(i) = n.toDouble
      i += 1
    }
    effectives(i) = gamma
    // sample Beta
    var Beta = dirichletSample(effectives)
    var betaU : Double = Beta(numClusters)

    var initial = new ArrayBuffer[(Array[Double], Double)]
    
      for(i <- 0 to centers.length-1){

        initial += centers(i).toArray -> Beta(i)
        
      }
      
    var initialClusters = new ArrayBuffer[GlobalCluster]
    
    var globalClusters = new ArrayBuffer[((Array[Double], Double), ArrayBuffer[(Int, Int)])]
       
    // initialze mean and sigma2 for each global cluster
    var parameters = new ArrayBuffer[(Array[Double], Array[Array[Double]])]
    for(i <-0 to numClusters - 1)    
    
    // initialize workers
    workers = workers.mapPartitions{ w => {Iterator(w.next.initializeWithClustering(initial, betaU, sigma : Double, beta : Double))}}.persist
    
    time += System.currentTimeMillis() - beginTime
             
      /* println
      println("Initialisation time : "  + time)
      println */
      
      /* writerResult.write("Initialisation time : "  + time)
      writerResult.write("\n") */

    time = 0
    
    beginTime = System.currentTimeMillis()
   
    for(j <- 0 to nbDist - 1){
      
      if(j != 0){
        
        beginTime = System.currentTimeMillis()
       
        effectives = new Array[Double](numClusters + 1) 
        var i = 0
        for (c <- initialClusters){
          effectives(i) = c.size.toDouble
          i += 1
        }
        effectives(i) = gamma
        // sample Beta
        Beta = dirichletSample(effectives)
        betaU = Beta(numClusters)

        globalClusters = new ArrayBuffer[((Array[Double], Double), ArrayBuffer[(Int, Int)])]

        parameters = new ArrayBuffer[(Array[Double], Array[Array[Double]])]

        for(cluster <- initialClusters){

          globalClusters += cluster.phi -> Beta(cluster.id) -> cluster.subIds
          // post mean and sigma2 for each global cluster
          parameters += cluster.mean -> cluster.SigmaC

        }
        // initialize workers with the new global clusters
        workers = workers.mapPartitions{ w => {Iterator(w.next.startWithClustering(globalClusters, betaU))}}.persist
        // Estimate parameters
        param = workers.mapPartitions{ w => {Iterator(w.next.parameters)}}.collect
    
        distParam = DistributionsParameters(param, n)

        beta = distParam._1
        beta0 = distParam._1
        sigma = distParam._2
        sigma0 = distParam._2
        
        
     
        time += System.currentTimeMillis() - beginTime
        
        /* if(withReal == 1){
        
          val contingencyTables = workers.mapPartitions{ w => {Iterator(w.next.contingencyTable(globalClusters.length, numRealClusters))}}.collect

          val ARI = adjustedRandIndex(contingencyTables, nbWorkers, globalClusters.length, numRealClusters, n)

          println("ARI : " + ARI)

          writerResult.write("ARI : " + ARI +" ")
        
        } 
        
        val RSS = workers.mapPartitions{ w => {Iterator(w.next.localRSS)}}.collect.sum
        
        println("RSS : " + RSS)
        
        
        writerResult.write("RSS : " + RSS)
        
        writerResult.write("\n")
        
      
        if(j == nbDist - 1){
          
          val c = workers.mapPartitions{ w => {Iterator(w.next.getC)}}.collect
          for(cj <- c)
            for(ci <- cj){

              writer.write(ci._1 + " " + ci._2)
              writer.write("\n")

            }
          
          for(clust <- initialClusters){           
            for(j <- 0 to dim - 1)
              writerPhis.write(clust.phi(j)+ " ")
            writerPhis.write("\n")
          }
        
        } */
        
        beginTime = System.currentTimeMillis()
        
      }
      
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1){
        val s = (Math.pow(sigma0, 2) / (2 * beta0)) * Math.exp(- beta0 * Math.pow(tStar(i)-tStar(j), 2))
        if (s < Math.pow(10, -3) )
          Sigma0(i)(j) = 0
        else
          Sigma0(i)(j) = s
      }
  
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1){
        val s = (Math.pow(sigma, 2) / (2 * beta)) * Math.exp(- beta * Math.abs(tStar(i)-tStar(j))) 
        if (s < Math.pow(10, -3) )
          Sigma(i)(j) = 0
        else
          Sigma(i)(j) = s
      }
      
      var resultClustering = workers.mapPartitions{ w => {Iterator(w.next.gibbsSampling(sigma * mult, beta, sigma0 * mult, beta0))}}.collect
      
      val result = arrayToArrayBuffer(unionBuffers(resultClustering))
      
      var crp = new Master(sigma, beta, sigma0, beta0, n, gamma, parameters, dim, tStar)
      crp.updateExistingClasses(numClusters, result)
      
      crp.initialize      
      
      initialClusters.clear
       
      initialClusters = crp.gibbsSampling(sigma, beta, sigma0, beta0)
       
      numClusters = initialClusters.length
      
      gamma = crp.getGamma
       
      time += System.currentTimeMillis() - beginTime
             
      println
      println("distribution : " + j + " nb Clusters : " + numClusters + " time : "  + time)
      println
      
      //writerResult.write("distribution : " + j + " nb Clusters : " + numClusters + " time : "  + time + " , sigma : " + sigma + " beta : " + beta)

      time = 0
      
    }
    
    /* writer.close
    writerPhis.close
    writerResult.close */
    
  }
 
  def Cn2(n : Int) : Double = n*(n-1)/2
  
  def adjustedRandIndex(contingencyTables : Array[Array[Array[Int]]], nbPart : Int, nbGlobalClusters : Int, numRealClusters : Int, n: Int)= {
    
    var globalContingencyTable = Array.fill[Int](nbGlobalClusters, numRealClusters)(0)
        
        for(m <- 0 to nbPart - 1)
          for(i <- 0 to nbGlobalClusters - 1)
            for(j <- 0 to numRealClusters - 1)
              globalContingencyTable(i)(j) += contingencyTables(m)(i)(j)
        
        var a  = Array.fill[Int](nbGlobalClusters)(0)
        var b  = Array.fill[Int](numRealClusters)(0)
        
        for(i <- 0 to nbGlobalClusters - 1)
          for(j <- 0 to numRealClusters - 1){
            a(i) += globalContingencyTable(i)(j)
            b(j) += globalContingencyTable(i)(j)
          }
        
        var index = 0.0
        var sumCai2 = 0.0
        var sumCbj2 = 0.0
        
        for(i <- 0 to nbGlobalClusters - 1){
          sumCai2 += Cn2(a(i))
          for(j <- 0 to numRealClusters - 1){
            index += Cn2(globalContingencyTable(i)(j))
            if(i == 0)
              sumCbj2 += Cn2(b(j))
          }
          
        }
        
    val expectedIndex = sumCai2 * sumCbj2 / Cn2(n)
        
    (index - expectedIndex) / (((sumCai2 + sumCbj2) / 2) - expectedIndex)
    
  }
  
  def DistributionsParameters(distParameters : Array[Tuple3[Double, Double, Int]], n : Int)={
    
    var c = 0
    var B = 0.0
    var D = 0.0
    
    for(t <- distParameters){
  
      B += t._1
      D += t._2
      c += t._3
      
    }
    
    B /= 2 * n - 2 * c
    D /= 2 * n - 2 * c
    
    val beta = Math.log(Math.abs(D / B))
    val sigma = Math.sqrt(2 * D * Math.log(Math.abs(D / B))) 
    
    if(B <= 0)
      println("B : " + B + " beta : " + beta + " sigma : " + sigma)
    
    (beta, sigma)

    
  }
  
}
