# HD4C

This is a Distributed Clustering with Spark based on Dirichlet Process Mixture adapted for high dimensional data, this approach is described in the following paper:

Khadidja Meguelati, Bénédicte Fontez, Nadine Hilgert, Florent Masseglia. [High Dimensional Data Clustering by means of Distributed Dirichlet Process Mixture Models](https://hal-lirmm.ccsd.cnrs.fr/lirmm-02364411). IEEE Big Data : IEEE International Conference on Big Data, Dec 2019, Los Angeles, CA, USA.

Please kindly cite our paper if the code helps you. Thank you.

```
@inproceedings{meguelati:lirmm-02364411,
  TITLE = {{High Dimensional Data Clustering by means of Distributed Dirichlet Process Mixture Models}},
  AUTHOR = {Meguelati, Khadidja and Fontez, B{\'e}n{\'e}dicte and Hilgert, Nadine and Masseglia, Florent},
  URL = {https://hal-lirmm.ccsd.cnrs.fr/lirmm-02364411},
  BOOKTITLE = {{IEEE International Conference on Big Data (IEEE BigData)}},
  ADDRESS = {Los-Angeles, United States},
  YEAR = {2019},
  MONTH = Dec,
  KEYWORDS = {Gaussian random process ; Dirichlet Process Mixture Model ; Clustering ; Parallelism ; Reproducing Kernel Hilbert Space},
  PDF = {https://hal-lirmm.ccsd.cnrs.fr/lirmm-02364411/file/IEEE_BigData_2019__HAL_.pdf},
  HAL_ID = {lirmm-02364411},
  HAL_VERSION = {v1},
}
```

## Requirements
HD4C works with [Apache Spark](http://spark.apache.org). In order to run it you must download and install [Spark Release 2.0.0](https://spark.apache.org/releases/spark-release-2-0-0.html).
The code is written in [Scala](https://www.scala-lang.org/), install [Scala 2.11.6](https://www.scala-lang.org/download/2.11.6.html)

## Building
We use maven to build it, Use the given [pom.xml](https://github.com/khadidjaM/HD4C/blob/master/pom.xml) file to build an executable jar containing all the dependencies.

## Use
To execute HD4C use the following command :
```
$SPARK_HOME/bin/spark-submit HD4C-jar-with-dependencies.jar <dimensions> <number of workers> <number of distributions> <target to data file> <number of clusters for Kmeans> <number of real clusters> <1 if the real number of clusters is known>
```
### Necessary parameters 
1. **dimensions:** the number of dimensions
2. **number of workers:**
3. **number of distributions:** in each distribution we perform several iterations of Gibbs Sampling on each worker and a synchronsation at the master level  
4. **target to data file:** The data file should be as follow :
  * each data in a line
  * values are seperated by space " "
  * if the ground truth is known, the data file should contain the label of the real cluster for each data in the last column.
5. **number of clusters for Kmeans:** the initialization Of DPM is done by a K-means step, you should indicate the number of clusters for Kmeans initialization
6. **number of real clusters:** if the ground truth is known, indicate the number of real clusters, else you can enter 0 
7. **real clusters are known:** if the ground truth is known, enter 1 else you can enter 0
