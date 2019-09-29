package org.ddahl.sdols
package clustering

import org.apache.commons.math3.util.FastMath.log
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import scala.annotation.tailrec
import scala.concurrent.{Await, Future}
import scala.concurrent.duration._
import scala.concurrent.ExecutionContext.Implicits.global
import scala.collection.parallel.immutable.ParVector

object ClusteringSummary {

  implicit val ordering = CrossCompatibility.doubleOrdering

  // Assumes that 'clusterings' is a non-empty sequence of clusterings.
  def expectedPairwiseAllocationMatrix[A](clusterings: Seq[Clustering[A]]): Array[Array[Double]] = {
    val nItems = clusterings(0).nItems
    val x = Array.ofDim[Int](nItems,nItems)
    clusterings.foreach { clustering =>
      clustering.foreach { cluster =>
        val y = cluster.toArray
        var i = 0
        while ( i < y.length ) {
          val ii = y(i)
          var j = i+1
          while ( j < y.length ) {
            val jj = y(j)
            x(ii)(jj) += 1
            x(jj)(ii) += 1
            j += 1
          }
          x(ii)(ii) += 1
          i += 1
        }
      }
    }
    val cl: Double = clusterings.length
    x.map(_.map(_/cl))
  }

  // Assumes that 'clusterings' is a non-zero-length array of arrays of equal lengths.
  def expectedPairwiseAllocationMatrix(clusterings: Array[Array[Int]]): Array[Array[Double]] = {
    val n = clusterings(0).length
    val x = Array.ofDim[Int](n,n)
    var k = 0
    while ( k < clusterings.length ) {
      val p = clusterings(k)
      var i = 0
      while ( i < n ) {
        val xi = x(i)
        val pi = p(i)
        var j = 0
        while ( j < n ) {
          if ( pi == p(j) ) xi(j) += 1
          j += 1
        }
        i += 1
      }
      k += 1
    }
    val cl: Double = clusterings.length
    x.map(_.map(_/cl))
  }

  def confidenceComputations(clustering: Array[Int], pam: Array[Array[Double]]): (Array[Double], Array[Array[Double]], Array[Int], Array[Int], Array[Int]) = {
    val clustering2 = Clustering(i => i, clustering)
    val confidence = (0 until clustering2.nItems).map { i =>
      val pami = pam(i)
      val cluster = clustering2.clusterFor(i)
      cluster.foldLeft(0.0) { (sum, j) => sum + pami(j) } / cluster.size
    }.toArray
    val sizes = clustering2.map( cluster => (cluster,cluster.size) ).toMap
    val clusters = clustering2.toList.sortWith( (c1,c2) => sizes(c1) > sizes(c2) )
    val confidenceMatrixLabels = clusters.map(_.parameter).toArray
    val confidenceMatrix = clusters.map { cluster1 =>
      clusters.map { cluster2 =>
        cluster1.foldLeft(0.0) { (sum1, i) =>
          val pami = pam(i)
          sum1 + cluster2.foldLeft(0.0) { (sum2, j) =>
            sum2 + pami(j)
          }
        } / ( cluster1.size * cluster2.size )
      }.toArray
    }.toArray
    val exemplar = clusters.map { cluster =>
      cluster.toList.sortWith { (i, j) =>
        if (confidence(i) > confidence(j)) true
        else if (confidence(i) == confidence(j)) i < j
        else false
      }
    }.map(_.head).toArray
    val orderList = clusters.zip(exemplar).flatMap { x =>
      @tailrec
      def findBest(candidates: List[Int], members: List[Int]): List[Int] = {
        if ( ! candidates.isEmpty ) {
          val best = candidates.map { y =>
            (members.foldLeft(0.0) { (sum, j) => sum + pam(y)(j) }, y)
          }.sortWith(_._1 > _._1).head._2
          findBest(candidates.filterNot(_ == best), best :: members)
        } else members
      }
      findBest((x._1.x - x._2).toList, List[Int](x._2))
    }.toArray
    (confidence, confidenceMatrix, confidenceMatrixLabels, orderList, exemplar)
  }

  def lowerBoundVariationOfInformation[A](clustering: Clustering[A], pam: Array[Array[Double]]): Double = {
    val nItems = clustering.nItems
    var sum1 = 0.0
    var i = 0
    while ( i < nItems ) {
      var sum2 = 0
      var sum3 = 0.0
      var sum4 = 0.0
      var j = 0
      while ( j < nItems ) {
        sum2 += (if (clustering.paired(i, j)) 1 else 0)
        sum3 += pam(i)(j)
        sum4 += (if (clustering.paired(i, j)) pam(i)(j) else 0.0)
        j += 1
      }
      sum1 += log(2, sum2) + log(2, sum3) - 2*log(2, sum4)
      i += 1
    }
    sum1/nItems
  }

  private def lowerBoundVariationOfInformationEngine[A](clustering: Clustering[A], pam: Array[Array[Double]]): Double = {
    var sum = 0.0
    val clusters = clustering.toArray
    var k = 0
    while ( k < clusters.length ) {
      val cluster = clusters(k).toArray
      sum += cluster.length * log(2, cluster.length)
      var sum2 = 0.0
      var ii = 0
      while ( ii < cluster.length ) {
        var sum3 = 0.0
        val pamii = pam(cluster(ii))
        var jj = 0
        while ( jj < cluster.length ) {
          sum3 += pamii(cluster(jj))
          jj += 1
        }
        sum2 += log(2,sum3)
        ii += 1
      }
      sum -= 2*sum2
      k += 1
    }
    sum
  }

  def sumOfAbsolutesSlow[A](clustering: Clustering[A], pam: Array[Array[Double]]): Double = {
    new Array2DRowRealMatrix(clustering.pairwiseAllocationMatrix.map(_.map(_.toDouble)),false).subtract(new Array2DRowRealMatrix(pam,false)).getDataRef.foldLeft(0.0)(_+_.foldLeft(0.0)((a,b) => a + b.abs))
  }

  def sumOfSquaresSlow[A](clustering: Clustering[A], pam: Array[Array[Double]]): Double = {
    new Array2DRowRealMatrix(clustering.pairwiseAllocationMatrix.map(_.map(_.toDouble)),false).subtract(new Array2DRowRealMatrix(pam,false)).getDataRef.foldLeft(0.0)(_+_.foldLeft(0.0)((a,b) => a + b*b))
  }

  private def binderOffset[A](pam: Array[Array[Double]], f: Double => Double) = {
    var offset = 0.0
    var i = 1
    while ( i < pam.length ) {
      val pi = pam(i)
      var j = 0
      while ( j < i ) {
        offset += f(pi(j))
        j += 1
      }
      i += 1
    }
    2*offset
  }

  def sumOfAbsolutes[A](clustering: Clustering[A], pam: Array[Array[Double]]): Double = {
    val pamTransform = pam.map(_.map(x => 2-4*x))
    binderOffset(pam,(x: Double) => x.abs) + binderEngine(clustering,pamTransform)
  }

  def sumOfSquares[A](clustering: Clustering[A], pam: Array[Array[Double]]): Double = {
    val pamTransform = pam.map(_.map(x => 2-4*x))
    binderOffset(pam,(x: Double) => x*x) + binderEngine(clustering,pamTransform)
  }

  private def binderEngine[A](clustering: Clustering[A], pamTransform: Array[Array[Double]]): Double = {
    var sum = 0.0
    clustering.foreach { cluster =>
      val y = cluster.toArray
      var i = 0
      while ( i < y.length ) {
        val xx = pamTransform(y(i))
        var j = i+1
        while ( j < y.length ) {
          sum += xx(y(j))
          j += 1
        }
        i += 1
      }
    }
    sum
  }

  private val emptyTuple = (Cluster.empty[Null](null),0.0)

  private def binderShortcutEngine(i: Int, clusteringWithoutI: Clustering[Null], maxSize: Int, pamTransform: Array[Array[Double]]): List[(Cluster[Null], Double)] = {
    val pamTransformi = pamTransform(i)
    val candidates = clusteringWithoutI.map(cluster => (cluster, cluster.foldLeft(0.0) { (sum, j) => sum + pamTransformi(j) })).toList
    if ((maxSize <= 0) || (clusteringWithoutI.size < maxSize)) emptyTuple :: candidates
    else candidates
  }

  private def lowerBoundVariationOfInformationShortcutEngine(i: Int, clusteringWithoutI: Clustering[Null], maxSize: Int, pam: Array[Array[Double]]): List[(Cluster[Null], Double)] = {
    val candidates = clusteringWithoutI.map(cluster => {
      (cluster,lowerBoundVariationOfInformationEngine(clusteringWithoutI.add(i,cluster),pam))
    }).toList
    if ((maxSize <= 0) || (clusteringWithoutI.size < maxSize)) (Cluster.empty[Null](null),lowerBoundVariationOfInformationEngine(clusteringWithoutI.add(i,Cluster.empty[Null](null)),pam)) :: candidates
    else candidates
  }

/*
  def minAmongDraws(candidates: Array[Array[Int]], maxSize: Int, multicore: Boolean, loss: String, pamOption: Option[Array[Array[Double]]] = None): Clustering[Null] = {
    minAmongDraws(candidates.map(Clustering.apply),maxSize,multicore,loss,pamOption)
  }
*/

  def minAmongDraws[A](candidates: Seq[Clustering[A]], maxSize: Int, multicore: Boolean, loss: String, pamOption: Option[Array[Array[Double]]]): Clustering[A] = {
    if ( candidates.isEmpty ) throw new IllegalArgumentException("'candidates' cannot be empty.")
    val pam = pamOption.getOrElse(expectedPairwiseAllocationMatrix(candidates))
    val (lossEngine, shortcutEngine, pamTransform) = getLoss[A](maxSize, loss, pam)
    val iter = (if ( multicore ) ParVector(candidates:_*) else candidates).iterator
    iter.minBy { clustering =>
      if ( ( maxSize > 0 ) && ( clustering.size > maxSize ) ) Double.PositiveInfinity
      else lossEngine(clustering, pamTransform)
    }
  }

  private def sequentiallyAllocatedLatentStructureOptimizationOld(initial: Clustering[Null], maxSize: Int, permutation: List[Int], pamTransform: Array[Array[Double]], lossEngine: (Clustering[Null],Array[Array[Double]]) => Double): Clustering[Null] = {
    var clustering = initial
    for ( i <- permutation ) {
      val candidates = clustering.map(cluster => clustering.add(i,cluster)).toList
      val candidates2 = if ( ( maxSize <= 0 ) || ( clustering.size < maxSize ) ) {
        clustering.add(Cluster(null,i)) :: candidates
      } else candidates
      clustering = candidates2.minBy(lossEngine(_,pamTransform))
    }
    clustering
  }

  def mkShortcutEngineBinder(maxSize: Int, pamTransform: Array[Array[Double]]) = (i: Int, clusteringWithoutI: Clustering[Null]) => {
    binderShortcutEngine(i, clusteringWithoutI, maxSize, pamTransform)
  }

  def mkShortcutEngineLowerBoundVariationOfInformation(maxSize: Int, pam: Array[Array[Double]]) = (i: Int, clusteringWithoutI: Clustering[Null]) => {
    lowerBoundVariationOfInformationShortcutEngine(i, clusteringWithoutI, maxSize, pam)
  }

  private def sequentiallyAllocatedLatentStructureOptimization(initial: Clustering[Null], maxScans: Int, permutation: List[Int], engine: (Int,Clustering[Null]) => List[(Cluster[Null], Double)]): (Clustering[Null], Int) = {
    var clustering = initial
    var firstPass = true
    var notDone = true
    var scanCounter = -1
    while ( firstPass || notDone ) {
      scanCounter += 1
      val previousClustering = clustering
      for (i <- permutation) {
        if ( ! firstPass ) clustering = clustering.remove(i)
        val candidates = engine(i,clustering)
        clustering = clustering.add(i, candidates.minBy(_._2)._1)
      }
      if ( firstPass ) firstPass = false
      notDone = ( clustering != previousClustering ) && ( scanCounter < maxScans )
    }
    (clustering, scanCounter)
  }

  private def getLoss[A](maxSize: Int, loss: String, pam: Array[Array[Double]]) = {
    loss match {
      case "binder" | "squaredError" | "absoluteError" => {
        val pamTransform = pam.map(_.map(x => 0.5-x))
        (binderEngine[A] _, mkShortcutEngineBinder(maxSize, pamTransform), pamTransform)
      }
      case "lowerBoundVariationOfInformation" => {
        (lowerBoundVariationOfInformationEngine[A] _, mkShortcutEngineLowerBoundVariationOfInformation(maxSize, pam), pam)
      }
    }
  }

  // Oops, when useOldAlgorithm, we don't do any sweetening scans.
  def sequentiallyAllocatedLatentStructureOptimization(nCandidates: Int, budgetInSeconds: Double, pam: Array[Array[Double]], maxSize: Int, maxScans: Int, multicore: Boolean, loss: String, useOldImplementation: Boolean): (Clustering[Null], Int, Int) = {
    val (lossEngine, shortcutEngine, pamTransform) = getLoss[Null](maxSize, loss, pam)
    val rng = new scala.util.Random()   // Thread safe!
    val nItems = pam.length
    val ints = List.tabulate(nItems) { identity }
    val empty = Clustering.empty[Null]()
    val nCores = if ( multicore ) Runtime.getRuntime.availableProcessors else 1
    val nCandidatesPerThread = (nCandidates - 1) / nCores + 1
    val budgetInMillis = if ( budgetInSeconds <= 0 ) Long.MaxValue else budgetInSeconds * 1000L
    val start = System.currentTimeMillis
    val futures = List.fill(nCores) { Future {
      var counter = 0
      var minScore = Double.PositiveInfinity
      var best: (Clustering[Null], Int) = (null,-1)
      while ( (counter < nCandidatesPerThread) && ( System.currentTimeMillis - start <= budgetInMillis ) ) {
        counter += 1
        val permutation = rng.shuffle(ints)
        val candidate = if ( ! useOldImplementation ) sequentiallyAllocatedLatentStructureOptimization(empty,maxScans,permutation,shortcutEngine)
        else (sequentiallyAllocatedLatentStructureOptimizationOld(empty,maxSize,permutation,pamTransform,lossEngine),-1)
        val score = lossEngine(candidate._1,pamTransform)
        if ( score < minScore ) {
          minScore = score
          best = candidate
        }
      }
      (best._1, best._2, minScore, counter)
    }}
    val seq = Await.result(Future.sequence(futures), Duration.Inf)
    val nCandidatesInPractice = seq.foldLeft(0)( (sum,tuple) => sum + tuple._4 )
    val best = seq.minBy(_._3)
    (best._1, best._2, nCandidatesInPractice)
  }

}

