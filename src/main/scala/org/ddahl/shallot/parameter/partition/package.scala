package org.ddahl.shallot
package parameter

package object partition {

  type Subset[A] = org.ddahl.sdols.clustering.Cluster[A]
  val Subset = org.ddahl.sdols.clustering.Cluster

  type Partition[A] = org.ddahl.sdols.clustering.Clustering[A]
  val Partition = org.ddahl.sdols.clustering.Clustering

}

