package org.ddahl.sdols
package clustering

final class Cluster[A] private[clustering](override val size: Int, val x: Set[Int], val parameter: A) extends Iterable[Int] {

  private val checks = false

  if (checks) {
    if (size != x.size) throw new RuntimeException("Internal error")
  }

  protected[clustering] def add(i: Int) = {
    if (checks) {
      if (contains(i)) throw new IllegalArgumentException("Cluster already contains " + i + ".")
    }
    new Cluster(size + 1, x + i, parameter)
  }

  protected[clustering] def add(cluster: Int*) = {
    if (checks) {
      if (cluster.exists(x.contains)) throw new IllegalArgumentException("Cluster already contains at least one of these items.")
    }
    new Cluster(size + cluster.size, x ++ cluster, parameter)
  }

  protected[clustering] def remove(i: Int) = {
    if (checks) {
      if (!contains(i)) throw new IllegalArgumentException("Cluster does not contain " + i + ".")
    }
    new Cluster(size - 1, x - i, parameter)
  }

  protected[clustering] def remove(cluster: Int*) = {
    if (checks) {
      if (!cluster.forall(x.contains)) throw new IllegalArgumentException("Cluster does not contain all of these items.")
    }
    new Cluster(size - cluster.size, x -- cluster, parameter)
  }

  protected[clustering] def replace(parameter: A) = {
    new Cluster(size, x, parameter)
  }

  def contains(i: Int) = x.contains(i)

  def contains(cluster: Int*) = cluster.forall(x.contains)

  override def equals(other: Any) = other match {
    case that: Cluster[A] =>
      if (that.size != size) false
      else that.x == x
    case _ => false
  }

  override def hashCode: Int = x.hashCode

  def iterator = x.iterator

  override def toString = {
    "({" + toList.sortWith(_ < _).mkString(",") + "}," + parameter + ")"
  }

  def toStringTerse = {
    "{" + toList.sortWith(_ < _).mkString(",") + "}"
  }

  def write(objOutputStream: java.io.ObjectOutputStream) = {
    objOutputStream.writeInt(size)
    iterator.foreach(objOutputStream.writeInt)
    objOutputStream.writeObject(parameter)
  }

}

object Cluster {

  def empty[A](parameter: A): Cluster[A] = new Cluster(0, Set[Int](), parameter)

  def apply[A](parameter: A, i: Int): Cluster[A] = new Cluster(1, Set[Int](i), parameter)

  def apply[A](parameter: A, i: Iterable[Int]): Cluster[A] = {
    val set = i.toSet
    new Cluster(set.size, set, parameter)
  }

  def apply[A](parameter: A, i: Set[Int]): Cluster[A] = new Cluster(i.size, i, parameter)

  def apply[A](parameter: A, i: Int*): Cluster[A] = apply(parameter, i)

  def apply[A](objInputStream: java.io.ObjectInputStream): Cluster[A] = {
    val nItems = objInputStream.readInt()
    val seq = Seq.fill(nItems) { objInputStream.readInt() }
    val parameter = objInputStream.readObject().asInstanceOf[A]
    val set = seq.toSet
    new Cluster(set.size,set,parameter)
  }
  
}

