name := "shallot"

organization := "org.ddahl"

version := "0.4.6.1"
//version := "0.4.6.1-SNAPSHOT"

scalaVersion := "2.12.8"

crossScalaVersions := Seq("2.11.12", "2.12.8")

scalacOptions ++= Seq( "-deprecation", "-unchecked", "-feature" )

mainClass in (Compile,run) := Some("org.ddahl.shallot.example.Main")

libraryDependencies ++= Seq(
  "org.ddahl" %% "sdols" % "1.7.3.1",
  "org.ddahl" %% "commonsmath" % "1.2.2.4",
  "org.apache.commons" % "commons-math3" % "3.6.1" withSources(),
  "org.scalatest" %% "scalatest" % "3.0.5" % "test"
)

resolvers += Resolver.bintrayRepo("dahl", "maven")

licenses := List(("Apache-2.0",url("https://www.apache.org/licenses/LICENSE-2.0")))

publishMavenStyle := true

