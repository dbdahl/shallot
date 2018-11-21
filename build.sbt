name := "shallot"

organization := "org.ddahl"

//version := "0.4.6"
version := "0.4.6-SNAPSHOT"

scalaVersion := "2.12.7"

crossScalaVersions := Seq("2.11.12", "2.12.7")

scalacOptions ++= Seq( "-deprecation", "-unchecked", "-feature" )

libraryDependencies ++= Seq(
  "org.ddahl" %% "sdols" % "1.8-SNAPSHOT",
  "org.apache.commons" % "commons-math3" % "3.6.1" withSources(),
  "org.scalatest" %% "scalatest" % "3.0.5" % "test"
)

retrieveManaged := true

mainClass in (Compile,run) := Some("org.ddahl.shallot.example.Main")

