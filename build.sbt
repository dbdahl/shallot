name := "shallot"

organization := "org.ddahl"

//version := "0.4.5"
version := "0.4.4-SNAPSHOT"

scalaVersion := "2.12.6"

crossScalaVersions := Seq("2.11.12", "2.12.6")

scalacOptions ++= Seq( "-deprecation", "-unchecked", "-feature" )

libraryDependencies ++= Seq(
  "org.ddahl" %% "sdols" % "1.6-SNAPSHOT",
  "org.apache.commons" % "commons-math3" % "3.6.1" withSources(),
  "org.scalatest" %% "scalatest" % "3.0.1" % "test"
)

retrieveManaged := true

mainClass in (Compile,run) := Some("org.ddahl.shallot.example.Main")

