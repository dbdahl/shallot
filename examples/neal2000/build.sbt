name := "neal2000"

version := "0.0.2"

scalaVersion := "2.12.0"

scalacOptions ++= Seq( "-deprecation", "-unchecked", "-feature" )

libraryDependencies ++= Seq(
  "org.apache.commons" % "commons-math3" % "3.6.1"
)

mainClass in (Compile,run) := Some("neal2000.Main")

