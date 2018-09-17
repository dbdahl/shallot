\name{shallot.jar}
\alias{shallot.jar}
\title{
Initialize the \sQuote{shallot} Package
}
\description{
This function provides the path to the shallot JAR files.
}
\usage{
shallot.jar(version="", others=TRUE)
}
\arguments{
    \item{version}{If \code{""}, the package JAR files for all Scala versions are returned. Or, if \code{"2.10.*"} or \code{"2.11.*"} where \code{*} is a placeholder, the JAR file for that specific major version of Scala is returned.}
    \item{others}{Should the path to the dependencies also be returned?}
}
\examples{
shallot.jar()
shallot.jar('2.11')
shallot.jar('2.11',others=FALSE)
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}

