\name{estimate.partition}
\alias{estimate.partition}
\title{Estimate Partition}
\description{
This function returns a partition that summarizes the partition distribution using the least-square clustering method (Dahl 2006), with extensions to
perform greedy optimization and limit the number of subsets.}
\usage{
estimate.partition(x, pairwise.probabilities = NULL,
    max.subsets = 0, max.scans = 0, parallel = TRUE)
}
\arguments{
\item{x}{An object from the \code{\link{sample.partitions}} function.}
\item{pairwise.probabilities}{An object of class \code{shallot.pairwiseProbability} obtained from \code{\link{pairwise.probabilities}}.  If not supplied, it will be computed from \var{x}.}
\item{max.subsets}{An integer limiting the number of subsets.  Defaults to \code{0}, which does not impose a constraint on the number of subsets.}
\item{max.scans}{An integer controlling the greedy search.  Defaults to \code{0}, which disables the greedy search.}
\item{parallel}{Should all of the CPU cores should be used?  Defaults to \code{TRUE}.}
}
\value{
A partition as a vector of cluster labels.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{process.partitions}},
\code{\link{sample.partitions}},
\code{\link{plot.partition}},
\code{\link{adj.rand.index}}
}
\references{
\href{https://doi.org/10.1080/01621459.2016.1165103}{Dahl, D. B., Day, R., and Tsai, J. (2017), Random Partition Distribution Indexed by Pairwise Information, \emph{Journal of the American Statistical Association}, accepted. <DOI:10.1080/01621459.2016.1165103>}

Dahl, D. B. (2006), Model-Based Clustering for Expression Data via a Dirichlet Process Mixture Model, in \emph{Bayesian Inference for Gene Expression and Proteomics}, Kim-Anh Do, Peter Mueller, Marina Vannucci (Eds.), Cambridge University Press.
}
\examples{
\dontrun{
example(shallot)
}
}
