\name{plotting}
\alias{plotting}
\alias{pairwise.probabilities}
\alias{print.shallot.pairwiseProbability}
\alias{as.matrix.shallot.pairwiseProbability}
\alias{confidence}
\alias{plot.shallot.confidence}
\title{Plotting}
\description{
  Several functions that prepare the shallot results for plot charts.}
\usage{
pairwise.probabilities(x, parallel = TRUE)
\method{print}{shallot.pairwiseProbability}(x, ...)
\method{as.matrix}{shallot.pairwiseProbability}(x, ...)
confidence(pairwise.probabilities, partition)
\method{plot}{shallot.confidence}(x, partition = NULL, data = NULL,
show.labels = length(x$partition) <= 50, ...)
}
\arguments{
\item{x}{An object of class \code{shallot.samples.raw}.}
\item{pairwise.probabilities}{An object of class \code{shallot.pairwiseProbability}.  Defaults to \code{NULL}.}
\item{partition}{An array listing which subset each element is assigned to.}
\item{data}{The dataset itself.}
\item{parallel}{A boolean object indicating whether the partition is parallel.  Defaults to \code{TRUE}.}
\item{show.labels}{An integer indicating the number of labels for the plot.}
\item{...}{An unused argument.}
}
\details{
These functions prepare the \var{partition} for plotting purposes.  The \code{\link{pairwise.probabilities}} function calculates the distance between each observation in the dataset and each of the other observations in the dataset.

The \code{\link{print.shallot.pairwiseProbability}} function prints a message prompting the user to use the \code{\link{as.matrix.shallot.pairwiseProbability}} function to see the chart of pairwise probabilities.

The \code{\link{as.matrix.shallot.pairwiseProbability}} function converts the pairwise probabilities to matrix form and prints the matrix.

The \code{\link{confidence}} function calculates the amount of confidence by which each each element is placed in its assigned subset in the \var{partition}.

The \code{\link{plot.shallot.confidence}} function prints the confidence both with and without the original dataset itself.  Plotted without, the result is a heat chart.  Plotted with, the result is a diagram plotting the relative distances between each element in the dataset and which centroid they each gravitate towards.
}
\value{
The \code{\link{pairwise.probabilities}} function outputs an object of class \code{shallot.pairwiseProbability} containing the distances between the observations in the dataset.

The \code{\link{confidence}} function outputs an object of class \code{shallot.confidence} containing the \var{partition} itself, the confidence in array form, the confidence matrix comparing each element to each of the other elements, the exemplar indicating what the largest number of elements across the partition that gravitate towards each subset (although the numbers are high enough that there must be a very large amount of overlap), and the order of appearance for the dataset elements.

The \code{\link{plot.shallot.confidence}} function displays a chart that illustrates the results.  Without the original dataset itself included, the chart takes the form of a heat chart.  With the dataset included, the chart is a diagram that shows the relative distances between each element and each other and also indicates where the centroids are.  This chart examines these distances with only two attributes emphasized to compare how strongly each attribute affects the data.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{sample.partition}}
\code{\link{estimate.partition}}
}
\examples{
\dontrun{
## Demonstrate pairwise.probabilities, as.matrix, confidence, and plot.
data <- Theoph[, c("Wt", "Dose", "Time", "conc", "Subject")]
truth <- as.numeric(levels(data[, ncol(data)]))[data[, ncol(data)]]
distance <- dist(scale(data[, -ncol(data)]))
decay <- decay.exponential(9.0, fixed = TRUE)
attract <- attraction(distance(as.matrix(distance)),
    permutation(n.items = nrow(data), fixed = FALSE), decay)
mass <- mass(1.0, fixed = TRUE)
discount <- discount(0.1, fixed = TRUE)
distribution <- ewens.pitman.attraction(mass, discount, attract)
samples <- collect(distribution, 500)
est <- estimate.partition(samples, max.subsets = 12, max.scans = 10)
pp <- pairwise.probabilities(samples)
print(pp)
print(as.matrix(pp))
conf <- confidence(pp, est)
print(conf)
plot(conf)
plot(conf, data = data)
}
}

