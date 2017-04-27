\name{adj.rand.index}
\alias{adj.rand.index}
\title{Adjusted Rand Index}
\description{
    Function calculating adjusted Rand index of two possible distributions of a given dataset.}
\usage{
adj.rand.index(c1, c2)
}
\arguments{
\item{c1}{A named array containing the estimated distribution in integer form of the observations of the dataset when cutree is executed at a given number of subsets or when cutreeDynamic is executed.}
\item{c2}{An array containing the distribution in integer form of the observations of the dataset based on the attribute that is selected to be the so-called truth of how the dataset should be distributed.}
}
\details{
The \code{\link{adj.rand.index}} function takes as its input two possible distributions of a given dataset and compares them to see how similar they are.  That is, the function compares how many times each observation is classified into the same-numbered subset.  For example, suppose \var{c1} indicates that the dataset should be distributed as \code{1 1 2 1 1 2 3}, while \var{c2} indicates that the distribution should be \code{1 1 2 1 2 3 1}.  Only the first four observations would match up.  As a result, the Rand index for this example would be approximately \code{57} percent.  However, this function factors out any expected matchups between the two distributions and only computes the adjusted Rand index.  As a result, the adjusted Rand index for this example would actually be \code{0.1428571}.

Although any two possible distributions can be used as the input, the second distribution is invariably the dataset attribute that has been selected to represent the so-called truth of how the dataset's observations are distributed.  All calculations determining the partition of the dataset only use the other attributes and leave out the truth to provide a pure comparison.  The first input distribution represents the results of these calculations courtesy of the \code{\link{estimate.partition}} function, which builds on a chain of other functions in this package that are designed to create this possible distribution while allowing the user flexibility in determining the values of several variables that can significantly alter the distribution.  By altering the distribution, the user can therefore alter the adjusted Rand index.
}
\value{
A number, representing the adjusted Rand index between the two distributions.  The bigger the number, the more similar the distributions truly are.
}
\note{
This function is not an improvement over the arandi function in the mcclust package.  Arandi is designed to be flexible between returning the adjusted Rand index and the unadjusted Rand index, whereas \code{\link{adj.rand.index}} only returns the adjusted Rand index.  The internal coding is similar to that of arandi and does not offer any significant improvement.
}
\author{David B. Dahl \email{dahl@stat.byu.edu}}
\seealso{
\code{\link{estimate.partition}}
}
\examples{
## Load example dataset.
data <- Theoph[, c("Wt", "Dose", "Time", "conc", "Subject")]
truth <- as.numeric(levels(data[, ncol(data)]))[data[, ncol(data)]]
distance <- dist(scale(data[, -ncol(data)]))

## Demonstrate adj.rand.index on average hclust method.
adj.rand.index(cutree(hclust(distance, method = "average"), k = length(unique(truth))), truth)

## Demonstrate adj.rand.index on centroid hclust method.
adj.rand.index(cutree(hclust(distance, method = "centroid"), k = length(unique(truth))), truth)

## Demonstrate adj.rand.index on complete hclust method.
adj.rand.index(cutree(hclust(distance, method = "complete"), k = length(unique(truth))), truth)

## Demonstrate adj.rand.index on median hclust method.
adj.rand.index(cutree(hclust(distance, method = "median"), k = length(unique(truth))), truth)

## Demonstrate adj.rand.index on single hclust method.
adj.rand.index(cutree(hclust(distance, method = "single"), k = length(unique(truth))), truth)

## Demonstrate adj.rand.index on ward.D hclust method.
adj.rand.index(cutree(hclust(distance, method = "ward.D"), k = length(unique(truth))), truth)

## Demonstrate adj.rand.index on ward.D2 hclust method.
adj.rand.index(cutree(hclust(distance, method = "ward.D2"), k = length(unique(truth))), truth)
}

