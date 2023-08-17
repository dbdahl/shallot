# shallot: Random Partition Distribution Indexed by Pairwise Information

We propose a random partition distribution indexed by pairwise similarity information such that partitions compatible with the similarities are given more probability. The use of pairwise similarities, in the form of distances, is common in some clustering algorithms (e.g., hierarchical clustering), but we show how to use this type of information to define a prior partition distribution for flexible Bayesian modeling. A defining feature of the distribution is that it allocates probability among partitions within a given number of subsets, but it does not shift probability among sets of partitions with different numbers of subsets. Our distribution places more probability on partitions that group similar items yet keeps the total probability of partitions with a given number of subsets constant. The distribution of the number of subsets (and its moments) is available in closed-form and is not a function of the similarities. Our formulation has an explicit probability mass function (with a tractable normalizing constant) so the full suite of MCMC methods may be used for posterior inference. We compare our distribution with several existing partition distributions, showing that our formulation has attractive properties. We provide three demonstrations to highlight the features and relative performance of our distribution.

## Installation

In R, install the package by executing:

```R
install.packages("commonsMath")
install.packages("remotes")
remotes::install_github("dbdahl/rscala/R/rscala")
remotes::install_github("dbdahl/shallot/R/shallot")
```

## Resources

* [Paper](https://doi.org/10.1080/01621459.2016.1165103) describes the methods implemented in this package.

```R
library(shallot)
example(shallot)
help(shallot)
```


