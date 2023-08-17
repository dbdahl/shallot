# aibd: Attraction Indian Buffet Distribution

We propose the attraction Indian buffet distribution (AIBD), a distribution for binary feature matrices influenced by pairwise similarity information. Binary feature matrices are used in Bayesian models to uncover latent variables (i.e., features) that explain observed data. The Indian buffet process (IBP) is a popular exchangeable prior distribution for latent feature matrices. In the presence of additional information, however, the exchangeability assumption is not reasonable or desirable. The AIBD can incorporate pairwise similarity information, yet it preserves many properties of the IBP, including the distribution of the total number of features. Thus, much of the interpretation and intuition that one has for the IBP directly carries over to the AIBD. A temperature parameter controls the degree to which the similarity information affects feature-sharing between observations. Unlike other nonexchangeable distributions for feature allocations, the probability mass function of the AIBD has a tractable normalizing constant, making posterior inference on hyperparameters straight-forward using standard MCMC methods. A novel posterior sampling algorithm is proposed for the IBP and the AIBD. We demonstrate the feasibility of the AIBD as a prior distribution in feature allocation models and compare the performance of competing methods in simulations and an application.

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


