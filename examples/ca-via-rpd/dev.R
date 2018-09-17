options("rscala.command.line.options"="-J-Xmx16g")
library(shallot)

x <- sample.partition(ewens(mass(1,fixed=TRUE),10),100)
x <- sample.partition(ewens.pitman(mass(1,fixed=FALSE),discount(0.0,fixed=TRUE),10),100,parallel=TRUE)

d <- dist(scale(USArrests))
p <- permutation(1:attr(d,"Size"),fixed=TRUE)
decay <- decay.exponential(temperature=temperature(),distance=d)

attr <- attraction(p,decay)

x <- sample.partition(ewens.attraction(mass(1,fixed=FALSE),attr),100,parallel=TRUE)
x <- sample.partition(ewens.pitman.attraction(mass(1,fixed=FALSE),discount(fixed=FALSE),attr),1000,parallel=TRUE)

y <- process.partition(x)
y <- process.partition(x,as.matrix=TRUE,sample.parameter=function() rnorm(1))

w <- pairwise.probabilities(x)
as.matrix(w)

e <- estimate.partition(x,w,max.subsets=3)
table(e)

cf <- confidence(w,e)

plot(cf)
plot(cf,data=USArrests)

