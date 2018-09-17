library("shallot")

shallot.load(java.heap.maximum="2G")

mass(3.2)
mass(3.2,fixed=FALSE)
mass(2.9,0.3,fixed=FALSE)

discount(0.2)
discount(2.9,0.3,fixed=FALSE)

covariates <- scale(iris[,-5])
d <- distance(as.matrix(dist(covariates)))

permutation(4,3,2,1)
permutation(sample(150))
permutation(n.items=10,fixed=FALSE)
permutation(1,2,3,fixed=FALSE)
permutation(1,2,3,fixed=TRUE)

decay.exponential(3.2)
decay.reciprocal(1,4,fixed=FALSE)
decay.subtraction(fixed=FALSE)
decay <- decay.exponential(fixed=TRUE)

attraction(d,permutation(n.items=d$n.items,fixed=FALSE),decay)
attraction(d,permutation(1:d$n.items,fixed=TRUE),decay)
attraction(d,permutation(n.items=d$n.items,fixed=FALSE),decay)
a <- attraction(d,permutation(sample(d$n.items)),decay)
a
apply(as.matrix(a),1,sum)

pd <- ewens(mass(),30)
pd

pd <- ewens.pitman(mass(),discount(),30)
pd

pd <- ewens.attraction(mass(),a)
pd

pd <- ewens.pitman.attraction(mass(),discount(),a)
pd

pd <- ddcrp(mass(),a)
pd

ppmx.similarity.continuous(as.matrix(iris[,-5]),mu0=0,n0=1,alpha=3,beta=2)
ppmx.similarity.categorical(as.matrix(iris[,5]),alpha=0.4)
ppmx.similarity.categorical(as.matrix(iris[,5]),levels=c("setosa","versicolor","virginica"),alpha=0.5)
ppmx.similarity.categorical(as.matrix(iris[,5]),levels=c("setosa","versicolor","virginica","others"))

s1 <- ppmx.similarity.continuous(as.matrix(iris[,-5])); s1
s2 <- ppmx.similarity.categorical(as.matrix(iris[,5])); s2  # In practice, you wouldn't use the clustering truth!
s <- ppmx.similarity.composite(s1,s2); s

pd <- ppmx(mass(),s)
pd

samples <- collect(pd,n.draws=100)
samples <- collect(samples,n.draws=100)

y <- process(samples)
str(y)
plot(y$entropy,type="l")

mean(y$n.subsets)
mean(apply(y$partitions,1,max))

pp <- pairwise.probabilities(samples)
pp
as.matrix(pp)[1:5,1:5]

e <- estimate.partition(samples,pp,max.subsets=3)
e

conf <- confidence(pp,e)
plot(conf)



