d <- dist(scale(USArrests))
as.matrix(d)

options("rscala.command.line.options"="-J-Xmx16g")
library(shallot)

mass()
mass(2.3)
mass(1,fixed=FALSE)
mass(3,4,fixed=FALSE)

discount()
discount(0.1)
discount(0.1,fixed=FALSE)
discount(3,4,fixed=FALSE)

adj.rand.index(c(1,1,1,2,2,2,3,3,3),c(1,1,1,1,1,1,2,2,2))

permutation(n.items=30,fixed=FALSE)
permutation(2,3,4,1)
permutation(2,3,4,1,fixed=FALSE)

temperature()
t1 <- temperature(2,4,fixed=FALSE)
t2 <- temperature(2)
t3 <- temperature(Inf)

decay.reciprocal(t1,d)
decay.reciprocal(t2,d)
decay.reciprocal(t3,d)
decay.exponential(t1,d)
decay.exponential(t2,d)
decay.exponential(t3,d)
decay.subtraction(t1,d)
decay.subtraction(t2,d)
decay.subtraction(t3,d)

a <- attraction(10)
a <- attraction(permutation(n.items=attr(d,"Size"),fixed=FALSE),decay.exponential(t1,d))
a <- attraction(permutation(n.items=attr(d,"Size"),fixed=FALSE),decay.exponential(t2,d))
a <- attraction(permutation(1:attr(d,"Size"),fixed=TRUE),decay.exponential(t1,d))
a <- attraction(permutation(1:attr(d,"Size"),fixed=TRUE),decay.exponential(t2,d))
a

ewens(mass(fixed=FALSE),50)
ewens.pitman(mass(fixed=FALSE),discount(fixed=FALSE),50)
ewens.attraction(mass(fixed=FALSE),a)
ewens.pitman.attraction(mass(fixed=FALSE),discount(fixed=FALSE),a)
ddcrp(mass(2.3),a)

mean(nsubsets.random(ewens(mass(1,fixed=TRUE),50),100000))
mean(nsubsets.random(ewens(mass(1,fixed=FALSE),50),100000))
mean(nsubsets.random(ewens.pitman(mass(1,fixed=FALSE),discount(fixed=FALSE),50),100000))
mean(nsubsets.random(ewens.attraction(mass(fixed=FALSE),a),100000))
mean(nsubsets.random(ewens.pitman.attraction(mass(fixed=FALSE),discount(fixed=FALSE),a),100000))
mean(nsubsets.random(ewens.pitman.attraction(mass(1.0,fixed=TRUE),discount(0.05,fixed=TRUE),a),100000))

nsubsets.probability(ewens.pitman.attraction(mass(n.items=a$n.items,fixed=TRUE),discount(0.1,fixed=TRUE),a),4)
nsubsets.probability(ewens.pitman.attraction(mass(n.items=a$n.items,fixed=TRUE),discount(0.1,fixed=TRUE),a),6)

nsubsets.average(ewens.pitman.attraction(mass(2,fixed=TRUE),discount(0.1,fixed=TRUE),a))
nsubsets.variance(ewens.attraction(mass(2,fixed=TRUE),a))
tryCatch(nsubsets.variance(ewens.pitman.attraction(mass(2,fixed=TRUE),discount(0.1,fixed=TRUE),a)), error=function(e) e)



system.time(a <- sample.partition(ewens(mass(1,fixed=TRUE),10),100000))
system.time(a <- sample.partition(ewens.pitman(mass(1,fixed=TRUE),discount(fixed=TRUE),10),100000))

system.time(a <- sample.partition(ewens(mass(1,fixed=FALSE),10),100000))

system.time(a <- sample.partition(ewens.pitman(mass(1,fixed=FALSE),discount(fixed=FALSE),10),1000000,parallel=FALSE))
system.time(a <- sample.partition(ewens.pitman(mass(1,fixed=FALSE),discount(fixed=FALSE),10),1000000,parallel=TRUE))

# system.time(a <- sample.partition(ewens.pitman(mass(1,fixed=FALSE),discount(0.0,fixed=TRUE),10),100,sample.parameter=function() rnorm(1),parallel=TRUE))
# 
# attr <- attraction(10)
# attr <- attraction(d,permutation(n.items=attr(d,"Size"),fixed=FALSE),decay.exponential(distance=d,fixed=FALSE))
# attr <- attraction(d,permutation(n.items=attr(d,"Size"),fixed=FALSE),decay.exponential(distance=d,fixed=TRUE))
# attr <- attraction(d,permutation(1:attr(d,"Size"),fixed=TRUE),decay.exponential(distance=d,fixed=FALSE))
# attr <- attraction(d,permutation(1:attr(d,"Size"),fixed=TRUE),decay.exponential(distance=d,fixed=TRUE))
# system.time(a <- sample.partition(ewens.attraction(mass(1,fixed=FALSE),attr),100,sample.parameter=function() rnorm(1),as.matrix=TRUE,parallel=TRUE))
# system.time(a <- sample.partition(ewens.pitman.attraction(mass(1,fixed=FALSE),discount(fixed=FALSE),attr),100,sample.parameter=function() rnorm(1),parallel=TRUE))

