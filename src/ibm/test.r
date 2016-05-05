# map a gaussian distribution to a uniform
library(pracma)
library(lattice)


t <- 1:200

p <- 0.05

et1 <- .5 + .5 * sin(t)

et2 <-  runif(t)

print(xyplot(et1 + et2 ~ t,type="l"))
