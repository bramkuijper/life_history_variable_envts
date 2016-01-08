# map a gaussian distribution to a uniform
library(pracma)
library(lattice)
k <- rnorm(1000,mean=0,sd=1)

z <- .5 * erfc(-k / sqrt(2))

print(hist(z))
