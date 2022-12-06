library(fnets)

set.seed(123)
n <- 500
p <- 50
common <- sim.restricted(n, p)
idio <- sim.var(n, p)
x <- common$data * apply(idio$data, 1, sd) / apply(common$data, 1, sd) + idio$data

abc <- factor.number(x, fm.restricted = TRUE, do.plot = TRUE)
abc

er <- factor.number(x, method = "er", fm.restricted = TRUE, do.plot = TRUE)
er
