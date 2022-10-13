library(fnets)

set.seed(123)
n <- 500
p <- 50
common <- sim.unrestricted(n, p)
idio <- sim.var(n, p)
x <- common$data * apply(idio$data, 1, sd) / apply(common$data, 1, sd) + idio$data

hl <- hl.factor.number(x, q.max = NULL, mm = floor(4 * (n / log(n))^(1 / 3)), do.plot = TRUE)
hl$q.hat
