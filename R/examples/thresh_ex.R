\dontrun{
library(fnets)

set.seed(123)
n <- 500
p <- 20
common <- sim.unrestricted(n, p)
idio <- sim.var(n, p)
x <- common$data + idio$data
out <- fnets(x,
   var.args = list(n.cores = 2)
)
# Granger-causal network
th1 <- threshold(out$idio.var$beta)
plot(th1)
print(th1)
# Partial correlations
th2 <- threshold(out$lrpc$pc)
# Long-run partial correlations
th3 <- threshold(out$lrpc$lrpc)
}
