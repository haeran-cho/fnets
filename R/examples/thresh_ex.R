\dontrun{
library(fnets)
out <- fnets(data.unrestricted,
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
