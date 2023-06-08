library(fnets)
set.seed(123)
n <- 500
p <- 50
common <- sim.unrestricted(n, p)
idio <- sim.var(n, p)
x <- common$data + idio$data
out <- fnets(x, q = NULL, var.method = "lasso", do.lrpc = FALSE)

plrpc <- par.lrpc(out, tuning.args = list(n.folds = 1, path.length = 10),
                  n.cores = 1)

test_that("plrpc executes", {
  expect_equal(is.null(plrpc), FALSE)
})

out$lrpc <- plrpc
out$do.lrpc <- TRUE

test_that("plot executes", {
plot(out, type = "pc", display = "network")
plot(out, type = "lrpc", display = "heatmap")
})
