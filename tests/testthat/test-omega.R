test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
library(fnets)
set.seed(123)
n <- 500
p <- 50
common <- sim.unrestricted(n, p)
idio <- sim.var(n, p)
x <- common$data + idio$data
out <- fnets(x, q = NULL, var.method = "lasso", do.lrpc = FALSE)

plrpc <- par.lrpc(out, x, tuning.args = list(n.folds = 1, path.length = 10, do.plot = TRUE),
                  n.cores = 1)

test_that("plrpc executes", {
  expect_equal(is.null(plrpc), FALSE)
})

out$lrpc <- plrpc
out$do.lrpc <- TRUE

test_that("plot executes", {
plot(out, type = "pc", display = "network", threshold = .05)
plot(out, type = "lrpc", display = "heatmap", threshold = .05)
})
