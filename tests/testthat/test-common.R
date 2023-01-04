
library(fnets)
set.seed(123)
n <- 500
p <- 50
common <- sim.unrestricted(n, p)
idio <- sim.var(n, p)
x <- common$data + idio$data
out <- fnets(x, q = NULL, var.order = 1, var.method = "lasso", do.lrpc = FALSE)


test_that("ic works", {
  cpre <- common.predict(out, x, h = 1, r = NULL)
})

test_that("er works", {
  cpre <- common.predict(out, x, h = 1, r = "er")
})

test_that("fc.restricted works", {
  cpre <- common.predict(out, x, h = 1, r = 2, fc.restricted = FALSE)
})
