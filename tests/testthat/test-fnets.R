# commented out for CRAN
library(fnets)
set.seed(123)
n <- 500
p <- 20
common <- sim.unrestricted(n, p)
idio <- sim.var(n, p)
x <- common$data + idio$data

out <- fnets(
  x,
  q = NULL,
  var.order = 1,
  var.method = "lasso",
  do.threshold = TRUE,
  do.lrpc = TRUE,
  tuning.args = list(
    tuning = "cv",
    n.folds = 1,
    path.length = 10
  ),
  var.args = list(n.cores = 2)
)

test_that("fnets executes", {
  skip_on_cran()
  expect_equal(attr(out, "class"), "fnets")
})



test_that("predict executes", {
  skip_on_cran()
  pre <- predict(out, common.method = "unrestricted")
  pre <- predict(out, common.method = "restricted")
  pre <- predict(out, common.method = "unrestricted", n.ahead = 10)
})

test_that("plot executes", {
  skip_on_cran()
  plot(out, type = "granger", display = "network")
  plot(out, type = "lrpc", display = "network")
  plot(out, type = "pc", display = "network")
  plot(out, type = "granger", display = "heatmap")
  plot(out, type = "lrpc", display = "heatmap")
  plot(out, type = "pc", display = "heatmap")
  plot(out, display = "tuning")
})

test_that("network executes", {
  skip_on_cran()
  network(out, type = "granger")
  network(out, type = "pc")
  network(out, type = "lrpc")
})

test_that("print executes", {
  skip_on_cran()
  print(out)
})

test_that("fnets.factor.model restricted executes", {
  skip_on_cran()
  out <- fnets.factor.model(x, fm.restricted = TRUE)
  expect_equal(attr(out, "class"), "fm")
})

test_that("fnets.factor.model unrestricted executes", {
  skip_on_cran()
  out <- fnets.factor.model(x, fm.restricted = FALSE)
  expect_equal(attr(out, "class"), "fm")
})


test_that("q=0", {
  out <- fnets(
    x,
    q = 0,
    var.order = 1,
    var.method = "lasso",
    do.threshold = TRUE,
    do.lrpc = TRUE,
    tuning.args = list(
      tuning = "cv",
      n.folds = 1,
      path.length = 10
    ),
    var.args = list(n.cores = 2)
  )
  predict(out, n.ahead = 10)
  predict(out, newdata = x, n.ahead = 10)
})
