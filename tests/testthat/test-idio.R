library(fnets)

set.seed(123)
n <- 500
p <- 20
idio <- sim.var(n, p)
x <- idio$data


test_that("var ds executes", {
  skip_on_cran()
fv <- fnets.var(x,
                center = TRUE, method = "ds", var.order = 1,
                tuning.args = list(tuning = "cv", n.folds = 1, path.length = 10),
                n.cores = 1
)
plot(fv)
plot(fv, display = 'tuning')
plot(fv, display = 'heatmap')
predict(fv)
par.lrpc(fv, n.cores = 1)
expect_equal(attr(fv, "class"), "fnets")
})
test_that("var high order", {
  skip_on_cran()
  fv <- fnets.var(x,
                  center = TRUE, method = "lasso", var.order = 5,
                  tuning.args = list(tuning = "cv", n.folds = 1, path.length = 10),
                  n.cores = 1
  )
  plot(fv)
  plot(fv, display = 'tuning')
  plot(fv, display = 'heatmap')
  predict(fv)
  predict(fv, n.ahead = 10)
  predict(fv, newdata = x, n.ahead = 10)
  expect_equal(attr(fv, "class"), "fnets")
})
test_that("threshold", {
  skip_on_cran()
  fv <- fnets.var(x,
                  center = TRUE, method = "lasso", var.order = 1, do.threshold = TRUE,
                  tuning.args = list(tuning = "cv", n.folds = 1, path.length = 10),
                  n.cores = 1
  )
  th <- threshold(fv$beta)
  th
  plot(th)
})
test_that("var cv executes", {
  skip_on_cran()
  fv <- fnets.var(x,
                  center = TRUE, method = "lasso", var.order = 1:2,
                  tuning.args = list(tuning = "cv", n.folds = 1, path.length = 10),
                  n.cores = 1
  )
  expect_equal(attr(fv, "class"), "fnets")
})
test_that("var bic executes", {
  skip_on_cran()
  fv <- fnets.var(x,
                  center = TRUE, method = "lasso", var.order = 1:2,
                  tuning.args = list(tuning = "bic", n.folds = 1, path.length = 10),
                  n.cores = 1
  )
  plot(fv, display = 'tuning')
  plot(fv, display = "heatmap", groups = rep(c(1, 2), each = p/2), group.colours = c("red", "blue"))
  expect_equal(attr(fv, "class"), "fnets")
})
