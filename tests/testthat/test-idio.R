library(fnets)

set.seed(123)
n <- 500
p <- 50
idio <- sim.var(n, p)
x <- idio$data

test_that("var ds executes", {
  skip_on_cran()
fv <- fnets.var(x,
                center = TRUE, method = "ds", var.order = 1,
                tuning.args = list(tuning = "cv", n.folds = 1, path.length = 10, do.plot = TRUE),
                n.cores = 1
)
expect_equal(attr(fv, "class"), "fnets")
})
test_that("var bic executes", {
  skip_on_cran()
  fv <- fnets.var(x,
                  center = TRUE, method = "lasso", var.order = 1,
                  tuning.args = list(tuning = "bic", n.folds = 1, path.length = 10, do.plot = TRUE),
                  n.cores = 1
  )
  expect_equal(attr(fv, "class"), "fnets")
})
