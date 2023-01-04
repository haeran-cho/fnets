# library(fnets)
# set.seed(123)
# n <- 500
# p <- 50
# common <- sim.unrestricted(n, p)
# idio <- sim.var(n, p)
# x <- common$data + idio$data
#
# out <- fnets(
#   x,
#   q = NULL,
#   var.order = 1,
#   var.method = "lasso",
#   do.threshold = TRUE,
#   do.lrpc = TRUE,
#   tuning.args = list(
#     tuning = "cv",
#     n.folds = 1,
#     path.length = 10,
#     do.plot = TRUE
#   ),
#   var.args = list(n.cores = 2)
# )
#
# test_that("fnets executes", {
#   skip_on_cran()
#   expect_equal(attr(out, "class"), "fnets")
# })
#
# test_that("predict executes", {
#   skip_on_cran()
#   pre <- predict(out, x, h = 1, common.method = "unrestricted")
#   pre <- predict(out, x, h = 1, common.method = "restricted")
# })
#
# test_that("plot executes", {
#   skip_on_cran()
#   plot(out, type = "granger", display = "network")
#   plot(out, type = "lrpc", display = "heatmap")
# })
#
# test_that("fnets.factor.model restricted executes", {
#   skip_on_cran()
#   out <- fnets.factor.model(x, fm.restricted = TRUE)
#   expect_equal(attr(out, "class"), "fm")
# })
#
# test_that("fnets.factor.model unrestricted executes", {
#   skip_on_cran()
#   out <- fnets.factor.model(x, fm.restricted = FALSE)
#   expect_equal(attr(out, "class"), "fm")
# })
