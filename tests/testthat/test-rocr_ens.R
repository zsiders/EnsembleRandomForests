context("Test the ROCR function")

test_that("roc generates key performance values", {
  x <- rocr_ens(simData$samples$prob, factor(simData$samples$obs,levels=c(0,1)))
  expect_type(x$auc, "double")
  expect_type(x$tss, "double")
  expect_type(x$rmse, "double")
})
#> Test passed 🎉
test_that("roc generates performance class", {
  x <- rocr_ens(simData$samples$prob, factor(simData$samples$obs,levels=c(0,1)))
  expect_s4_class(x$tpr, "performance")
  expect_s4_class(x$tnr, "performance")
  expect_s4_class(x$phi, "performance")
  expect_s4_class(x$acc, "performance")
  expect_s4_class(x$err, "performance")
  expect_s4_class(x$fpr, "performance")
  expect_s4_class(x$sens, "performance")
  expect_s4_class(x$spec, "performance")
})
#> Test passed 🎉