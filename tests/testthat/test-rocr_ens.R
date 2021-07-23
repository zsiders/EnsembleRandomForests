context("Test the ROCR function")

test_that("roc generates values", {
  x <- rocr_ens(simData$samples$prob, factor(simData$samples$obs,levels=c(0,1)))
  expect_is(x$auc, "numeric")
})
#> Test passed 🎉
