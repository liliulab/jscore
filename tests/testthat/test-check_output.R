#
library(jScore)


#Test whether the input truth and pred lengths are different
test_that("Output is numeric", {
  pred <- c(1,1,1,1,1,3,4,5)
  truth <- c(1,1,1,3,3,3,3,3)
  j <- jscore(truth, pred)
  expect_is(j, "numeric")
})

