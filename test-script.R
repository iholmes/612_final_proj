library(testthat)
source("deterministic.R")
test_dir("./")


expect_that(hostMat, is_a("matrix"))

expect_that(sum(hostMat[1,]), equals(1))

expect_that(sum(nextSeed(c(1000, 2000))), equals(1000))

expect_that(length(internalLV(c(0.5, 0.5))), equals(2))

expect_that(sum(paramTest(c(0.1, 0.1), hostMat)), equals(1000))

expect_that(par1est<1, equals(TRUE))

expect_that(par2est<1, equals(TRUE))

expect_that(length(weightedDist), equals(length(kParams)/2))

expect_that(nrow(make.input(0, 1, 0, 1, 1, 10, 0.2)), equals(2160))

testPt=300
expect_that(max(apply(tmat[,5:6], 1, pairwiseDist, testPoint=testPt))<=max(apply(rbind(c(min(tmat[,5]),min(tmat[,6])), 
  c(min(tmat[,5]), max(tmat[,6])), c(max(tmat[,5]), min(tmat[,6])), c(max(tmat[,5]), max(tmat[,6]))), 1,
  pairwiseDist, testPoint=testPt)), equals(TRUE))

expect_that(checkPt(300, 4), is_a("numeric"))

expect_that(checkPt(300, 4), equals(c(0.5, 0.25, 0.2, 2)))


