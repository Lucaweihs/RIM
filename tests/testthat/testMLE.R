library(RIM)
context("Testing MLE")
source("HelperTestFunctions.R")

testMLEWithTree = function(a, samples = 40000, tol = 2 * 10^-2) {
  A = rRIM(samples, a)
  expected = a
  found = thetaMLERIM(a, averageDiscMatrix(A))
  expect_equal(expected[,-4], found[,-4])
  expect_true(all(abs(expected[,4] - found[,4]) <  tol))
}

testStructByDP = function(trueTree, refRanking, samples = 40000,
                          tol = 2 * 10^-2, makeCanonical = F) {
  A = rRIM(samples, trueTree)
  expected = trueTree
  found = structByDP(averageDiscMatrix(A), refRanking, makeCanonical)
  expect_equal(expected[,-4], found[,-4])
  expect_true(all(abs(expected[,4] - found[,4]) <  tol))
}

test_that("Test thetaMLERIM", {
  set.seed(2039)
  a = treeListToMat(list(
    c(0, 0, 1, 0, 0)
  ))
  testMLEWithTree(a)

  a = treeListToMat(list(
    c(1, 2, 0, -.9, 0),
    c(0, 0, 1, 0, 0),
    c(0, 0, 1, 0, 1)
  ))
  testMLEWithTree(a)

  a = treeListToMat(list(
    c(1, 2, 0, -.1, 0),
    c(3, 4, 0, .8, 0),
    c(5, 6, 0, 1.6, 0),
    c(0, 0, 1, 0, 0),
    c(0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 2),
    c(0, 0, 1, 0, 3)
  ))
  testMLEWithTree(a)

  a = treeListToMat(list(
    c(1, 2, 0, 0, 0),
    c(3, 4, 0, 0, 0),
    c(5, 6, 0, 0, 0),
    c(0, 0, 1, 0, 0),
    c(0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 2),
    c(0, 0, 1, 0, 3)
  ))
  testMLEWithTree(a)

  a = treeListToMat(list(
    c(1, 2, 0, .1, 0),
    c(3, 4, 0, -.2, 0),
    c(5, 6, 0, 0, 0),
    c(7, 8, 0, .3, 0),
    c(0, 0, 1, 0, 0),
    c(9, 10, 0, .15, 0),
    c(0, 0, 1, 0, 5),
    c(0, 0, 1, 0, 3),
    c(0, 0, 1, 0, 2),
    c(0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 4)
  ))
  testMLEWithTree(a, tol = 3 * 10^-2)

  a = treeListToMat(list(
    c(1, 2, 0, 0, 0),
    c(3, 4, 0, 0, 0),
    c(5, 6, 0, 0, 0),
    c(7, 8, 0, 0, 0),
    c(0, 0, 1, 0, 2),
    c(9, 10, 0, 0, 0),
    c(0, 0, 1, 0, 0),
    c(0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 3),
    c(0, 0, 1, 0, 4),
    c(0, 0, 1, 0, 5)
  ))
  testMLEWithTree(a, tol = 3 * 10^-2)
})

test_that("Test structByDP", {
  set.seed(2039)
  a = treeListToMat(list(
    c(0, 0, 1, 0, 0)
  ))
  testStructByDP(a, 0)

  a = treeListToMat(list(
    c(1, 2, 0, -.9, 0),
    c(0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 0)
  ))
  testStructByDP(a, c(1,0))

  a = treeListToMat(list(
    c(1, 2, 0, -.1, 0),
    c(3, 4, 0, .8, 0),
    c(5, 6, 0, 1.6, 0),
    c(0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 0),
    c(0, 0, 1, 0, 3),
    c(0, 0, 1, 0, 2)
  ))
  testStructByDP(a, c(1,0,3,2))

  a = treeListToMat(list(
    c(1, 2, 0, 1.1, 0),
    c(3, 4, 0, -1, 0),
    c(5, 6, 0, .2, 0),
    c(0, 0, 1, 0, 0),
    c(0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 2),
    c(0, 0, 1, 0, 3)
  ))
  testStructByDP(a, 0:3)

  a = treeListToMat(list(
    c(1, 2, 0, .1, 0),
    c(3, 4, 0, .2, 0),
    c(5, 6, 0, 1, 0),
    c(7, 8, 0, .3, 0),
    c(0, 0, 1, 0, 0),
    c(9, 10, 0, .15, 0),
    c(0, 0, 1, 0, 5),
    c(0, 0, 1, 0, 3),
    c(0, 0, 1, 0, 2),
    c(0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 4)
  ))
  testStructByDP(a, c(5,4,1, 0,3,2), makeCanonical = T)

  a = treeListToMat(list(
    c(1, 2, 0, .35, 0),
    c(3, 4, 0, .39, 0),
    c(5, 6, 0, -.1, 0),
    c(7, 8, 0, -.15, 0),
    c(0, 0, 1, 0, 2),
    c(9, 10, 0, 1, 0),
    c(0, 0, 1, 0, 0),
    c(0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 3),
    c(0, 0, 1, 0, 4),
    c(0, 0, 1, 0, 5)
  ))
  testStructByDP(a, c(1,3,2, 4,5,0))
})
