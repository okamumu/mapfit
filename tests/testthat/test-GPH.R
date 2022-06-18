test_that("create GPH", {
  ph <- GPH$new(alpha=c(0.5, 0.4), Q=rbind(c(-3, 2), c(0, -2)), xi=c(1, 2))
  expect_equal(ph$alpha, c(0.5, 0.4))
})

test_that("moment", {
  ph <- ph(alpha=c(0.5, 0.4), Q=rbind(c(-3, 2), c(0, -2)), xi=c(1, 2))
  expected <- c(0.5333333,
                0.5888889,
                0.9388889,
                1.9518519,
                5.0030864,
                15.2561728,
                53.9727366,
                217.4272977,
                983.0318930,
                4930.5229767)
  result <- ph$moment(10)
  expect_equal(result, expected)
})
