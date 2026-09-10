test_that("density", {
  library(Matrix)
  tt <- seq(0, 10, length.out=10)
  alpha <- c(0.6, 0.4)
  rate <- c(3, 4)

  p <- cf1(alpha=alpha, rate=rate)
  f <- function(t, alpha, xi, Q) {
    (alpha %*% Matrix::expm(Q*t) %*% xi)[1]
  }
  
  result <- dphase(tt, p)
  expected <- sapply(tt, f, alpha=p$alpha(), Q=p$Q(), xi=p$xi())
  expect_equal(result, expected)
})

test_that("estimate_point", {
  RNGkind(kind = "Mersenne-Twister")
  set.seed(1234)
  wsample <- rweibull(100, shape=2, scale=1)
  tres <- system.time(result <- phfit.point(ph=cf1(5), x=wsample))
  print(result)
  print(tres)
})

test_that("estimate_group", {
  RNGkind(kind = "Mersenne-Twister")
  set.seed(1234)
  wsample <- rweibull(100, shape=2, scale=1)
  h.res <- hist(wsample, breaks="fd", plot=FALSE)
  tres <- system.time(result <- phfit.group(ph=cf1(5), counts=h.res$counts, breaks=h.res$breaks))
  print(result)
  print(tres)
})

test_that("cf1 does not destroy its arguments", {
  ## phase_cf1_sort sorts in place through Rcpp; cf1() must not touch the caller
  rate <- c(3, 1, 2)
  alpha <- c(0.2, 0.3, 0.5)
  invisible(cf1(alpha=alpha, rate=rate))
  expect_equal(rate, c(3, 1, 2))
  expect_equal(alpha, c(0.2, 0.3, 0.5))
})

test_that("the model is consistent after estimation", {
  RNGkind(kind = "Mersenne-Twister")
  set.seed(1234)
  wsample <- rweibull(n=100, shape=2, scale=1)
  m <- phfit.point(ph=cf1(3), x=wsample, cf1.verbose=FALSE)$model

  ## the exit rate must follow the estimated rates
  expect_equal(m$xi()[m$size()], m$rate()[m$size()])
  expect_equal(m$xi()[-m$size()], rep(0, m$size()-1))

  ## and hence dphase must be a probability density function
  expect_equal(integrate(function(x) dphase(x, m), 0, 50)$value, 1, tolerance=1e-6)

  ## rebuilding the model from the estimated parameters gives the same density
  tt <- c(0.3, 0.8, 1.5)
  expect_equal(dphase(tt, m),
               dphase(tt, cf1(alpha=as.vector(m$alpha()), rate=as.vector(m$rate()))))
})

test_that("an unsupported data class is rejected", {
  expect_error(cf1(3)$emfit(data.frame.map.time(time=c(1,2,3)), emoptions()),
               "data class")
})
