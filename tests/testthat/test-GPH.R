test_that("create GPH", {
  ph <- GPH$new(alpha=c(0.6, 0.4), Q=rbind(c(-3, 2), c(0, -2)), xi=c(1, 2))
  expect_equal(ph$alpha(), c(0.6, 0.4))
})

test_that("density", {
  tt <- seq(0, 10, length.out=10)
  alpha <- c(0.6, 0.4)
  Q <- rbind(c(-3, 2), c(0, -2))
  xi <- c(1, 2)
  
  f <- function(t, alpha, xi, Q) {
    (alpha %*% expm(Q*t) %*% xi)[1]
  }
  
  result <- dphase(tt, ph(alpha=alpha, Q=Q, xi=xi))
  expected <- sapply(tt, f, alpha=alpha, Q=Q, xi=xi)
  expect_equal(result, expected)
})

test_that("ccdf", {
  tt <- seq(0, 10, length.out=10)
  alpha <- c(0.6, 0.4)
  Q <- rbind(c(-3, 2), c(0, -2))
  xi <- c(1, 2)
  
  f <- function(t, alpha, xi, Q) {
    vone <- rep(1, length(alpha))
    (alpha %*% expm(Q*t) %*% vone)[1]
  }
  
  result <- pphase(tt, ph(alpha=alpha, Q=Q, xi=xi), lower.tail = F)
  expected <- sapply(tt, f, alpha=alpha, Q=Q, xi=xi)
  expect_equal(result, expected)
})

test_that("estimate_point", {
  RNGkind(kind = "Mersenne-Twister")
  set.seed(1234)
  wsample <- rweibull(100, shape=2, scale=1)
  tres <- system.time(result <- phfit.point(ph=ph(5), x=wsample))
  print(result)
  print(tres)
})

test_that("estimate_group", {
  RNGkind(kind = "Mersenne-Twister")
  set.seed(1234)
  wsample <- rweibull(100, shape=2, scale=1)
  h.res <- hist(wsample, breaks="fd", plot=FALSE)
  tres <- system.time(result <- phfit.group(ph=ph(5), counts=h.res$counts, breaks=h.res$breaks))
  print(result)
  print(tres)
})
