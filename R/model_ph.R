#' General phase-type distribution
#'
#' A continuous distribution dominated by a continuous-time Markov chain.
#' A random time is given by an absorbing time.
#' 
GPH <- R6::R6Class("GPH",
    private = list(
      matclass = "dgCMatrix",
      P = NULL
    ),
    public = list(
      #' @field alpha A vector of initial probability
      alpha = NULL,
      
      #' @field Q A matrix of infinitesimal generator
      Q = NULL,
      
      #' @field xi An exit vector of transient states to the absorbing state
      xi = NULL,
    
      #' @description
      #' Create a GPH
      #' @param alpha A vector of initial probability
      #' @param Q An infinitesimal generator
      #' @param xi An exit rate vector
      #' @return An instance of GPH
      initialize = function(alpha, Q, xi) {
        self$alpha <- alpha
        self$Q <- as(Q, private$matclass)
        self$xi <- xi
        private$P <- as(Q, private$matclass)
      },
      
      #' @description 
      #' The number of phases
      #' @return The number of phases
      size = function() {
        length(self$alpha)
      },

      #' @description 
      #' Degrees of freedom
      #' @return The degrees of freedom
      df = function() {
        zero <- 1.0e-8
        sum(self$alpha > zero) - 1 + sum(abs(self$Q) > zero) +
          sum(self$xi > zero) - size + sum(abs(Matrix::diag(selfQ)) < zero)
      },
      
      #' @description 
      #' Moments of GPH
      #' @param k A value to indicate the degrees of moments. k-th moment
      #' @param ... Others
      #' @return A vector of moments from 1st to k-th moments
      moment = function(k, ...) {
        tmp <- self$alpha
        A <- t(-as.matrix(self$Q))
        tmp2 <- 1.0
        res <- numeric(0)
        for (i in 1:k) {
          tmp <- solve(a=A, b=tmp)
          tmp2 <- tmp2 * i
          res <- c(res, tmp2 * sum(tmp))
        }
        res
      },
      
      #' @description 
      #' Print
      #' @param ... Others
      print = function(...) {
        cat(gettextf("Size : %d\n", self$size()))
        cat("Initial : ", self$alpha, "\n")
        cat("Exit    : ", self$xi, "\n")
        cat("Infinitesimal generator : \n")
        print(self$Q)
      },
      
      #' @description 
      #' PDF
      #' @param x A vector of points
      #' @param poisson.eps A value of tolerance error for uniformization
      #' @param ufactor A value of uniformization factor
      #' @param ... Others
      #' @return A vector of densities.
      #' @note
      #' This function provides the values of p.d.f. for PH distribution with
      #' the uniformization technique.
      pdf = function(x, poisson.eps = 1.0e-8, ufactor = 1.01, ...) {
        inv <- order(order(x))
        sx <- sort(x)
        dt <- c(sx[1], diff(sx))
        y <- phase_dist_pdf(dt, max(dt), self$alpha, self$Q, self$xi,
                            poisson.eps, ufactor, private$P)
        y[inv]
      },
      
      #' @description 
      #' CDF
      #' @param x A vector of points
      #' @param poisson.eps A value of tolerance error for uniformization
      #' @param ufactor A value of uniformization factor
      #' @param ... Others
      #' @return A vector of probabilities
      #' @note
      #' This function provides the values of c.d.f. for PH distribution with
      #' the uniformization technique.
      cdf = function(x, poisson.eps = 1.0e-8, ufactor = 1.01, ...) {
        inv <- order(order(x))
        sx <- sort(x)
        dt <- c(sx[1], diff(sx))
        y <- phase_dist_ccdf(dt, max(dt), self$alpha, self$Q,
                             poisson.eps, ufactor, private$P)
        1 - y[inv]
      },

      #' @description 
      #' Complementary CDF
      #' @param x A vector of points
      #' @param poisson.eps A value of tolerance error for uniformization
      #' @param ufactor A value of uniformization factor
      #' @param ... Others
      #' @return A vector of probabilies
      #' @note
      #' This function provides the values of complementary c.d.f. for
      #' PH distribution with the uniformization technique.
      ccdf = function(x, poisson.eps = 1.0e-8, ufactor = 1.01, ...) {
        inv <- order(order(x))
        sx <- sort(x)
        dt <- c(sx[1], diff(sx))
        y <- phase_dist_ccdf(dt, max(dt), self$alpha, self$Q,
                            poisson.eps, ufactor, private$P)
        y[inv]
      },
      
      #' @description 
      #' Make a sample
      #' @param ... Others
      #' @return A sample of GPH
      sample = function(...) {
        s <- which(as.vector(rmultinom(n=1, size=1, prob=c(self$alpha,0)))==1)
        t <- 0
        while (s != self$size()+1) {
          x <- c(self$Q[s,], self$xi[s])
          r <- -x[s]
          p <- x / r
          p[s] <- p[s] + 1
          t <- t + rexp(n=1, rate=r)
          s <- which(as.vector(rmultinom(n=1, size=1, prob=p))==1)
        }
        t
      }
    )
)

#' Create GPH distribution
#' 
#' Create an instance of GPH
#' 
#' @param size An integer for the number of phases
#' @param alpha A vector of initial probability
#' @param Q An infinitesimal generator
#' @param xi An exit rate vector
#' @return An instance of GPH
#' 
#' @note
#' This function can omit several patterns of arguments. For example, `ph(5)`
#' omit the arguments `alpha`, `Q` and `xi`. In this case, the default values are
#' assigned to them.
#' 
#' @examples
#' ## create a PH (full matrix) with 5 phases
#' (param1 <- ph(5))
#'
#' ## create a PH (full matrix) with 5 phases
#' (param1 <- ph(size=5))
#'
#' ## create a PH with specific parameters
#' (param2 <- ph(alpha=c(1,0,0),
#'               Q=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-1)),
#'               xi=c(2,2,0)))
#'
#' @export

ph <- function(size, alpha, Q, xi) {
  if (missing(size)) {
    if (missing(alpha) || missing(Q) || missing(xi)) {
      stop("alpha, Q and xi are needed.")
    } else {
      size <- length(alpha)
    }
  } else {
    if (!missing(alpha) && !missing(Q) && !missing(xi)) {
      warning("size is ignored.")
      size <- length(alpha)
    } else {
      if (!missing(alpha) || !missing(Q) || !missing(xi)) {
        warning("alpha, Q and xi are ignored.")
      }
      alpha <- rep(1.0/size, size)
      Q <- matrix(1.0, size, size)
      diag(Q) <- rep(-size, size)
      xi <- rep(1.0, size)
    }
  }
  if (missing(xi)) {
    xi = -apply(Q, 1, sum)
  }
  GPH$new(alpha=alpha, Q=Q, xi=xi)
}


#' Create a bi-diagonal PH distribution
#' 
#' Create an instance of bi-diagonal PH distribution.
#' 
#' @param size An integer for the number of phases
#' @return An instance of bi-diagonal PH distribution
#' 
#' @note
#' Bi-diagonal PH distribution is the PH distribution whose infinitesimal
#' generator is given by a upper bi-diagonal matrix. This is similar to
#' canonical form 1. But there is no restriction on the order for diagonal
#' elements. 
#' 
#' @examples
#' ## create a bidiagonal PH with 5 phases
#' (param1 <- ph.bidiag(5))
#'
#' @export

ph.bidiag <- function(size) {
  if (size <= 1) {
    ph(size)
  } else {
    alpha <- rep(1/size,size)
    xi <- rep(0, size)
    Q <- matrix(0, size, size)
    for (i in 1:(size-1)) {
      Q[i,i] <- -1
      Q[i,i+1] <- 1
    }
    Q[size,size] <- -1
    xi[size] <- 1
    ph(alpha=alpha, Q=Q, xi=xi)
  }
}

#' Create a tri-diagonal PH distribution
#' 
#' Create an instance of tri-diagonal PH distribution.
#' 
#' @param size An integer for the number of phases
#' @return An instance of tri-diagonal PH distribution
#' 
#' @note
#' Tri-diagonal PH distribution is the PH distribution whose infinitesimal
#' generator is given by a tri-diagonal matrix (band matrix).
#' 
#' @examples
#' ## create a tridiagonal PH with 5 phases
#' (param1 <- ph.tridiag(5))
#'
#' @export

ph.tridiag <- function(size) {
  if (size <= 2) {
    ph(size)
  } else {
    alpha <- rep(1/size,size)
    xi <- rep(0, size)
    Q <- matrix(0, size, size)
    Q[1,1] <- -1
    Q[1,2] <- 1
    for (i in 2:(size-1)) {
      Q[i,i] <- -2
      Q[i,i-1] <- 1
      Q[i,i+1] <- 1
    }
    Q[size,size-1] <- 1
    Q[size,size] <- -2
    xi[size] <- 1
    ph(alpha=alpha, Q=Q, xi=xi)
  }
}

#' Probability density function of PH distribution
#' 
#' Compute the probability density function (p.d.f.) for a given PH distribution
#' 
#' @param x A numeric vector of quantiles.
#' @param ph An instance of PH distribution.
#' @param log logical; if TRUE, densities y are returned as log(y)
#' @param ... Others.
#' @return A vector of densities.
#' @examples 
#' ## create PH instance (ex. full PH with 5 phases)
#' (phdist <- ph(5))
#' 
#' x <- runif(10) # points of quantiles
#' dphase(x, ph)
#' 
#' @export

dphase <- function(x, ph, log = FALSE, ...) {
  y <- ph$pdf(x, ...)
  if (log) {
    log(y)
  } else {
    y
  }
}

#' Distribution function of PH distribution
#' 
#' Compute the cumulative distribution fuction (c.d.f.) for a given PH distribution
#' 
#' @param q A numeric vector of quantiles.
#' @param ph An instance of PH distribution.
#' @param lower.tail logical; if TRUE, probabilities are P(X <= x), otherwise P(X > x).
#' @param log.p logical; if TRUE, probabilities p are returned as log(p).
#' @param ... Others
#' @return A vector of densities
#' @examples 
#' ## create PH instance (ex. full PH with 5 phases)
#' (phdist <- ph(5))
#' 
#' x <- runif(10) # points of quantiles
#' pphase(x, ph)
#' 
#' @export

pphase <- function(q, ph, lower.tail = TRUE, log.p = FALSE, ...) {
  if (lower.tail) {
    y <- ph$cdf(q, ...)
    if (log.p) {
      log(y)
    } else {
      y
    }
  } else {
    y <- ph$ccdf(q, ...)
    if (log.p) {
      log(y)
    } else {
      y
    }
  }
}

#' Sampling of PH distributions
#' 
#' Generate a sample from a given PH distribution.
#' 
#' @param n An integer of the number of samples.
#' @param ph An instance of PH distribution.
#' @param ... Ohters
#' @return A vector of samples.
#' @examples 
#' ## create PH instance (ex. full PH with 5 phases)
#' (phdist <- ph(5))
#'
#' ## Generate 10 samples
#' rphase(10, ph)
#' 
#' @export

rphase <- function(n, ph, ...) {
  sapply(1:n, function(i) ph$sample(...))
}

#' Moments of PH distribution
#' 
#' Generate moments up to k-th moments for a given PH distribution.
#' 
#' @param k An integer for the moments to be computed
#' @param ph An instance of PH distribution
#' @param ... Others
#' @return A vector of moments
#' @examples
#' ## create PH instance (ex. full PH with 5 phases)
#' (phdist <- ph(5))
#'
#' ph.moment(10, ph)
#' 
#' @export

ph.moment <- function(k, ph, ...) {
  ph$moment(k, ...)
}

#' Mean of PH distribution
#' 
#' Compute the mean of a given PH distribution.
#' 
#' @param ph An instance of PH distribution
#' @param ... Others
#' @return A value of mean
#' @examples
#' ## create PH instance (ex. full PH with 5 phases)
#' (phdist <- ph(5))
#'
#' ph.mean(ph)
#' 
#' @export

ph.mean <- function(ph, ...) {
  ph$moment(1, ...)
}

#' Variance of PH distribution
#' 
#' Compute the variance of a given PH distribution.
#' 
#' @param ph An instance of PH distribution
#' @param ... Others
#' @return A value of variance
#' @examples
#' ## create PH instance (ex. full PH with 5 phases)
#' (phdist <- ph(5))
#'
#' ph.var(ph)
#' 
#' @export

ph.var <- function(ph, ...) {
  x <- ph$moment(2, ...)
  x[2] - x[1]^2
}

#' setMethod("emfit.init", signature(model = "ph"),
#'   function(model, data, verbose = list(), ...) {
#'     ph.param.random(size=model@size, data=data,
#'       skelpi=model@alpha, skelQ=as.matrix(model@Q),
#'       skelxi=model@xi, verbose=verbose, class=class(model@Q))
#'   }
#' )
#' 
#' ph.param.random <- function(size, data, skelpi, skelQ, skelxi, verbose, class) {
#'   if (missing(size)) size <- length(skelpi)
#'   if (missing(skelpi)) skelpi <- rep(1, size)
#'   if (missing(skelQ)) skelQ <- matrix(1, size, size)
#'   if (missing(skelxi)) skelxi <- rep(1, size)
#' 
#'   mean <- mapfit.mean(data)
#' 
#'   alpha <- skelpi * runif(size)
#'   alpha <- alpha / sum(alpha)
#' 
#'   diag(skelQ) <- 0
#'   Q <- skelQ * matrix(runif(size*size), size, size)
#'   xi <- skelxi * runif(size)
#'   diag(Q) <- -(apply(Q, 1, sum) + xi)
#' 
#'   p <- ph(alpha=alpha, Q=Q, xi=xi, class=class)
#'   m <- ph.mean(p) / mean
#'   p@Q <- as(as.matrix(p@Q * m), class)
#'   p@xi <- p@xi * m
#'   p
#' }
#' 
#' #### estep
#' 
#' setMethod("emfit.estep", signature(model = "ph", data = "phdata.wtime"),
#'   function(model, data, ufact = 1.01, eps = sqrt(.Machine$double.eps), ...) {
#'     res <- .Call('phfit_estep_gen_wtime', PACKAGE='mapfit', model, data, eps, ufact)
#'     list(eres=list(etotal=res[[1]], eb=res[[2]], ey=res[[3]], ez=res[[4]], en=res[[5]]), llf=res[[6]])
#'   })
#' 
#' setMethod("emfit.estep", signature(model = "ph", data = "phdata.group"),
#'   function(model, data, ufact = 1.01, eps = sqrt(.Machine$double.eps), ...) {
#' 
#'   data@data$instant[is.na(data@data$counts)] <- 0
#'   data@data$counts[is.na(data@data$counts)] <- -1
#'   l <- data@size
#'   if (is.infinite(data@data$time[l])) {
#'     gdatlast <- data@data$counts[l]
#'     data@data <- data@data[-l,]
#'     data@size <- data@size - 1
#'   } else {
#'     gdatlast <- 0
#'   }
#' 
#'   ba <- msolve(alpha=1.0, A=-as.matrix(model@Q), x=model@alpha, transpose=TRUE)
#' 
#'   res <- .Call('phfit_estep_gen_group', PACKAGE='mapfit', model, ba, data, gdatlast, eps, ufact)
#'   list(eres=list(etotal=res[[1]], eb=res[[2]], ey=res[[3]], ez=res[[4]], en=res[[5]]), llf=res[[6]])
#'   })
#' 
#' #### mstep
#' 
#' setMethod("emfit.mstep", signature(model = "ph"),
#'   function(model, eres, data, ...) {
#'     res <- .Call('phfit_mstep_gen', PACKAGE='mapfit', model, eres, data)
#'     model@alpha <- res[[1]]
#'     model@xi <- res[[2]]
#'     model@Q@x <- res[[3]]
#'     model
#'   })
#' 
