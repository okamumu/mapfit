#' Canonical phase-type distribution
#'
#' A continuous distribution dominated by a continuous-time Markov chain.
#' A random time is given by an absorbing time. In the CF1 (canonical form 1),
#' the infinitesimal generator is given by a bi-diagonal matrix, and whose order
#' is determiend by the ascending order.
#' 
CF1 <- R6::R6Class(
  "CF1",
  inherit = GPH,
  public = list(
    #' @field rate A vector of rates
    rate = NULL,
    
    #' @description
    #' Create a CF1
    #' @param alpha A vector of initial probability
    #' @param rate A vector of rates
    #' @return An instance of CF1
    initialize = function(alpha, rate) {
      phase_cf1_sort(alpha, rate)
      self$rate <- rate
      size <- length(alpha)
      xi <- numeric(size)
      xi[size] <- rate[size]
      if (size >= 2) {
        i <- c(1:size, 1:(size-1))
        j <- c(1:size, 2:size)
        x <- c(-rate, rate[1:(size-1)])
        Q <- sparseMatrix(i=i, j=j, x=x)
      } else {
        Q <- matrix(-rate[1],1,1)
      }
      super$initialize(alpha, Q, xi)
    },

    #' @description 
    #' Print
    #' @param ... Others
    print = function(...) {
      cat(gettextf("Size : %d\n", self$size()))
      cat("Initial : ", self$alpha, "\n")
      cat("Rate    : ", self$rate, "\n")
    }
  )
)

#' Create CF1
#' 
#' Create an instance of CF1.
#' 
#' @param size An integer of the number of phases
#' @param alpha A vector of initial probabilities
#' @param rate A vector of rates
#' @return An instance of CF1.
#' @examples
#' ## create a CF1 with 5 phases
#' (param1 <- cf1(5))
#' 
#' ## create a CF1 with 5 phases
#' (param1 <- cf1(size=5))
#' 
#' ## create a CF1 with specific parameters
#' (param2 <- cf1(alpha=c(1,0,0), rate=c(1.0,2.0,3.0)))
#' 
#' @export

cf1 <- function(size, alpha, rate) {
  if (missing(size)) {
    if (missing(alpha) || missing(rate)) {
      stop("alpha and rate are needed.")
    } else {
      size <- length(alpha)
    }
  } else {
    if (!missing(alpha) && !missing(rate)) {
      warning("size is ignored.")
      size <- length(alpha)
    } else {
      if (!missing(alpha) || !missing(rate)) {
        warning("alpha and rate are ignored.")
      }
      alpha <- rep(1.0/size, size)
      rate <- rep(1.0, size)      
    }
  }
  CF1$new(alpha, rate)
}
