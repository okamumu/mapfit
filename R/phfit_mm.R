#' PH fitting based on three moments
#' 
#' Fits a phase-type (PH) distribution by matching the first three moments.
#' 
#' @param m1 The first moment (mean).
#' @param m2 The second moment (mean of squares).
#' @param m3 The third moment (mean of cubes).
#' @param method The moment-matching method to be used. Choices are \code{"Osogami06"} (default) or \code{"Bobbio05"}.
#' @param max.phase The maximum number of phases allowed in the \code{"Osogami06"} method.
#' @param epsilon The precision tolerance for the \code{"Osogami06"} method.
#' 
#' @return An object of class \code{GPH} representing the fitted PH distribution.
#' 
#' @note
#' The \code{"Osogami06"} method checks whether there exists a PH distribution matching the given first three moments.
#' If no feasible PH distribution exists, the \code{"Bobbio05"} method may fail or return an error.
#' 
#' @references
#' Osogami, T. and Harchol-Balter, M. (2006).  
#' Closed-form solutions for mapping general distributions to minimal PH distributions.  
#' \emph{Performance Evaluation}, \bold{63}(6), 524--552.
#' 
#' Bobbio, A., Horvath, A., and Telek, M. (2005).  
#' Matching three moments with minimal acyclic phase-type distributions.  
#' \emph{Stochastic Models}, \bold{21}(2-3), 303--326.
#' 
#' @examples
#' ## Fit a PH distribution to match moments of a Weibull(shape=2, scale=1)
#' ## Moments: (mean = 0.886227, second moment = 1.0, third moment = 1.32934)
#' result1 <- phfit.3mom(0.886227, 1.0, 1.32934)
#' 
#' ## Fit using the Bobbio05 method
#' result2 <- phfit.3mom(0.886227, 1.0, 1.32934, method = "Bobbio05")
#' 
#' ## Calculate mean
#' ph.mean(result1)
#' ph.mean(result2)
#' 
#' ## Calculate variance
#' ph.var(result1)
#' ph.var(result2)
#' 
#' ## Calculate moments up to order 5
#' ph.moment(5, result1)
#' ph.moment(5, result2)
#' 
#' @export

phfit.3mom <- function(m1, m2, m3, method = c("Osogami06", "Bobbio05"),
	max.phase = 50, epsilon = sqrt(.Machine$double.eps)) {
  method <- match.arg(method)
  switch(
    method,
    "Osogami06" = {
      res <- matching3PH(c(m1, m2, m3), epsilon=epsilon, max.phase=max.phase)
      ph(alpha=res$tau, Q=res$T, xi=res$xi)
    },
    "Bobbio05" = {
      mm.bobbio05(m1, m2, m3)
    },
    stop("The method should be chosen from Osogami06 and Bobbio05.")
  )
}
