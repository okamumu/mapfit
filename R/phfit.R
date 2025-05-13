#' PH fitting with point data
#' 
#' Fits a phase-type (PH) distribution to point data via maximum likelihood estimation.
#' 
#' @param ph An object of R6 class for PH distributions. The estimation algorithm is selected based on this class.
#' @param x A numeric vector of observed point data (e.g., event times).
#' @param weights A numeric vector of weights associated with each point. If omitted, equal weights are assumed.
#' @param ... Additional options for the fitting algorithm.
#' 
#' @return A list of class \code{phfit.result} with the following components:
#' \item{model}{An object for the estimated PH distribution.}
#' \item{llf}{The maximized log-likelihood value.}
#' \item{df}{The degrees of freedom of the fitted model.}
#' \item{aic}{The Akaike information criterion (AIC) value.}
#' \item{iter}{The number of iterations performed.}
#' \item{convergence}{A logical value indicating whether the algorithm converged.}
#' \item{ctime}{The computation time (user time).}
#' \item{data}{An object containing the input data.}
#' \item{aerror}{The absolute error of the log-likelihood at the final iteration.}
#' \item{rerror}{The relative error of the log-likelihood at the final iteration.}
#' \item{options}{A list of options used for the fitting.}
#' \item{call}{The matched function call.}
#' 
#' @examples
#' ## Generate a sample
#' wsample <- rweibull(n = 100, shape = 2, scale = 1)
#' 
#' ## Fit a general PH distribution
#' result1 <- phfit.point(ph = ph(2), x = wsample)
#' 
#' ## Fit a Canonical Form 1 (CF1) distribution
#' result2 <- phfit.point(ph = cf1(2), x = wsample)
#' 
#' ## Fit a hyper-Erlang distribution
#' result3 <- phfit.point(ph = herlang(3), x = wsample)
#' 
#' ## Calculate statistics
#' ph.mean(result1$model)
#' ph.var(result2$model)
#' ph.moment(5, result3$model)
#' 
#' @export

phfit.point <- function(ph, x, weights, ...) {
  call <- match.call()
  ph <- ph$copy()
  
  options <- emoptions()
  con <- list(...)
  nmsC <- names(options)
  options[(namc <- names(con))] <- con
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  data <- data.frame.phase.time(x=x, weights=weights)
  if (options$initialize == TRUE) {
    ph$init(data, options)
  }
  tres <- system.time(result <- ph$emfit(data, options, ...))
  result <- c(result, list(model=ph, aic=-2*(result$llf - ph$df()), df=ph$df(),
                           data=data, ctime=tres[1], options=options, call=call))
  class(result) <- "phfit.result"
  result
}

#' PH fitting with grouped data
#' 
#' Fits a phase-type (PH) distribution to grouped data via maximum likelihood estimation.
#' 
#' @param ph An object of R6 class for PH distributions. The estimation algorithm is selected based on this class.
#' @param counts A numeric vector of counts in each interval.
#' @param breaks A numeric vector of boundaries for the intervals. Equivalent to \code{c(0, cumsum(intervals))}.
#' If omitted, defaults to \code{0:length(counts)}.
#' @param intervals A numeric vector of interval lengths. Equivalent to \code{diff(breaks)}.
#' If omitted, defaults to \code{rep(1, length(counts))}.
#' @param instants A numeric vector indicating whether an instantaneous sample is drawn at the end of each interval.
#' If \code{instants[i] = 1}, a sample is drawn at the right boundary of the \code{i}-th interval.
#' If omitted, defaults to \code{rep(0L, length(counts))} (no instantaneous samples).
#' @param ... Additional options for the EM algorithm.
#' 
#' @return A list of class \code{phfit.result} with the following components:
#' \item{model}{An object for the estimated PH distribution.}
#' \item{llf}{The maximized log-likelihood value.}
#' \item{df}{The degrees of freedom of the fitted model.}
#' \item{aic}{The Akaike information criterion (AIC) value.}
#' \item{iter}{The number of iterations performed.}
#' \item{convergence}{A logical value indicating whether the algorithm converged.}
#' \item{ctime}{The computation time (user time).}
#' \item{data}{An object containing the input data.}
#' \item{aerror}{The absolute error of the log-likelihood at the final iteration.}
#' \item{rerror}{The relative error of the log-likelihood at the final iteration.}
#' \item{options}{A list of options used for the fitting.}
#' \item{call}{The matched function call.}
#' 
#' @note
#' This method allows handling of truncated data by using \code{NA} and \code{Inf}:
#' 
#' \code{phfit.group(ph = cf1(5), counts = c(countsdata, NA), breaks = c(breakdata, +Inf))}
#' 
#' Here, \code{NA} indicates missing count data, and \code{Inf} represents an open-ended last interval \code{[last break point, infinity)}.
#' 
#' @examples
#' ## Generate a sample
#' wsample <- rweibull(n = 100, shape = 2, scale = 1)
#' wgroup <- hist(x = wsample, breaks = "fd", plot = FALSE)
#' 
#' ## Fit a general PH distribution
#' result1 <- phfit.group(ph = ph(2), counts = wgroup$counts, breaks = wgroup$breaks)
#' 
#' ## Fit a Canonical Form 1 (CF1) distribution
#' result2 <- phfit.group(ph = cf1(2), counts = wgroup$counts, breaks = wgroup$breaks)
#' 
#' ## Fit a hyper-Erlang distribution
#' result3 <- phfit.group(ph = herlang(3), counts = wgroup$counts, breaks = wgroup$breaks)
#' 
#' ## Calculate statistics
#' ph.mean(result1$model)
#' ph.var(result2$model)
#' ph.moment(5, result3$model)
#' 
#' @export

phfit.group <- function(ph, counts, breaks, intervals, instants, ...) {
  call <- match.call()
  ph <- ph$copy()
  
  options <- emoptions()
  con <- list(...)
  nmsC <- names(options)
  options[(namc <- names(con))] <- con
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  data <- data.frame.phase.group(counts=counts, breaks=breaks,
                                 intervals=intervals, instants=instants)
  
  if (options$initialize == TRUE) {
    ph$init(data, options)
  }
  tres <- system.time(result <- ph$emfit(data, options, ...))
  result <- c(result, list(model=ph, aic=-2*(result$llf - ph$df()), df=ph$df(),
                           data=data, ctime=tres[1], options=options, call=call))
  class(result) <- "phfit.result"
  result
}

#' PH fitting with a density function
#' 
#' Fits a phase-type (PH) distribution to a given density function via maximum likelihood estimation.
#' 
#' @param ph An object of R6 class for PH distributions. The estimation algorithm is selected based on this class.
#' @param f A function object representing a density function. The function should have the form \code{f(x, ...)}, where \code{x} is the first argument.
#' @param deformula An object specifying the numerical integration formula. 
#' Defaults to \code{deformula.zeroinf}, suitable for functions defined on \code{[0, +Inf)}.
#' @param weight.zero A threshold below which weights in numerical integration are regarded as zero. Defaults to \code{1e-12}.
#' @param weight.reltol The relative tolerance for numerical integration. Defaults to \code{1e-8}.
#' @param start.divisions The initial number of divisions used in numerical integration. Defaults to \code{8}.
#' @param max.iter The maximum number of iterations allowed for increasing the number of divisions. Defaults to \code{12}.
#' @param ... Additional arguments passed to the density function \code{f}, or control options for the EM algorithm.
#' 
#' @return A list of class \code{phfit.result.density} with the following components:
#' \item{model}{An object for the estimated PH distribution.}
#' \item{llf}{The maximized log-likelihood value (equivalent to the negative cross-entropy).}
#' \item{df}{The degrees of freedom of the fitted model.}
#' \item{KL}{The estimated Kullback-Leibler divergence between the target density and the fitted PH distribution.}
#' \item{iter}{The number of iterations performed.}
#' \item{convergence}{A logical value indicating whether the algorithm converged.}
#' \item{ctime}{The computation time (user time).}
#' \item{data}{An object containing the generated data points and weights.}
#' \item{aerror}{The absolute error of the log-likelihood at the final iteration.}
#' \item{rerror}{The relative error of the log-likelihood at the final iteration.}
#' \item{options}{A list of options used for the fitting.}
#' \item{call}{The matched function call.}
#' 
#' @note
#' Any proper density function can be used for the argument \code{f}. 
#' 
#' The first argument of \code{f} must be the variable of integration, while additional parameters can be set via \code{...}.
#' 
#' Truncated densities and densities on \code{[0, +Inf)} are both supported.
#' 
#' @examples
#' ## Fit PH distribution to a normal density truncated on [0, +Inf)
#' 
#' ## General PH distribution
#' result1 <- phfit.density(ph = ph(2), f = dnorm, mean = 3, sd = 1)
#' 
#' ## Canonical Form 1 (CF1) distribution
#' result2 <- phfit.density(ph = cf1(2), f = dnorm, mean = 3, sd = 1)
#' 
#' ## Hyper-Erlang distribution
#' result3 <- phfit.density(ph = herlang(3), f = dnorm, mean = 3, sd = 1)
#' 
#' ## Calculate statistics
#' ph.mean(result1$model)
#' ph.var(result2$model)
#' ph.moment(5, result3$model)
#' 
#' @export

phfit.density <- function(
    ph, f, deformula = deformula.zeroinf, weight.zero = 1.0e-12,
    weight.reltol = 1.0e-8, start.divisions = 8, max.iter = 12,
    ...) {
  call <- match.call()
  ph <- ph$copy()
  
  x <- deformula(f, ..., zero.eps = weight.zero,
                 rel.tol = weight.reltol,
                 start.divisions = start.divisions, max.iter = max.iter)
  ll <- sum(x$w * log(f(x$x, ...)))
  data <- data.frame.phase.time(x=x$x, weights=x$w)
  
  options <- emoptions()
  con <- list(...)
  nmsC <- names(options)
  options[(namc <- names(con))] <- con
  # if (length(noNms <- namc[!namc %in% nmsC])) 
  #   warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  if (options$initialize == TRUE) {
    ph$init(data, options)
  }
  tres <- system.time(result <- ph$emfit(data, options, ...))
  result <- c(result, list(model=ph, KL=x$h * (ll-result$llf), df=ph$df(),
                           data=data, ctime=tres[1], options=options, call=call))
  class(result) <- "phfit.result.density"
  result
}

#' PH fitting with left-truncated and right-censored data
#' 
#' Fits a phase-type (PH) distribution to left-truncated and right-censored survival data via maximum likelihood estimation.
#' 
#' @param ph An object of R6 class for PH distributions. The estimation algorithm is selected based on this class.
#' @param x A numeric vector of observed times (either event times or censoring times).
#' @param delta A numeric vector of event indicators.
#' If \code{delta = 1}, the observation is an \strong{event time} (i.e., the event of interest occurred).
#' If \code{delta = 0}, the observation is \strong{right-censored}.
#' @param tau A numeric vector of left-truncation times.
#' If omitted, all observations are assumed to have no left truncation (\code{NA}).
#' @param ... Additional options for the fitting algorithm.
#' 
#' @return A list of class \code{phfit.result} with the following components:
#' \item{model}{An object for the estimated PH distribution.}
#' \item{llf}{The maximized log-likelihood value.}
#' \item{df}{The degrees of freedom of the fitted model.}
#' \item{aic}{The Akaike information criterion (AIC) value.}
#' \item{iter}{The number of iterations performed.}
#' \item{convergence}{A logical value indicating whether the algorithm converged.}
#' \item{ctime}{The computation time (user time).}
#' \item{data}{An object containing the input data.}
#' \item{aerror}{The absolute error of the log-likelihood at the final iteration.}
#' \item{rerror}{The relative error of the log-likelihood at the final iteration.}
#' \item{options}{A list of options used for the fitting.}
#' \item{call}{The matched function call.}
#'
#' @note
#' \code{herlang} cannot be directly used with \code{phfit.surv} yet.
#' Use general PH (\code{ph}) or Canonical Form 1 (\code{cf1}) models.
#'
#' @examples
#' ## Generate a survival dataset
#' set.seed(123)
#' true_times <- rexp(100, rate = 0.5)  # true event times
#' censoring_times <- rexp(100, rate = 0.3)  # censoring times
#' 
#' x <- pmin(true_times, censoring_times)
#' delta <- as.integer(true_times <= censoring_times)  # 1 = event, 0 = censored
#' tau <- rep(NA, length(x))  # no left truncation
#' 
#' ## Fit a general PH distribution
#' result1 <- phfit.surv(ph = ph(2), x = x, delta = delta, tau = tau)
#' 
#' ## Fit a Canonical Form 1 (CF1) distribution
#' result2 <- phfit.surv(ph = cf1(2), x = x, delta = delta, tau = tau)
#' 
#' ## Calculate statistics
#' ph.mean(result1$model)
#' ph.var(result2$model)
#' ph.moment(5, result1$model)
#' ph.moment(5, result2$model)
#' 
#' @export

phfit.surv <- function(ph, x, delta, tau, ...) {
  call <- match.call()
  ph <- ph$copy()
  
  options <- emoptions()
  con <- list(...)
  nmsC <- names(options)
  options[(namc <- names(con))] <- con
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  data <- data.frame.phase.surv(x=x, delta=delta, tau=tau)
  if (options$initialize == TRUE) {
    ph$init(data, options)
  }
  tres <- system.time(result <- ph$emfit(data, options, ...))
  result <- c(result, list(model=ph, aic=-2*(result$llf - ph$df()), df=ph$df(),
                           data=data, ctime=tres[1], options=options, call=call))
  class(result) <- "phfit.result"
  result
}

#' PH fitting with interval data
#' 
#' Fits a phase-type (PH) distribution to interval data via maximum likelihood estimation.
#' 
#' @param ph An object of R6 class for PH distributions. The estimation algorithm is selected based on this class.
#' @param data A list consisting of point data and interval data.
#' @param weights A numeric vector of weights.
#' @param ... Additional options for the fitting algorithm.
#' 
#' @return A list of class \code{phfit.result} with the following components:
#' \item{model}{An object for the estimated PH distribution.}
#' \item{llf}{The maximized log-likelihood value.}
#' \item{df}{The degrees of freedom of the fitted model.}
#' \item{aic}{The Akaike information criterion (AIC) value.}
#' \item{iter}{The number of iterations performed.}
#' \item{convergence}{A logical value indicating whether the algorithm converged.}
#' \item{ctime}{The computation time (user time).}
#' \item{data}{An object containing the input data.}
#' \item{aerror}{The absolute error of the log-likelihood at the final iteration.}
#' \item{rerror}{The relative error of the log-likelihood at the final iteration.}
#' \item{options}{A list of options used for the fitting.}
#' \item{call}{The matched function call.}
#'
#' @note
#' \code{herlang} cannot be directly used with \code{phfit.interval} yet.
#' Use general PH (\code{ph}) or Canonical Form 1 (\code{cf1}) models.
#' 
#' @export

phfit.interval <- function(ph, data, weights, ...) {
  call <- match.call()
  ph <- ph$copy()
  
  options <- emoptions()
  con <- list(...)
  nmsC <- names(options)
  options[(namc <- names(con))] <- con
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  data <- data.frame.phase.interval(data=data, weights=weights)
  if (options$initialize == TRUE) {
    ph$init(data, options)
  }
  tres <- system.time(result <- ph$emfit(data, options, ...))
  result <- c(result, list(model=ph, aic=-2*(result$llf - ph$df()), df=ph$df(),
                           data=data, ctime=tres[1], options=options, call=call))
  class(result) <- "phfit.result"
  result
}

#' @aliases phfit.point phfit.group phfit.surv phfit.interval
#' @export

print.phfit.result <- function (x, ...) {
  cat("\n")
  cat(sprintf("Maximum LLF: %f\n", x$llf))
  cat(sprintf("DF: %d\n", x$df))
  cat(sprintf("AIC: %f\n", x$aic))
  cat(sprintf("Iteration:  %d / %d\n", x$iter, x$options$maxiter))
  cat(sprintf("Computation time (user): %f\n", x$ctime))
  cat(sprintf("Convergence: %s\n", x$convergence))
  cat(sprintf("Error (abs): %e (tolerance %e)\n", x$aerror, x$options$abstol))
  cat(sprintf("Error (rel): %e (tolerance %e)\n", x$rerror, x$options$reltol))
  cat("\n")
  x$model$print(...)
  cat("\n\n")
  invisible(x)
}

#' @aliases phfit.density
#' @export

print.phfit.result.density <- function (x, ...) {
  cat("\n")
  cat(sprintf("Maximum LLF: %f\n", x$llf))
  cat(sprintf("DF: %d\n", x$df))
  cat(sprintf("KL: %f\n", x$KL))
  cat(sprintf("Iteration:  %d / %d\n", x$iter, x$options$maxiter))
  cat(sprintf("Computation time (user): %f\n", x$ctime))
  cat(sprintf("Convergence: %s\n", x$convergence))
  cat(sprintf("Error (abs): %e (tolerance %e)\n", x$aerror, x$options$abstol))
  cat(sprintf("Error (rel): %e (tolerance %e)\n", x$rerror, x$options$reltol))
  cat("\n")
  x$model$print(...)
  cat("\n\n")
  invisible(x)
}
