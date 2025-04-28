#' MAP fitting with point data
#' 
#' Fits a Markovian Arrival Process (MAP) to point data via maximum likelihood estimation.
#' 
#' @param map An object of R6 class representing a MAP. The estimation algorithm is selected based on this class.
#' @param x A numeric vector of observed point data (event times).
#' @param intervals A numeric vector of interval lengths. Optional.
#' @param ... Additional options for the fitting algorithm.
#' 
#' @return A list of class \code{mapfit.result} containing the following components:
#' \item{model}{An object representing the estimated MAP.}
#' \item{llf}{The maximized log-likelihood value.}
#' \item{df}{The degrees of freedom of the fitted model.}
#' \item{aic}{The Akaike information criterion (AIC) value.}
#' \item{iter}{The number of iterations performed.}
#' \item{convergence}{A logical value indicating whether the algorithm converged.}
#' \item{ctime}{The computation time (user time).}
#' \item{data}{An object containing the input data set.}
#' \item{aerror}{The absolute error of the log-likelihood at the final iteration.}
#' \item{rerror}{The relative error of the log-likelihood at the final iteration.}
#' \item{options}{A list of options used for the fitting.}
#' \item{call}{The matched function call.}
#'
#' @examples
#' ## Load trace data
#' data(BCpAug89)
#' BCpAug89s <- head(BCpAug89, 50)
#' 
#' ## Fit a general MAP
#' result1 <- mapfit.point(map = map(2), x = cumsum(BCpAug89s))
#' 
#' ## Fit a Markov Modulated Poisson Process (MMPP)
#' result2 <- mapfit.point(map = mmpp(2), x = cumsum(BCpAug89s))
#' 
#' ## Fit an Erlang Renewal Hidden Markov Model (ER-HMM)
#' result3 <- mapfit.point(map = erhmm(3), x = cumsum(BCpAug89s))
#' 
#' ## Marginal moments
#' map.mmoment(k = 3, map = result1$model)
#' map.mmoment(k = 3, map = result2$model)
#' map.mmoment(k = 3, map = result3$model)
#' 
#' ## Joint moments
#' map.jmoment(lag = 1, map = result1$model)
#' map.jmoment(lag = 1, map = result2$model)
#' map.jmoment(lag = 1, map = result3$model)
#' 
#' ## Lag-k autocorrelation
#' map.acf(map = result1$model)
#' map.acf(map = result2$model)
#' map.acf(map = result3$model)
#' 
#' @export

mapfit.point <- function(map, x, intervals, ...) {
  call <- match.call()
  map <- map$copy()

  options <- emoptions()
  con <- list(...)
  nmsC <- names(options)
  options[(namc <- names(con))] <- con
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  data <- data.frame.map.time(time=x, intervals=intervals)
  if (options$initialize == TRUE) {
    map$init(data, options)
  }
  tres <- system.time(result <- map$emfit(data, options, ...))
  result <- c(result, list(model=map, aic=-2*(result$llf - map$df()), df=map$df(),
                           data=data, ctime=tres[1], options=options, call=call))
  class(result) <- "mapfit.result"
  result
}

#' MAP fitting with grouped data
#' 
#' Fits a Markovian Arrival Process (MAP) to grouped data via maximum likelihood estimation.
#' 
#' @param map An object of R6 class representing a MAP. The estimation algorithm is selected based on this class.
#' @param counts A numeric vector of counts within each interval.
#' @param breaks A numeric vector specifying the boundaries of the intervals. Equivalent to \code{c(0, cumsum(intervals))}.
#' If omitted, defaults to \code{0:length(counts)}.
#' @param intervals A numeric vector of interval lengths. Equivalent to \code{diff(breaks)}.
#' If omitted, defaults to \code{rep(1, length(counts))}.
#' @param instants A numeric vector indicating whether an instantaneous event is observed at the end of each interval.
#' If \code{instants[i] = 1}, a sample is observed at the right boundary of the \code{i}-th interval.
#' Defaults to \code{rep(0L, length(counts))} (no instantaneous events).
#' @param ... Additional options for the fitting algorithm.
#' 
#' @return A list of class \code{mapfit.result} containing the following components:
#' \item{model}{An object representing the estimated MAP.}
#' \item{llf}{The maximized log-likelihood value.}
#' \item{df}{The degrees of freedom of the fitted model.}
#' \item{aic}{The Akaike information criterion (AIC) value.}
#' \item{iter}{The number of iterations performed.}
#' \item{convergence}{A logical value indicating whether the algorithm converged.}
#' \item{ctime}{The computation time (user time).}
#' \item{data}{An object containing the input data set.}
#' \item{aerror}{The absolute error of the log-likelihood at the final iteration.}
#' \item{rerror}{The relative error of the log-likelihood at the final iteration.}
#' \item{options}{A list of options used for the fitting.}
#' \item{call}{The matched function call.}
#'
#' @examples
#' ## Load trace data
#' data(BCpAug89)
#' BCpAug89s <- head(BCpAug89, 50)
#' 
#' ## Create grouped data
#' BCpAug89.group <- hist(cumsum(BCpAug89s),
#'                        breaks = seq(0, 0.15, 0.005),
#'                        plot = FALSE)
#' 
#' ## Fit a general MAP
#' result1 <- mapfit.group(map = map(2),
#'                         counts = BCpAug89.group$counts,
#'                         breaks = BCpAug89.group$breaks)
#' 
#' ## Fit a Markov Modulated Poisson Process (MMPP)
#' result2 <- mapfit.group(map = mmpp(2),
#'                         counts = BCpAug89.group$counts,
#'                         breaks = BCpAug89.group$breaks)
#' 
#' ## Fit an approximate MMPP (G-MMPP)
#' result3 <- mapfit.group(map = gmmpp(2),
#'                         counts = BCpAug89.group$counts,
#'                         breaks = BCpAug89.group$breaks)
#' 
#' ## Marginal moments
#' map.mmoment(k = 3, map = result1$model)
#' map.mmoment(k = 3, map = result2$model)
#' map.mmoment(k = 3, map = result3$model)
#' 
#' ## Joint moments
#' map.jmoment(lag = 1, map = result1$model)
#' map.jmoment(lag = 1, map = result2$model)
#' map.jmoment(lag = 1, map = result3$model)
#' 
#' ## Lag-k autocorrelation
#' map.acf(map = result1$model)
#' map.acf(map = result2$model)
#' map.acf(map = result3$model)
#' 
#' @export

mapfit.group <- function(map, counts, breaks, intervals, instants, ...) {
  call <- match.call()
  map <- map$copy()
  
  options <- emoptions()
  con <- list(...)
  nmsC <- names(options)
  options[(namc <- names(con))] <- con
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  data <- data.frame.map.group(counts=counts, breaks=breaks,
                               intervals=intervals, instants=instants)

  if (options$initialize == TRUE) {
    map$init(data, options)
  }
  tres <- system.time(result <- map$emfit(data, options, ...))
  result <- c(result, list(model=map, aic=-2*(result$llf - map$df()), df=map$df(),
                           data=data, ctime=tres[1], options=options, call=call))
  class(result) <- "mapfit.result"
  result
}

#' @aliases mapfit.point mapfit.group
#' @export

print.mapfit.result <- function (x, ...) {
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
