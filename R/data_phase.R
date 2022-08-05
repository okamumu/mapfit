#' Create data for phase with weighted sample
#' 
#' Provide a data.frame with weighted samples.
#' 
#' @param x A vector of point (quantiles)
#' @param weights A vector of weights
#' @return A dataframe
#' 
#' @note 
#' The point time is sorted and their differences are stored as the column of `time`
#' 
#' @examples
#' x <- runif(10)
#' w <- runif(10)
#' 
#' dat <- data.frame.phase.time(x=x, weights=w)
#' print(dat)
#' mean(dat)
#' 
#' @export

data.frame.phase.time <- function(x, weights) {
  if (missing(weights)) {
    weights <- array(1, length(x))
  }
  tt <- diff(c(0, sort(x)))

  # check
  if (! (length(tt) == length(weights))) {
    stop(sprintf("The length of time and weights should be same. time=%d, weights=%d",
                 length(time), length(weights)))
  }

  data <- list(time=tt, weights=weights[order(x)],
               maxtime=max(tt))
  class(data) <- "phase.time"
  data
}

#' @aliases data.frame.phase.time
#' @export
print.phase.time <- function(x, ...) {
  print(data.frame(x=x$time, weights=x$weights))
}

#' @aliases data.frame.phase.time
#' @export
mean.phase.time <- function(x, ...) {
  sum(cumsum(x$time) * x$weights) / sum(x$weights)
}

#' Create group data for phase
#' 
#' Provide the data.frame for group data.
#' 
#' @param counts A vector of the number of samples
#' @param breaks A vector of break points
#' @param intervals A vector of differences of time
#' @param instants A vector meaning whether a sample is observed at the end of break.
#' @return A dataframe
#' @examples
#' dat <- data.frame.phase.group(counts=c(1,2,1,1,0,0,1,4))
#' print(dat)
#' mean(dat)
#' 
#' @export

data.frame.phase.group <- function(counts, breaks, intervals, instants) {
  # replace na to -1
  counts[is.na(counts)] <- -1
  
  if (missing(breaks)) {
    if (missing(intervals)) {
      breaks <- 0:length(counts)
    } else {
      breaks <- c(0,cumsum(intervals))
    }
  }
  # check for glast
  if (is.infinite(breaks[length(breaks)])) {
    glast <- counts[length(counts)]
    counts <- counts[-length(counts)]
    breaks <- breaks[-length(breaks)]
  } else {
    glast <- 0
  }
  # check instants
  if (missing(instants)) {
    instants <- array(0, length(counts))
  }
  # check for left-truncation
  if (breaks[1] != 0) {
    breaks <- c(0, breaks)
    counts <- c(NA, counts)
    instants <- c(0, instants)
  }
  dt <- diff(breaks)
  
  # check
  if (! (length(dt) == length(counts) && length(dt) == length(instants))) {
    stop(sprintf("The length of time, counts, indicators should be same. intervals=%d, counts=%d, instants=%d",
         length(dt), length(counts), length(instants)))
  }
  
  data <- list(
    counts = counts,
    intervals = dt,
    instants = instants,
    maxinterval = max(dt),
    maxcount = max(counts),
    lastcount = glast)
  class(data) <- "phase.group"
  data
}

#' @aliases data.frame.phase.group
#' @export
print.phase.group <- function(x, ...) {
  print(data.frame(counts=c(x$counts,x$last), intervals=c(x$intervals,+Inf), instants=c(x$instants,NA)))
}

#' @aliases data.frame.phase.group
#' @export
mean.phase.group <- function(x, ...) {
  i <- is.finite(x$intervals) & is.finite(x$counts) & is.finite(x$instants)
  time <- cumsum(x$intervals[i])
  ind <- x$instants[i] == 1
  counts <- x$counts[i]
  m <- sum(time * counts) + sum(time[ind])
  m / (sum(counts) + sum(ind))
}

#' Create left-truncated and right-censored data for phase
#' 
#' Provide the data.frame for left-truncated and right-censored data.
#' 
#' @param x A vector of time points
#' @param delta A vector of indicators whether x is censoring time or not. If delta=1, the corresponding x is the censoring time
#' If delta=0, the corresponding x is the event time.
#' @param tau A vector of left-truncation time points. If tau is missing, all the left-truncation times are NA (no truncation).
#' @return A dataframe
#' @examples
#' dat <- data.frame.phase.surv(x=c(1.0, 4.5, 1.2), delta=c(1, 0, 0), tau=c(NA, 2.0, NA))
#' print(dat)
#' mean(dat)
#' 
#' @export

data.frame.phase.surv <- function(x, delta, tau) {
  if (missing(delta)) {
    delta <- rep(0, length(x))
  }
  if (missing(tau)) {
    tau <- rep(NA, length(x))
  }

  # length check
  if (! (length(x) == length(delta) && length(x) == length(tau))) {
    stop(sprintf("The length of time, counts, indicators should be same. x=%d, delta=%d, tau=%d",
                 length(x), length(delta), length(tau)))
  }
  # value check
  s <- is.finite(tau)
  if (all(x[s] >= tau[s]) != TRUE) {
    stop(sprintf("The left-truncation time should be less than or equal to x"))
  }
  if (all(delta == 1 | delta == 0) != TRUE) {
    stop(sprintf("The delta should be 0 or 1"))
  }
  
  # replace na to 0.0
  tau[is.na(tau)] <- 0.0
  d <- sapply(delta, function(v) ifelse(v==1, 0, 1))

  t <- c(x, tau)
  nu <- c(d, rep(3, length(x)))
  
  s <- order(t)
  dt <- diff(c(0, t[s]))
  nu <- nu[s]
  
  data <- list(
    intervals = dt,
    nu = nu,
    maxinterval = max(dt))
  class(data) <- "phase.surv"
  data
}

#' @aliases data.frame.phase.surv
#' @export
print.phase.surv <- function(x, ...) {
  print(data.frame(intervals=x$intervals, nu=x$nu))
}

#' @aliases data.frame.phase.surv
#' @export
mean.phase.surv <- function(x, ...) {
  s <- (x$nu == 0) | (x$nu == 1)
  t <- cumsum(x$intervals)[s]
  m <- sum(t)
  m / length(t)
}
