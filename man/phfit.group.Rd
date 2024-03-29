% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phfit.R
\name{phfit.group}
\alias{phfit.group}
\title{PH fitting with grouped data}
\usage{
phfit.group(ph, counts, breaks, intervals, instants, ...)
}
\arguments{
\item{ph}{An object of R6 class. The estimation algorithm is selected depending on this class.}

\item{counts}{A vector of the number of points in intervals.}

\item{breaks}{A vector for a sequence of points of boundaries of intervals.
This is equivalent to \code{c(0,cumsum(intervals))}.
If this is missing, it is assigned to \code{0:length(counts)}.}

\item{intervals}{A vector of time lengths for intervals.
This is equivalent to \code{diff(breaks)}).
If this is missing, it is assigned to \code{rep(1,length(counts))}.}

\item{instants}{A vector of integers to indicate whether sample is drawn at
the last of interval. If instant is 1, a sample is drawn at the last of interval.
If instant is 0, no sample is drawn at the last of interval.
By using instant, point data can be expressed by grouped data.
If instant is missing, it is given by \code{rep(0L,length(counts))}, i.e.,
there are no samples at the last of interval.}

\item{...}{Further options for EM steps.}
}
\value{
Returns a list with components, which is an object of S3 class \code{phfit.result};
\item{model}{an object for estimated PH class.}
\item{llf}{a value of the maximum log-likelihood.}
\item{df}{a value of degrees of freedom of the model.}
\item{aic}{a value of Akaike information criterion.}
\item{iter}{the number of iterations.}
\item{convergence}{a logical value for the convergence of estimation algorithm.}
\item{ctime}{computation time (user time).}
\item{data}{an object for data class}
\item{aerror}{a value of absolute error for llf at the last step of algorithm.}
\item{rerror}{a value of relative error for llf at the last step of algorithm.}
\item{options}{a list of options used in the fitting.}
\item{call}{the matched call.}
}
\description{
Estimates PH parameters from grouped data.
}
\note{
In this method, we can handle truncated data using \code{NA} and \code{Inf};
\code{phfit.group(ph=cf1(5), counts=c(countsdata, NA), breaks=c(breakdata, +Inf))}
\code{NA} means missing of count data at the corresponding interval, and \code{Inf} is allowed to put
the last of breaks or intervals which represents a special interval [the last break point,infinity).
}
\examples{
## make sample
wsample <- rweibull(n=100, shape=2, scale=1)
wgroup <- hist(x=wsample, breaks="fd", plot=FALSE)

## PH fitting for general PH
(result1 <- phfit.group(ph=ph(2), counts=wgroup$counts, breaks=wgroup$breaks))

## PH fitting for CF1
(result2 <- phfit.group(ph=cf1(2), counts=wgroup$counts, breaks=wgroup$breaks))

## PH fitting for hyper Erlang
(result3 <- phfit.group(ph=herlang(3), counts=wgroup$counts, breaks=wgroup$breaks))

## mean
ph.mean(result1$model)
ph.mean(result2$model)
ph.mean(result3$model)

## variance
ph.var(result1$model)
ph.var(result2$model)
ph.var(result3$model)

## up to 5 moments 
ph.moment(5, result1$model)
ph.moment(5, result2$model)
ph.moment(5, result3$model)

}
