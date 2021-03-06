\name{mapfit.point}
\alias{mapfit.point}
\title{MAP fitting with time point data}

\description{
  estimates MAP parameters from time point data.
}
\usage{
mapfit.point (map, x, intervals, stationary = TRUE,
  method = c("all", "increment"), lbound = 1, ubound = NULL,
  control = list(), verbose = list(), ...)
}
\arguments{
  \item{map}{an object of S4 class for MAP. The estimation algorithm is selected depending on thie class.}
  \item{x}{a vector for time sequence of arrivals. This is equivalent to \code{cumsum(intervals)}. Either time or difftime should be given.}
  \item{intervals}{a vector for the data for intrarrival time. This is equivalent to \code{diff(c(0,x)}). Either time or difftime should be given.}
  \item{stationary}{a logical value that determine whether initial probability is given by a stationary vector of underlying Markov process or not.}
  \item{method}{the name of estimation method for ER-HMM (\code{\linkS4class{erhmm}}).}
  \item{lbound}{a value for lower limit for the number of states in ER-HMM (\code{\linkS4class{erhmm}}).}
  \item{ubound}{a value for upper limit for the number of states in ER-HMM (\code{\linkS4class{erhmm}}).}
  \item{control}{a list of parameters for controlling the fitting process.}
  \item{verbose}{a list of parameters for displaying the fitting process.}
  \item{\dots}{further arguments for methods.}
}
\value{
  returns a list with components, which is an object of S3 class \code{mapfit.result};
  \item{model}{an object for estimated MAP class (\code{\linkS4class{map}}, \code{\linkS4class{erhmm}}).}
  \item{llf}{a value of the maximum log-likelihood.}
  \item{df}{a value of degrees of freedom of the model.}
  \item{aic}{a value of Akaike information criterion.}
  \item{iter}{the number of iterations.}
  \item{convergence}{a logical value for the convergence of estimation algorithm.}
  \item{ctime}{computation time (user time).}
  \item{stationary}{a logical value for the argument \code{stationary}.}
  \item{data}{an object for MAP data class}
  \item{aerror}{a value of absolute error for llf at the last step of algorithm.}
  \item{rerror}{a value of relative error for llf at the last step of algorithm.}
  \item{control}{a list of the argument of \code{control}.}
  \item{verbose}{a list of the argument of \code{verbose}.}
  \item{call}{the matched call.}
}
\seealso{
  \code{\link{mapfit.group}}, \code{\linkS4class{map}} and \code{\linkS4class{erhmm}}
}
\examples{
## load trace data
data(BCpAug89)
BCpAug89s <- head(BCpAug89, 50)

## MAP fitting for general MAP
(result1 <- mapfit.point(map=map(2), x=cumsum(BCpAug89s)))

## MAP fitting for MMPP
(result2 <- mapfit.point(map=mmpp(2), x=cumsum(BCpAug89s)))

## MAP fitting for ER-HMM
(result3 <- mapfit.point(map=erhmm(3), x=cumsum(BCpAug89s)))

## marginal moments for estimated MAP
map.mmoment(k=3, map=result1$model)
map.mmoment(k=3, map=result2$model)
map.mmoment(k=3, map=as(result3$model, "map"))

## joint moments for estimated MAP
map.jmoment(lag=1, map=result1$model)
map.jmoment(lag=1, map=result2$model)
map.jmoment(lag=1, map=as(result3$model, "map"))

## lag-k correlation
map.acf(map=result1$model)
map.acf(map=result2$model)
map.acf(map=as(result3$model, "map"))

}
