% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phfit.R
\name{phfit.density}
\alias{phfit.density}
\title{PH fitting with density function}
\usage{
phfit.density(
  ph,
  f,
  deformula = deformula.zeroinf,
  weight.zero = 1e-12,
  weight.reltol = 1e-08,
  start.divisions = 8,
  max.iter = 12,
  ...
)
}
\arguments{
\item{ph}{An object of R6 class. The estimation algorithm is selected depending on this class.}

\item{f}{A function object for a density function.}

\item{deformula}{An object for formulas of numerical integration.
It is not necessary to change it when the density function is defined on
the positive domain [0,infinity).}

\item{weight.zero}{A absolute value which is regarded as zero in numerical integration.}

\item{weight.reltol}{A value for precision of numerical integration.}

\item{start.divisions}{A value for starting value of divisions in deformula.}

\item{max.iter}{A value for the maximum number of iterations to increase divisions in deformula.}

\item{...}{Options for EM steps, which is also used to send the arguments to density function.}
}
\value{
Returns a list with components, which is an object of S3 class \code{phfit.result};
\item{model}{an object for estimated PH class.}
\item{llf}{a value of the maximum log-likelihood (a negative value of the cross entropy).}
\item{df}{a value of degrees of freedom of the model.}
\item{KL}{a value of Kullback-Leibler divergence.}
\item{iter}{the number of iterations.}
\item{convergence}{a logical value for the convergence of estimation algorithm.}
\item{ctime}{computation time (user time).}
\item{data}{an object for data class}
\item{aerror}{a value of absolute error for llf at the last step of algorithm.}
\item{rerror}{a value of relative error for llf at the last step of algorithm.}
\item{options}{a list of options.}
\item{call}{the matched call.}
}
\description{
Estimates PH parameters from density function.
}
\note{
Any of density function can be applied to the argument \code{f}, where
\code{f} should be defined \code{f <- function(x, ...)}.
The first argument of \code{f} should be an integral parameter.
The other parameters are set in the argument \code{...} of \code{phfit.density}.
The truncated density function can also be used directly.
}
\examples{
####################
##### truncated density
####################

## PH fitting for general PH
(result1 <- phfit.density(ph=ph(2), f=dnorm, mean=3, sd=1))

## PH fitting for CF1
(result2 <- phfit.density(ph=cf1(2), f=dnorm, mean=3, sd=1))

## PH fitting for hyper Erlang
(result3 <- phfit.density(ph=herlang(3), f=dnorm, mean=3, sd=1))

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
