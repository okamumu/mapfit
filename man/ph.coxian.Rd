% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_ph.R
\name{ph.coxian}
\alias{ph.coxian}
\title{Create a Coxian PH distribution}
\usage{
ph.coxian(size)
}
\arguments{
\item{size}{An integer for the number of phases}
}
\value{
An instance of coxian PH distribution
}
\description{
Create an instance of coxian PH distribution.
}
\note{
Coxian PH distribution is the PH distribution whose infinitesimal
generator is given by a upper bi-diagonal matrix. This is also called
canonical form 3.
}
\examples{
## create a Coxian PH with 5 phases
(param1 <- ph.coxian(5))

}
