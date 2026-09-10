#' EM Options
#' 
#' A list of options for EM
#' 
#' @return A list of options with default values
#' 
#' @export

emoptions <- function() {
  list(
    maxiter = 2000,
    abstol = +Inf,
    reltol = sqrt(.Machine$double.eps),
    em.verbose = FALSE,
    steps = 1,
    uniform.factor = 1.01,
    poisson.eps = sqrt(.Machine$double.eps),
    initialize = TRUE,
    cf1.diff.init = c(1, 4, 16, 64, 256, 1024),
    cf1.scale.init = c(0.5, 1.0, 2.0),
    cf1.verbose = TRUE,
    cf1.maxiter = 5,
    lbound = 1,
    ubound = NA,
    shape.verbose = TRUE,
    shape.method = "all",
    map.stationary = TRUE,
    inte.divide=30,
    inte.eps=1.0e-8,
    annealing = FALSE,
    annealing.temperature = seq(0.9, 1, length.out=10),
    annealing.iter = NULL
  )
}

# Coerce to a dense general Matrix (dgeMatrix).
#
# The EM kernels for MAP and ER-HMM read D0, D1 and P through
# S4matrix<DenseMatrixT>, so these slots have to be dgeMatrix. Coercing
# directly with as(x, "dgeMatrix") is deprecated in Matrix, and the coercions
# to the virtual classes "unpackedMatrix" and "CsparseMatrix" detect symmetry
# and triangularity: a symmetric argument would become a dsyMatrix and an
# argument with a zero triangle a dtrMatrix, neither of which has the dense
# general layout the C++ code assumes. This three-step coercion is the
# replacement documented in help("Matrix-deprecated") and always gives a
# dgeMatrix.
as.dge <- function(x) {
  as(as(as(x, "dMatrix"), "generalMatrix"), "unpackedMatrix")
}

#' Markov stationary
#' 
#' Compute the stationary vector with GTH
#' 
#' @param Q DTMC/CTMC kernel
#' @return The stationary vector of DTMC/CTMC
#' 
#' @export

ctmc.st <- function(Q) {
  A <- as.matrix(Q)
  if (dim(Q)[1] != dim(Q)[2]) {
    stop("Q should be a square matrix")
  }
  x <- numeric(dim(Q)[1])
  markov_gth_dense(Q, x)
}