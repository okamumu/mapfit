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
    herlang.lbound = 1,
    herlang.ubound = NA,
    herlang.verbose = TRUE,
    annealing = FALSE,
    annealing.temperature = seq(0.9, 1, length.out=10),
    annealing.iter = NULL
    )
}
