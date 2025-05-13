#' @export
Remstep_cf1_interval <- function(data, weights, params) {
    n <- length(params$alpha)
    alpha <- params$alpha
    rate <- params$rate
    options <- list(uniform.factor=1.01, poisson.eps=1.0e-8)
    dat <- data.frame.phase.interval(data=data, weights=weights)
    if (n >= 2) {
        i <- c(1:n, 1:(n-1))
        j <- c(1:n, 2:n)
        x <- c(-rate, rate[1:(n-1)])
        Q <- sparseMatrix(i=i, j=j, x=x, repr = "T")
    } else {
        Q <- matrix(-rate[1],1,1)
    }
    matclass <- "CsparseMatrix"
    emstep_cf1_interval(alpha, rate, dat, options,
        as(Q, matclass),
        as(Q, matclass),
        as(Q, matclass))
}
