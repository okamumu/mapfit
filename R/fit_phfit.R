

phfit.point <- function(ph, x, weights, method = c("all", "increment"),
  lbound = 1, ubound = NULL, control = list(), verbose = list(), ...) {
  data <- phfit.time.data.frame(time=x, weights=weights)
  switch(class(ph),
    "ph"=phfit.gen(ph=ph, data=data, control=control, verbose=verbose, ...),
    "cf1"=phfit.cf1(ph=ph, data=data, control=control, verbose=verbose, ...),
    "herlang"={
      phsize <- sum(ph@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
      phfit.herlang(phsize=phsize, data=data, method=method, lbound=lbound, ubound=ubound,
        control=control, verbose=verbose, ...)
    })
}

phfit.group <- function(ph, counts, breaks, intervals, instant,
 method = c("all", "increment"), lbound = 1, ubound = NULL, control = list(), verbose = list(), ...) {
  data <- phfit.group.data.frame(counts=counts, breaks=breaks, difftime=intervals, instant=instant)
  switch(class(ph),
    "ph"=phfit.gen(ph=ph, data=data, control=control, verbose=verbose, ...),
    "cf1"=phfit.cf1(ph=ph, data=data, control=control, verbose=verbose, ...),
    "herlang"={
      phsize <- sum(ph@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
    	phfit.herlang(phsize=phsize, data=data, method=method, lbound=lbound, ubound=ubound,
    		control=control, verbose=verbose, ...)
    })
}

phfit.density <- function(ph, f,
    method = c("all", "increment"), lbound = 1, ubound = NULL, 
    deformula = deformula.zeroinf, weight.zero = 1.0e-12,
    weight.reltol = 1.0e-8, start.divisions = 8, max.iter = 12,
    control = list(), verbose = list(), ...) {
  x <- deformula(f, ..., zero.eps = weight.zero, rel.tol = weight.reltol, start.divisions = start.divisions,
    max.iter = max.iter)
  ll <- x$h * sum(x$w * log(f(x$x, ...)))
  data <- phfit.time.data.frame(time=x$x, weights=x$w)
  res <- switch(class(ph),
    "ph"=phfit.gen(ph=ph, data=data, control=control, verbose=verbose, ...),
    "cf1"=phfit.cf1(ph=ph, data=data, control=control, verbose=verbose, ...),
    "herlang"={
      phsize <- sum(ph@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
      phfit.herlang(phsize=phsize, data=data, method=method, lbound=lbound, ubound=ubound,
        control=control, verbose=verbose, ...)
    })
  c(res, list(KL=ll - res$llf * x$h))
}
