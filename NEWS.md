# mapfit 1.0.1

- Fixed the handling of left-truncated grouped data. `breaks` not starting at 0
  add an interval whose count is unknown, and it was marked with `NA` after the
  conversion of missing counts to the internal sentinel. The estimation step
  then treated the interval as if it had been observed to contain no samples.
- Missing counts (`NA`) and left-truncated data are now rejected in
  `mapfit.group()`. The estimation algorithm for MAP uses the number of arrivals
  in an interval as an array index, so an unknown count read past the end of the
  array instead of being marginalised out. Missing counts remain supported in
  `phfit.group()`, which is what the documentation has always described.
- Fixed the exit rate of a CF1 after estimation. `phfit.point()` and friends
  updated the initial probabilities and the rates but left the exit rate at its
  initial value, so `dphase()` on an estimated CF1 returned a multiple of the
  density; the values did not integrate to 1. The degrees of freedom reported
  for a CF1 are recomputed from the estimated parameters as well.
- `cf1()` no longer sorts the vectors passed as `alpha` and `rate` in the
  caller's environment.
- Passing data to a model that cannot handle it now raises an error. In
  particular `mapfit.point()` with a `gmmpp` model, which is an algorithm for
  grouped data, silently returned a result object with no log-likelihood.

# mapfit 1.0.0

- Refactor with Rcpp
- Use traits for the implementation of vec-mat operations
- emfit is written by Rcpp
- Change the options for emfit

# mapfit 0.9.9

- Added a `NEWS.md` file to track changes to the package.
- Change the license
- Use Roxygen2



