## Resubmission of an archived package

mapfit was archived on 2026-02-15 because it imports deformula, which was archived on the
same day. deformula 0.1.3 has been submitted with the C++11 specification issue corrected;
this submission depends on it being back on CRAN.

This version also fixes the following bugs, none of which was reported by the CRAN checks:

* Grouped data whose breaks do not start at 0 were mishandled. The leading interval, whose
  number of samples is unknown, was marked with `NA` after the missing counts had already
  been converted to the internal sentinel, so the estimation step treated it as an interval
  observed to contain no samples.
* `mapfit.group()` now rejects missing counts (`NA`) and left-truncated data. The estimation
  algorithm for MAP uses the number of arrivals in an interval as an array index, so an
  unknown count read past the end of the array. Missing counts remain supported in
  `phfit.group()`, which is what the documentation describes.
* The exit rate of a CF1 was not updated by the estimation, so `dphase()` on an estimated
  CF1 returned a multiple of the density instead of the density itself.
* Passing data to a model that cannot handle it now raises an error instead of returning a
  result object with no log-likelihood.

The deprecated coercion `as(<matrix>, "dgeMatrix")` has been replaced with the three-step
coercion documented in `help("Matrix-deprecated")`.

## Test environments

* local Docker, Ubuntu 24.04, R 4.5.1 (release), R CMD check --as-cran
* win-builder, R-devel
* macOS builder, R-release
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1), macOS-latest and
  windows-latest (R-release)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission
  Package was archived on CRAN

  This is the resubmission of the archived package described above.
