# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`mapfit` is a CRAN R package (with an Rcpp/C++ core) for estimating parameters of phase-type
distributions (PH) and Markovian arrival processes (MAP) from point data, weighted point data,
grouped data (with missing/truncated bins), and theoretical density functions.
Estimation is done by EM algorithms (MLE) plus moment-matching (MM) methods.

## Commands

The host does not need R installed. Everything runs in the Docker image built from
`Dockerfile`; `Makefile` wraps the targets and mounts the repo as the invoking user so
generated files are not left owned by root.

```sh
make image      # build the image (once, and after editing Dockerfile)
make document   # roxygen2::roxygenise() -- regenerate NAMESPACE and man/
make test       # devtools::test()
make check      # R CMD build + R CMD check --as-cran
make readme     # devtools::build_readme() -- README.md is generated from README.Rmd
make win        # devtools::check_win_devel()    (result arrives by mail)
make mac        # devtools::check_mac_release()  (result arrives by mail)
make submit     # interactive R session for devtools::submit_cran()
make shell      # bash in the image
make clean
```

`make test` runs the whole suite; for a single file, `make shell` then
`Rscript -e 'devtools::test(filter = "CF1")'`. `Rcpp::compileAttributes()` regenerates
`src/RcppExports.cpp` and `R/RcppExports.R` after changing an `// [[Rcpp::export]]`.

CI (`.github/workflows/`) runs R-CMD-check on ubuntu (devel/release/oldrel-1), macOS and
Windows, plus codecov via `covr`.

Generated files that must never be hand-edited: `NAMESPACE`, `man/*.Rd`, `R/RcppExports.R`,
`src/RcppExports.cpp`, `README.md`.

### deformula

`deformula` is an Import that is not on CRAN at the moment: it was archived on 2026-02-15
(taking mapfit with it), and 0.1.3 was submitted on 2026-09-10. Until it is published, the
image and the CI workflow install it from `okamumu/deformula` on GitHub; both carry a
comment marking the line to drop once it is back.

## Architecture

### Three-layer structure

1. **User-facing S3 fitting functions** (`R/phfit.R`, `R/mapfit.R`, `R/phfit_mm*.R`):
   `phfit.point`, `phfit.group`, `phfit.density`, `phfit.3mom`, `mapfit.point`, `mapfit.group`.
   These are thin: they merge `...` into the defaults from `emoptions()` (`R/common.R`),
   build a data object, optionally call `model$init(data, options)`, call `model$emfit(...)`,
   and wrap the result in an S3 class (`phfit.result`, `mapfit.result`) carrying
   `model`, `llf`, `df`, `aic`, `iter`, `convergence`, `ctime`, `data`, `options`, `call`.
   The *algorithm is selected by the class of the model object passed in*, not by an argument.

2. **R6 model classes** (`R/model_*.R`) — one file per model family, each exposing
   parameter getters, `copy()`, `df()`, `emfit(data, options, ...)` and `init(data, ...)`:
   - `GPHClass` (`ph()`), `CF1Class` (`cf1()`, canonical form 1), `HErlangClass` (`herlang()`),
     plus `AHerlangClass`/`AERHMMClass` "aggregate" classes that search over Erlang shape
     vectors and keep the best fit.
   - `MAPClass` (`map()`, `mmpp()`), `ERHMMClass` (`erhmm()`), `GMMPPClass` (`gmmpp()`).
   Each `emfit` `switch`es on `class(data)` to pick the matching C++ entry point.

3. **C++ EM kernels** (`src/`), exposed through Rcpp attributes as `emfit_*` functions.
   `src/emfit.h` holds the generic `emfit(model, data, options, eres, work)` template loop:
   it repeatedly calls `estep(...)` / `mstep(...)`, which are overloaded per model in the
   headers (`phase_gen.h`, `phase_cf1.h`, `phase_herlang.h`, `map_gen.h`, `map_erhmm.h`,
   `map_gmmpp.h`). Adding a model means providing an `estep`/`mstep` overload plus a
   workspace/eres struct, and a thin `// [[Rcpp::export]]` wrapper in an `emfit_*.cpp`.
   `EMOptions` in `emfit.h` is the C++ mirror of the R `emoptions()` list — the two must
   stay in sync. Uniformization (`unif.h`), Poisson right-tail bounds (`poisson.h`),
   Gauss integration (`gauss_inte.h`) and GTH stationary solve (`gth.*`) are shared numerics.

Only the `.cpp` files listed in `src/Makevars` `OBJECTS` are compiled; `src/tests/*.cpp` are
standalone Rcpp harnesses for the C++ kernels and are deliberately **not** in that list.

### Data objects

`R/data_phase.R` and `R/data_map.R` build plain data frames tagged with S3 classes
`phase.time`, `phase.group`, `map.time`, `map.group` (constructors `data.frame.phase.time`
etc., with `print`/`mean` methods). These class tags are the dispatch key used inside every
`emfit` method, so a new data form needs both a constructor and a branch in each model's
`switch`.

### Shape search

`R/shape.R` enumerates Erlang shape vectors (`shape.all`, `shape.increment`, ...) for the
hyper-Erlang / ER-HMM aggregate classes; `options$shape.method`, `lbound`, `ubound` control it.
