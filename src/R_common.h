// // R library

#pragma once

#include <iostream>
#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "mexp.hpp"
#include "phfit.hpp"
#include "mapfit.hpp"

sci::matrix<double>* ccMatrix(const sci::matrix<double>& m, SEXP& v);

extern "C" {

	SEXP getListElement(SEXP list, const char *str);
	SEXP getSlot(SEXP obj, const char *str);
	sci::matrix<double>* createMatrix(SEXP m);
	sci::spmatrix<double>* createSpMatrix(SEXP m);

}

