/*
 erlang ph
*/

#pragma once

#include <cfloat>
#include "sci_spblas.hpp"

namespace mapfit {

	double erlang_pdf(int a, double b, double t);
	double erlang_lpdf(int a, double b, double t);
	double erlang_cdf(int a, double b, double t);
	double erlang_ccdf(int a, double b, double t);

}
