#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List data_frame_phase_interval(List data, NumericVector weights) {
  if (data.size() != weights.size()) {
    stop("The length of data and weights should be the same.");
  }

    std::vector<double> expanded_values;
    std::vector<double> weights_values;
    std::vector<int> index_values;
    for (int i = 0; i < data.size(); ++i) {
        NumericVector x = data[i];
        if (x.size() == 1) {
            expanded_values.push_back(x[0]);
            weights_values.push_back(weights[i]);
            index_values.push_back(i + 1);
        } else if (x.size() == 2) {
            double a = x[0];
            double b = x[1];
            if (a == 0) {
                expanded_values.push_back(b);
                weights_values.push_back(weights[i]);
                index_values.push_back(0);
            } else if (std::isinf(b)) {
                expanded_values.push_back(a);
                weights_values.push_back(weights[i]);
                index_values.push_back(-1);
            } else {
                expanded_values.push_back(a);
                expanded_values.push_back(b);
                weights_values.push_back(weights[i]);
                weights_values.push_back(weights[i]);
                index_values.push_back(i + 1);
                index_values.push_back(i + 1);
            }
        }
    }

    int m = expanded_values.size();
    std::vector<int> ord(m);
    for (int i = 0; i < m; ++i) {
        ord[i] = i;
        if (index_values[i] == -1) {
            index_values[i] = m + 1;
        }
    }
    std::sort(ord.begin(), ord.end(), [&](int i, int j) {
        return expanded_values[i] < expanded_values[j];
    });
    std::vector<int> inv_ord(m);
    for (int i = 0; i < m; ++i) {
        inv_ord[ord[i]] = i + 1;
    }

    std::vector<int> z(m);
    int prev_index = -1;
    for (int i = 0; i < index_values.size(); ++i) {
        if (index_values[i] == prev_index) {
            z[i - 1] = inv_ord[i];
            z[i] = inv_ord[i - 1];
        } else if (index_values[i] == 0) {
            z[i] = 0;
        } else if (index_values[i] == m + 1) {
            z[i] = m + 1;
        } else {
            z[i] = inv_ord[i];
        }
        prev_index = index_values[i];
    }

    NumericVector reordered_values(m);
    NumericVector reordered_weights(m);
    IntegerVector reordered_z(m);
    for (int i = 0; i < m; ++i) {
        reordered_values[i] = expanded_values[ord[i]];
        reordered_weights[i] = weights_values[ord[i]];
        reordered_z[i] = z[ord[i]];
    }
    for (int i = m-1; i >= 1; --i) {
        reordered_values[i] -= reordered_values[i - 1];
    }
    double max_t = max(reordered_values);

    return List::create(
        Named("intervals") = reordered_values,
        Named("z") = reordered_z,
        Named("weights") = reordered_weights,
        Named("maxinterval") = max_t
    );
}