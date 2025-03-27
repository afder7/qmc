#pragma once
#include <cmath>
#include <limits>
#include <vector>

// d-dimensional point type definition
using Point = std::vector<double>;

namespace DiscrepancyRates {
  // auxiliary functions
  namespace auxiliary {

    // Indicates if x is in anchored box [0;a)
    int indicator(int dim, Point& x, Point& a) {
      for (int j = 0; j < dim; ++j) {
        if (x[j] >= a[j]) {
          return 0;
        }
      }

      return 1;
    }

    // Generates pseudo random point from [0;1]
    double gen_random() { return static_cast<double>(std::rand()) / RAND_MAX; }
    
    // Local discrepancy delta(a).
    double local_discrepancy(int dim, Point& a, std::vector<Point>& points) {
      // todo (?) : optimize indicator calculation assuming a values are ascending
      int n = static_cast<int>(points.size());

      int ind_sum = 0;
      for (Point& x : points) {
        ind_sum += indicator(dim, x, a);
      }

      double vol_a = 1;
      for (double a_j : a) {
        vol_a *= a_j;
      }

      double result = (1.0l / n) * ind_sum - vol_a;
      return result;
    }
  }

  // Star discrepancy is calculated using Simple Monte Carlo.
  // Calculation is computed according to definition, see Owen (2023).
  // sample_size could be increased for higher precision or lowered for faster computation.
  double star_discrepancy(std::vector<Point>& points, int sample_size = 1000) {
    int dim = static_cast<int>(points[0].size());

    double result = 0;
    for (int i = 0; i < sample_size; ++i) {
      Point a(dim);
      for (int k = 0; k < dim; ++k) {
        a[k] = auxiliary::gen_random();
      }
      result = std::max(result,
                        std::abs(auxiliary::local_discrepancy(dim, a, points)));
    }

    return result;
  }

  // L2-star discrepancy is calculated with Warnock's (1972) formula.
  // Asymptotics: O(n^2 * d), so it is impractical for large n.
  double L2_star_discrepancy(std::vector<Point>& points) {
    int n = static_cast<int>(points.size());
    int d = static_cast<int>(points[0].size());

    double term1 = 0;
    for (int i = 0; i < n; ++i) {
      double product = 1;
      for (int j = 0; j < d; ++j) {
        product *= std::pow((1 - points[i][j]), 2) / 2;
      }
      term1 += product;
    }

    double term2 = 0;
    for (int i = 0; i < n; ++i) {
      double subsum = 0;
      for (int k = 0; k < n; ++k) {
        double product = 1;
        for (int j = 0; j < d; ++j) {
          product *= (1 - std::max(points[i][j], points[k][j]));
        }
        subsum += product;
      }
      term2 += subsum;
    }

    double result = std::pow(1.0f / 3, d);
    result -= (2.0f / n) * term1;
    result += (1.0f / std::pow(n, 2)) * term2;

    return std::sqrt(result);
  }
};
