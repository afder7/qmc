#include "discrepancy.h"
#include <iostream>
#include <iomanip>
#include <boost/random/sobol.hpp>
#include <boost/random/faure.hpp>
#include <string>

void get_data(int dim, int num,
              std::vector<Point>& points, std::string type) {
  std::cout << "Discrepancy of " << num << " " << type
            << " points in dim " << dim
            << "\n";
  std::cout << std::setprecision(17) << "Star:   " << DiscrepancyRates::star_discrepancy(points) << "\n";
  // L2-star discrepancy measurement method used (Warnock's formula)
  // takes too long to compute for large values of num.
  if (num < 1e4) {
    std::cout << std::setprecision(17) << "L2:     " << DiscrepancyRates::L2_star_discrepancy(points)
              << "\n";
  }
}

void run_non_uniform(int dim, int num) {
  std::vector<Point> points(num, Point(dim));
  Point point_template(dim);

  points = std::vector<Point>(num, Point(dim));
  for (int i = 0; i < num; ++i) {
    point_template = Point(dim, 0.1);
    for (int j = 0; j < dim; ++j) {
      points[i][j] = static_cast<double>(point_template[j]) / SIZE_MAX;
    }
  }
  get_data(dim, num, points, "Non-Uniform");

  std::cout << "\n";
}

void run_tests(int dim, int num,
               boost::random::sobol sobol_gen,
               boost::random::faure faure_gen) {
  std::vector<Point> points(num, Point(dim));
  Point point_template(dim);

  points = std::vector<Point>(num, Point(dim));
  for (int i = 0; i < num; ++i) {
    sobol_gen.generate(point_template.begin(), point_template.end());
    for (int j = 0; j < dim; ++j) {
      // boost::sobol generates Sobol numbers from [0;SIZE_MAX],
      // so we map them into [0;1]
      points[i][j] = static_cast<double>(point_template[j]) / SIZE_MAX;
    }
  }
  get_data(dim, num, points, "Sobol no discard");

  sobol_gen.discard(100000);  // skips first n = 10000 points
  points = std::vector<Point>(num, Point(dim));
  for (int i = 0; i < num; ++i) {
    sobol_gen.generate(point_template.begin(), point_template.end());
    for (int j = 0; j < dim; ++j) {
      points[i][j] = static_cast<double>(point_template[j]) / SIZE_MAX;
    }
  }
  get_data(dim, num, points, "Sobol with discard");

  points = std::vector<Point>(num, Point(dim));
  for (int i = 0; i < num; ++i) {
    faure_gen.generate(point_template.begin(), point_template.end());
    for (int j = 0; j < dim; ++j) {
      points[i][j] = point_template[j];
    }
  }
  get_data(dim, num, points, "Faure no discard");

  points = std::vector<Point>(num, Point(dim));
  for (int i = 0; i < num; ++i) {
    for (int j = 0; j < dim; ++j) {
      points[i][j] = DiscrepancyRates::auxiliary::gen_random();
    }
  }

  get_data(dim, num, points, "Random");

  std::cout << "\n";
}

int main() {
  int dim = 15;
  int num = 1024;

  boost::random::sobol sobol_gen(dim);
  boost::random::faure faure_gen(dim);

  dim = 15;
  num = 1024;

  run_tests(dim, num, sobol_gen, faure_gen);

  run_non_uniform(dim, num);

  run_non_uniform(100, 1024);

  dim = 5;
  num = 131072;

  sobol_gen = boost::random::sobol(dim);
  faure_gen = boost::random::faure(dim);

  run_tests(dim, num, sobol_gen, faure_gen);

  dim = 15;
  num = 100000;

  sobol_gen = boost::random::sobol(dim);
  faure_gen = boost::random::faure(dim);

  run_tests(dim, num, sobol_gen, faure_gen);

  dim = 150;
  num = 65536;

  sobol_gen = boost::random::sobol(dim);
  faure_gen = boost::random::faure(dim);

  run_tests(dim, num, sobol_gen, faure_gen);
}
