// Copyright 2023 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#pragma once

#include <random>

struct rnd_t {
  std::mt19937 rndgen;

  rnd_t() {
    std::random_device rd;
    std::mt19937 rndgen_t(rd());
    rndgen = rndgen_t;
  }

  explicit rnd_t(size_t seed) {
    std::mt19937 rndgen_t(seed);
    rndgen = rndgen_t;
  }

  std::uniform_real_distribution<float> unif_dist =
    std::uniform_real_distribution<float>(0.0f, 1.0f);

  int random_number(int n)    {
    if (n <= 1) return 0;
    return std::uniform_int_distribution<> (0, n - 1)(rndgen);
  }

  float uniform()    {
    return unif_dist(rndgen);
  }

  void set_seed(unsigned seed)    {
    std::mt19937 new_randomizer(seed);
    rndgen = new_randomizer;
  }

  bool bernouilli(double p) {
    std::bernoulli_distribution d(p);
    return(d(rndgen));
  }

  size_t binomial(int n, double p) {
    std::binomial_distribution<> d(n, p);
    return(d(rndgen));
  }

  int poisson(double lambda) {
    return std::poisson_distribution<int>(lambda)(rndgen);
  }

  double normal_positive(double m, double s) {
    std::normal_distribution<double> norm_dist(m, s);
    double  output = norm_dist(rndgen);
    while (output < 0.0) output = norm_dist(rndgen);
    return output;
  }
};
