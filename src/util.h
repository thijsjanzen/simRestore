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

#include <vector>
#include <Rcpp.h>

void vector_to_numericmatrix(const std::vector< std::vector< double >>& v,
                             Rcpp::NumericMatrix* m) {

  int n_rows = v.size();
  int n_cols = v[0].size();
  *m = Rcpp::NumericMatrix(n_rows, n_cols);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      (*m)(i, j) = v[i][j];
    }
  }
  return;
}

void numericmatrix_to_vector(Rcpp::NumericMatrix* m,
                             std::vector< std::vector< double >>* v) {

  *v = std::vector< std::vector< double> >((*m).nrow(),
                                           std::vector<double>((*m).ncol(), 0.0));
  for (size_t i = 0; i < (*m).nrow(); ++i) {
    std::vector<double> row((*m).ncol(), 0.0);
    for (size_t j = 0; j < (*m).ncol(); ++j) {
      row[j] = (*m)(i, j);
    }
    (*v)[i] = row;
  }
  return;
}
