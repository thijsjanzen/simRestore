#ifndef GENOME_TOOLS_HPP
#define GENOME_TOOLS_HPP

#include "Duck.h"
#include "rand_t.h"
#include "main_functions.h"
#include "analysis.h"

#include <Rcpp.h>
using namespace Rcpp;


std::vector<int> detect_ancestry(const double& G,
                                 const NumericVector& markers) {
  std::vector<int> output(markers.size(), G);
  return output;
}


std::vector<int> detect_ancestry(const std::vector< junction >& G,
                                 const NumericVector& markers) {
  std::vector<int> output(markers.size());

  int j = 0;
  for(int i = 0; i < markers.size(); ++i) {
    float focalPos = markers[i];
    for(; j <= (G.size()-1); ++j) {
      float left = G[j].pos;
      float right = G[j+1].pos;
      if(left <= focalPos && right >= focalPos) {
        output[i] = ((int)G[j].right);
        break;
      }
    }
    j-=5; //just in case
    if(j < 0) j = 0;
  }
  return output;
}



NumericMatrix convert_genomes(const std::vector< genome_data > genomes) {

  NumericMatrix output_genomes(genomes.size(), 6);
  if(!genomes.empty()) {
    for(int i = 0; i < genomes.size(); ++i) {
      output_genomes(i, 0) = genomes[i].t;
      output_genomes(i, 1) = genomes[i].r;
      output_genomes(i, 2) = genomes[i].ind;
      output_genomes(i, 3) = genomes[i].pos;
      output_genomes(i, 4) = genomes[i].allele_1;
      output_genomes(i, 5) = genomes[i].allele_2;
    }
  }
  return output_genomes;
}

NumericMatrix convert_genomes(const std::vector< Duck_emp>& v) {
  int nrow = v.size();
  int ncol = v[0].get_chromosome1().size() * 2; // both chromosomes are merged

  NumericMatrix output(nrow, ncol);

  for(int i = 0; i < v.size(); ++i) {
    for (int j = 0; j < v[i].get_chromosome1().size(); ++j) {
      output(i, j * 2) = v[i].get_chromosome1()[j];
      output(i, j * 2 +1) = v[i].get_chromosome2()[j];
    }
  }
  return output;
}

NumericVector get_sexes(const std::vector< Duck_emp>& v) {
  NumericVector output(v.size());
  for (int i = 0; i < v.size(); ++i) {
    output(i) = v[i].get_sex();
  }
  return(output);
}

std::vector< Duck_emp > NumericMatrix_to_emp_bird(const Rcpp::NumericMatrix& mat,
                                                  const emp_genome& empgen,
                                                  rnd_t& rndgen) {
  std::vector< Duck_emp > pop;
  for (int i = 0; i < mat.nrow(); i += 2) {
 //   Rcpp::Rcout << i << "\n";
    Rcpp::NumericVector c1 = mat(i, _);
    Rcpp::NumericVector c2 = mat(i + 1, _);
    std::vector< int > chrom1(c1.begin(), c1.end());
    std::vector< int > chrom2(c2.begin(), c2.end());
    Duck_emp to_add(std::make_tuple(chrom1, empgen),
                    std::make_tuple(chrom2, empgen),
                    rndgen);
    pop.push_back(to_add);
  }
  return(pop);
}

void numericmatrix_to_vector(const Rcpp::NumericMatrix& m,
                             std::vector< std::vector< double >>& v) {

  v = std::vector< std::vector< double>>(m.nrow());
  for (int i = 0; i < m.nrow(); ++i) {
    std::vector<double> row(m.ncol());
    for (int j = 0; j < m.ncol(); ++j) {
      row[j] = m(i, j);
    }
    v[i] = row;
  }
  return;
}

template <class bird>
std::array<double, 3> calc_freq_hawaii(std::vector< bird >& f,
                                       std::vector< bird >& m) {
  std::array<double, 3> avg_freq = {0.0, 0.0, 0.0};

  for (auto& i : f) {
    double freq_i = i.get_freq_hawaii();
    avg_freq[2] += freq_i;
    avg_freq[0] += freq_i;
  }
  for (auto& i : m) {
    double freq_i = i.get_freq_hawaii();
    avg_freq[1] += freq_i;
    avg_freq[0] += freq_i;
  }

  avg_freq[0] *= 1.0 / (f.size() + m.size());
  avg_freq[1] *= 1.0 / (m.size());
  avg_freq[2] *= 1.0 / (f.size());

  return avg_freq;
}


//' calculate ancestry
//' @param pop pop
//' @param anc_info_R anc info
//' @export
// [[Rcpp::export]]
double calc_anc_cpp(const Rcpp::NumericMatrix& pop,
                    const Rcpp::NumericMatrix& anc_info_R) {

//  Rcpp::Rcout << "loading anc_info_cpp\n"; force_output();
  std::vector< std::vector< double >> anc_info_cpp;

  numericmatrix_to_vector(anc_info_R, anc_info_cpp);

  emp_genome empgen(anc_info_cpp);

  rnd_t rndgen;
 // Rcpp::Rcout << "loading genomes\n"; force_output();
  std::vector<Duck_emp> population = NumericMatrix_to_emp_bird(pop,
                                                                empgen,
                                                                rndgen);

 // Rcpp::Rcout << "calculating ancestries\n"; force_output();
  double avg_anc = 0.0;
  for(auto& i : population) {
    avg_anc += i.get_freq_hawaii();
  }
  avg_anc *= 1.0/ (population.size());
  return(avg_anc);
}



#endif
