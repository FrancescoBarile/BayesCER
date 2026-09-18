#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <numeric>
#include <limits>
#include <math.h>


using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
double marg_lik_CER_gCER(double alpha,
                         int n,
                         int n_0,
                         const arma::vec& G_bar,
                         const arma::vec& W) {

  int M = G_bar.n_elem;

  double log_alpha = std::log(alpha);
  double log_1m_alpha = std::log1p(-alpha);

  double out = 0.0;

  for (int i = 0; i < M; ++i) {

    double e1 = n * G_bar[i];
    double e2 = n_0 * (1.0 - 2.0 * W[i]);
    double e3 = n * (1.0 - G_bar[i]) + e2;

    double a1 = e1 * log_alpha + e3 * log_1m_alpha;
    double a2 = e3 * log_alpha + e1 * log_1m_alpha;

    double m1 = std::max(a1, a2);
    double f1 = m1 + std::log(std::exp(a1 - m1) + std::exp(a2 - m1));

    double b1 = e2 * log_alpha;
    double b2 = e2 * log_1m_alpha;

    double m2 = std::max(b1, b2);
    double f2 = m2 + std::log(std::exp(b1 - m2) + std::exp(b2 - m2));

    out += f1 - f2;
  }

  return -out;
}

//[[Rcpp::export]]
double marg_lik_CER_gCER_par(double alpha,
                             int n,
                             int n_0,
                             arma::vec& G_bar,
                             arma::vec& W){

  int M = G_bar.n_elem;  // Use n_elem for armadillo vectors

  arma::vec likel_ij(M, arma::fill::zeros);  // Vector to store the likelihoods

  // std::vector<double> f1_values(M, 0.0);  // Vector to store f1 values
  // std::vector<double> f2_values(M, 0.0);  // Vector to store f2 values

  // Define a functor that will be used for parallel execution
  struct MyParallelTask : public RcppParallel::Worker {
    const arma::vec& G_bar;
    const arma::vec& W;
    arma::vec& likel_ij;  // Output vector to store the likelihoods
    const double alpha;
    const int n;
    const int n_0;

    // std::vector<double>& f1_values;  // Reference to store f1 values
    // std::vector<double>& f2_values;  // Reference to store f2 values

    // Constructor to initialize the member variables
    MyParallelTask(const arma::vec& G_bar, const arma::vec& W, double alpha, int n, int n_0, arma::vec& likel_ij)
      // std::vector<double>& f1_values, std::vector<double>& f2_values)
      : G_bar(G_bar), W(W), alpha(alpha), n(n), n_0(n_0), likel_ij(likel_ij) {}
    // f1_values(f1_values), f2_values(f2_values) {}

    // Parallel worker class
    void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {
        // Computation
        double e1 = n * G_bar[i];
        double e2 = n_0 * (1 - 2 * W[i]);
        double e3 = n * (1 - G_bar[i]) + e2;

        // double f1 = log( pow(alpha, e1) * pow(1 - alpha, e3) + pow(alpha, e3) * pow(1 - alpha, e1) );
        // double f2 = log( pow(alpha, e2) + pow(1 - alpha, e2) );

        // log(alpha^e1 * (1-alpha)^e3 + alpha^e3 * (1-alpha)^e1)

        double a1 = e1 * std::log(alpha) + e3 * std::log1p(-alpha);
        double a2 = e3 * std::log(alpha) + e1 * std::log1p(-alpha);

        double m1 = std::max(a1, a2);

        double f1 = m1 + std::log(std::exp(a1 - m1) + std::exp(a2 - m1));


        // log(alpha^e2 + (1-alpha)^e2)

        double b1 = e2 * std::log(alpha);
        double b2 = e2 * std::log1p(-alpha);

        double m2 = std::max(b1, b2);

        double f2 = m2 + std::log(std::exp(b1 - m2) + std::exp(b2 - m2));

        likel_ij[i] = f1 - f2;  // Store the result in the output vector

        // Store f1 and f2 values for later printing
        // f1_values[i] = f1;
        // f2_values[i] = f2;
      }
    }
  };

  // Create the worker for parallel execution
  MyParallelTask task(G_bar, W, alpha, n, n_0, likel_ij);
  // f1_values, f2_values);

  // Parallelize the loop from 0 to M (the length of the vectors)
  RcppParallel::parallelFor(0, M, task);

  // Print f1 and f2 values after parallel execution
  // for (int i = 0; i < M; ++i) {
  //   Rcpp::Rcout << "f1[" << i << "] = " << f1_values[i] << ", f2[" << i << "] = " << f2_values[i] << std::endl;
  // }

  // Use arma::sum() to calculate the sum of the log-likelihoods
  double result = - arma::sum(likel_ij);  // Sum of all elements in likel_ij (negative to maximize)

  return result;
}


// [[Rcpp::export]]
arma::imat rmvbern_CER(const int n,
                       const arma::ivec& C,
                       const double alpha) {

  if (alpha <= 0.0 || alpha >= 0.5)
    Rcpp::stop("'alpha' must be in (0,1/2).");

  arma::uword M = C.n_elem;

  arma::vec eta = 2.0 * arma::conv_to<arma::vec>::from(C) - 1.0;

  double odds = alpha/(1.0-alpha);

  arma::vec p = 1.0 /
    (1.0 + arma::exp(eta * std::log(odds)));

  arma::imat G(n, M);

  for (int i = 0; i < n; i++) {

    arma::vec U = arma::randu(M);

    U -= p;

    G.row(i) = arma::conv_to<arma::irowvec>::from(U < 0.0);
  }

  return G;
}
// [[Rcpp::export]]
arma::rowvec rmvbern_CER_sum(const int n,
                              const arma::ivec& C,
                              const double alpha) {

  if (alpha <= 0.0 || alpha >= 0.5)
    Rcpp::stop("'alpha' must be in (0,1/2).");

  arma::uword M = C.n_elem;

  // eta = 2C - 1
  arma::vec eta = 2.0 * arma::conv_to<arma::vec>::from(C) - 1.0;

  double odds = alpha/(1.0-alpha);

  arma::vec p = 1.0 /
    (1.0 + arma::exp(eta * std::log(odds)));

  // store only accumulated sums
  arma::vec G_sum(M, arma::fill::zeros);

  for (int i = 0; i < n; i++) {

    arma::vec U = arma::randu(M);

    G_sum += arma::conv_to<arma::vec>::from(U < p);
  }

  return G_sum.t();   // transpose: M x 1 --> 1 x M
}

// Useful to sample from a GCER
// [[Rcpp::export]]
arma::imat rmvbern(const int n, const arma::vec & p) {
  arma::mat U = arma::randu(n, p.n_elem);
  U.each_row() -= p.t();
  return arma::conv_to<arma::imat>::from(U < 0.0);
}

// [[Rcpp::export]]
double expected_transitivity_CER_approx_cpp(double alpha,
                                            const arma::mat& C_mat) {

  int N = C_mat.n_cols;

  // Degree vector
  arma::vec degree_vec = sum(C_mat, 0).t();

  // Number of edges
  double degree_C = accu(degree_vec) / 2.0;

  // Number of triangles
  double T_C = trace(C_mat * C_mat * C_mat) / 6.0;

  double choose3 = static_cast<double>(N) * (N - 1) * (N - 2) / 6.0;
  double eta = 1.0 - 2.0 * alpha;

  double E_T_G =
    std::pow(alpha, 3.0) * choose3 +
    alpha * alpha * eta * (N - 2.0) * degree_C +
    0.5 * alpha * eta * eta *
    (accu(square(degree_vec)) - 2.0 * degree_C) +
    std::pow(eta, 3.0) * T_C;

  double num = 3.0 * E_T_G;

  double den =
    3.0 * alpha * alpha * choose3 +
    eta * (2.0 * (N - 1.0) * alpha - 1.0) * degree_C +
    0.5 * eta * eta * accu(square(degree_vec));

  return num / den;
}
