#include <RcppArmadillo.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp14)]]

#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <numeric>
#include <limits>
#include <math.h>
#include <random>
#include <cstdint>
#include <cmath>
#include <algorithm>


using namespace Rcpp;
using namespace arma;


// [[Rcpp::export(rng = false)]]
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

// [[Rcpp::export(rng = false)]]
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
    MyParallelTask(const arma::vec& G_bar, const arma::vec& W, arma::vec& likel_ij, double alpha, int n, int n_0)
      // std::vector<double>& f1_values, std::vector<double>& f2_values)
      : G_bar(G_bar), W(W), likel_ij(likel_ij), alpha(alpha), n(n), n_0(n_0) {}
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
  MyParallelTask task(G_bar, W, likel_ij, alpha, n, n_0);
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
arma::imat rmvbern_cpp(const int n, const arma::vec & p) {
  arma::mat U = arma::randu(n, p.n_elem);
  U.each_row() -= p.t();
  return arma::conv_to<arma::imat>::from(U < 0.0);
}

// [[Rcpp::export]]
arma::imat rmvbern(const int n, const arma::vec& p) {
  return rmvbern_cpp(n, p);
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

// [[Rcpp::export(rng = false)]]
double triangles(const std::vector<uint64_t>& a,int N,int nw) {
  uint64_t total=0;
  for(int i=0;i<N;i++) for(int w=0;w<nw;w++) {
    uint64_t v=a[i*nw+w];
    while(v) {
      int j=64*w+__builtin_ctzll(v);v&=v-1;
      if(j>i && j<N) for(int k=0;k<nw;k++) total+=__builtin_popcountll(a[i*nw+k]&a[j*nw+k]);
    }
  }
  return double(total)/3.;
}

//' @export
// [[Rcpp::export(rng = false)]]
arma::vec exact_graph_stats(const arma::ivec& C) {

  int M = C.n_elem;
  int N = static_cast<int>((1.0 + std::sqrt(1.0 + 8.0 * M)) / 2.0);

  int nw = (N + 63) / 64;


  std::vector<uint64_t> a(N * nw, 0);
  std::vector<int> degree(N, 0);

  int e = 0;
  double edges = 0.0;
  double d2 = 0.0;

  for (int j = 0; j < N; j++) {
    for (int i = j + 1; i < N; i++, e++) {

      if (C[e] == 1) {

        a[i * nw + j / 64] |= uint64_t(1) << (j % 64);
        a[j * nw + i / 64] |= uint64_t(1) << (i % 64);

        degree[i]++;
        degree[j]++;
        edges++;
      }
    }
  }

  for (int d : degree)
    d2 += static_cast<double>(d) * d;

  arma::vec out(3);

  out(0) = edges;
  out(1) = d2;
  out(2) = triangles(a, N, nw);

  return out;
}

//' @export
 // [[Rcpp::export(rng = false)]]
double clustering(double a,double edges,double d2,double T,int N) {
  double eta=1-2*a,choose3=double(N)*(N-1)*(N-2)/6.;
  double et= a*a*a*choose3+a*a*eta*(N-2)*edges+.5*a*eta*eta*(d2-2*edges)+eta*eta*eta*T;
  double ey=3*a*a*choose3+eta*(2*(N-1)*a-1)*edges+.5*eta*eta*d2;
  return 3*et/ey;
}


// [[Rcpp::export(rng = false)]]
Rcpp::List fit_eb_cpp(const arma::ivec& counts, int n) {

  const double eps = 1e-9;
  const int M = counts.n_elem;

  // Histogram of edge counts
  arma::ivec h(n + 1, arma::fill::zeros);

  for (int i = 0; i < M; i++) {
    h[counts[i]]++;
  }

  // Objective function
  auto objective = [&](double alpha) {

    double la = std::log(alpha);
    double l1a = std::log1p(-alpha);

    double out = 0.0;

    for (int s = 0; s <= n; s++) {

      if (h[s] == 0) continue;

      double a = s * la + (n - s) * l1a;
      double b = (n - s) * la + s * l1a;

      double mx = std::max(a, b);

      double log_sum = mx +
        std::log1p(std::exp(-std::abs(a - b)));

      out -= h[s] * (log_sum - std::log(2.0));
    }

    return out;
  };

  // Golden-section optimization
  double lower = eps;
  double upper = 0.5 - eps;
  const double tol = 1e-11;
  const double gr = (std::sqrt(5.0) - 1.0) / 2.0;

  double x1 = upper - gr * (upper - lower);
  double x2 = lower + gr * (upper - lower);

  double f1 = objective(x1);
  double f2 = objective(x2);

  while (upper - lower > tol) {

    if (f1 < f2) {
      upper = x2;
      x2 = x1;
      f2 = f1;
      x1 = upper - gr * (upper - lower);
      f1 = objective(x1);
    } else {
      lower = x1;
      x1 = x2;
      f1 = f2;
      x2 = lower + gr * (upper - lower);
      f2 = objective(x2);
    }
  }

  double alpha_opt = 0.5 * (lower + upper);

  // Majority-rule center
  arma::ivec C_hat = arma::conv_to<arma::ivec>::from(2 * counts >= n);

  return Rcpp::List::create(
    Rcpp::_["alpha_hat"] = alpha_opt,
    Rcpp::_["C_hat"] = C_hat.t(),
    Rcpp::_["objective"] = objective(alpha_opt)
  );
}


// [[Rcpp::export]]
List optimized_gibbs(arma::ivec counts,
                     int n,
                     int keep = 1000,
                     int burn = 1000,
                     int thin = 10,
                     uint64_t seed = 123456
){

  int M = counts.n_elem;
  int N = static_cast<int>((1.0 + std::sqrt(1.0 + 8.0 * M)) / 2.0);

  int nw = (N + 63) / 64;
  std::int64_t baseS = 0;

  std::vector<std::vector<int>> groups(n + 1);
  std::vector<int> active;

  arma::ivec u(M), v(M), baseC(M);
  arma::ivec degree0(N, arma::fill::zeros);
  arma::ivec kflip(n + 1, arma::fill::zeros);

  // Bit-packed adjacency matrix
  std::vector<uint64_t> base(N * nw, 0);

  arma::vec cm(M, arma::fill::zeros);
  arma::vec alphas(keep);
  arma::vec ecc(keep);
  arma::vec density(keep);
  arma::vec nedges(keep);
  arma::vec d2chain(keep);
  arma::vec tri(keep);


  int incremental_draws = 0;
  int e = 0;
  double edges0 = 0.0;


  for (int j = 0; j < N; j++) {
    for (int i = j + 1; i < N; i++, e++) {

      u[e] = i;
      v[e] = j;

      int s = counts[e];

      groups[s].push_back(e);

      baseS += std::min(s, n - s);

      baseC[e] = (2 * s >= n);
      cm[e] = baseC[e] * keep;

      if (baseC[e]) {

        base[i * nw + j / 64] |= uint64_t(1) << (j % 64);
        base[j * nw + i / 64] |= uint64_t(1) << (i % 64);

        degree0[i]++;
        degree0[j]++;

        edges0++;
      }
    }
  }

  for (int s = 0; s <= n; s++) {
    if (!groups[s].empty())
      active.push_back(s);
  }

  double tri0 = triangles(base, N, nw);
  double alpha = 0.5;
  int retained = 0;


  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> unif01(0.0, 1.0);


  std::int64_t niter =
    static_cast<std::int64_t>(burn) +
    static_cast<std::int64_t>(thin) * keep;

  for (std::int64_t it = 1; it <= niter; it++){

    if (it % 100 == 0)
      checkUserInterrupt();

    double logratio = std::log1p(-alpha) - std::log(alpha);

    std::int64_t S = baseS;
    int flips = 0;

    for (int s : active) {

      int d = std::abs(n - 2 * s);

      double x = -d * logratio;

      double q = x >= 0.0
      ? 1.0 / (1.0 + std::exp(-x))
        : std::exp(x) / (1.0 + std::exp(x));

      std::binomial_distribution<int> binom(
          static_cast<int>(groups[s].size()), q
      );

      kflip[s] = binom(rng);


      S += static_cast<std::int64_t>(d) * kflip[s];
      flips += kflip[s];
    }

    double aa = 1. + S;
    double bb = 1. + double(n) * M - S;

    double logcdf = R::pbeta(0.5, aa, bb, true, true);

    alpha = R::qbeta(
      std::log(R::runif(0, 1)) + logcdf,
      aa, bb, true, true
    );


    if (it > burn && (it - burn) % thin == 0) {

      std::vector<uint64_t> a = base;
      arma::ivec degrees = degree0;

      double T = tri0;
      double edges = edges0;

      bool incremental = (flips < M / 8);

      if (incremental)
        incremental_draws++;

      for (int s : active) {

        auto &g = groups[s];

        int size = g.size();
        int k = kflip[s];

        for (int j = 0; j < k; j++) {

          int h = j + std::uniform_int_distribution<int>(0, size - j - 1)(rng);

          std::swap(g[j], g[h]);

          int edge = g[j];

          int i = u[edge];
          int l = v[edge];

          int sign = baseC[edge] ? -1 : 1;

          cm[edge] += sign;

          edges += sign;

          degrees[i] += sign;
          degrees[l] += sign;

          if (incremental) {

            int common = 0;

            for (int w = 0; w < nw; w++) {

              common += __builtin_popcountll(
                a[i * nw + w] & a[l * nw + w]
              );
            }

            T += sign * common;
          }

          // Update bit-packed adjacency matrix
          a[i * nw + l / 64] ^= uint64_t(1) << (l % 64);
          a[l * nw + i / 64] ^= uint64_t(1) << (i % 64);
        }
      }


      if (!incremental)
        T = triangles(a, N, nw);

      double d2 = 0.0;

      for (int d : degrees)
        d2 += double(d) * d;

      alphas[retained] = alpha;

      ecc[retained] = clustering(alpha, edges, d2, T, N);

      density[retained] =
        (M * alpha + (1 - 2 * alpha) * edges) / M;

      tri[retained] = T;
      nedges[retained] = edges;
      d2chain[retained] = d2;

      retained++;
    }
  }


  cm /= keep;

  return List::create(
    _["C_mean"] = cm.t(),
    _["alpha_chain"] = alphas,
    _["ecc_dist"] = ecc,
    _["density_dist"] = density,
    _["triangle_count"] = tri,
    _["degree_squared_sum"] = d2chain,
    _["center_edge_count"] = nedges
  );
}
