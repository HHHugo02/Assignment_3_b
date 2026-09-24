#include <Rcpp.h>
#include <cmath>
#include <vector>

// One EM step for the augmented complete data (X, Z, W), written in C++.
// Same algorithm as create_em_aug_step() in R/em-algorithms.R.

// [[Rcpp::export]]
Rcpp::NumericVector em_aug_step_cpp(const Rcpp::NumericVector& x,
                                    const Rcpp::NumericVector& par,
                                    const Rcpp::NumericVector& nu) {
  const int N = x.size();
  const double p = par[0], mu1 = par[1], mu2 = par[2];
  const double s1 = par[3], s2 = par[4];
  const double nu1 = nu[0], nu2 = nu[1];

  // log of the constants p * Gamma((nu+1)/2) / (Gamma(nu/2) sqrt(pi nu) sigma)
  const double c1 = std::log(p) + std::lgamma((nu1 + 1) / 2) -
                    std::lgamma(nu1 / 2) - 0.5 * std::log(M_PI * nu1) -
                    std::log(s1);
  const double c2 = std::log1p(-p) + std::lgamma((nu2 + 1) / 2) -
                    std::lgamma(nu2 / 2) - 0.5 * std::log(M_PI * nu2) -
                    std::log(s2);

  std::vector<double> a1(N), a2(N);
  double N1 = 0, sa1 = 0, sa1x = 0, sa2 = 0, sa2x = 0;

  for (int i = 0; i < N; ++i) {
    const double r1 = x[i] - mu1, r2 = x[i] - mu2;
    const double d1 = r1 * r1 / (s1 * s1), d2 = r2 * r2 / (s2 * s2);
    const double l1 = c1 - 0.5 * (nu1 + 1) * std::log1p(d1 / nu1);
    const double l2 = c2 - 0.5 * (nu2 + 1) * std::log1p(d2 / nu2);
    const double g = 1.0 / (1.0 + std::exp(l2 - l1)); // E-step
    a1[i] = g * (nu1 + 1) / (nu1 + d1);
    a2[i] = (1 - g) * (nu2 + 1) / (nu2 + d2);
    N1 += g;
    sa1 += a1[i];
    sa1x += a1[i] * x[i];
    sa2 += a2[i];
    sa2x += a2[i] * x[i];
  }

  const double mu1_new = sa1x / sa1, mu2_new = sa2x / sa2;
  double ss1 = 0, ss2 = 0;
  for (int i = 0; i < N; ++i) {
    ss1 += a1[i] * (x[i] - mu1_new) * (x[i] - mu1_new);
    ss2 += a2[i] * (x[i] - mu2_new) * (x[i] - mu2_new);
  }

  return Rcpp::NumericVector::create(N1 / N, mu1_new, mu2_new,
                                     std::sqrt(ss1 / N1),
                                     std::sqrt(ss2 / (N - N1)));
}
