#include <Rcpp.h>
#include <cmath>

// The building blocks of the log-scale class "tmix_log" written as plain C++
// loops, for the class "tmix_cpp". Each function makes one pass over the data
// and keeps only a few numbers per observation, so no vectors of length N are
// created. par = (p, mu1, mu2, sigma1, sigma2) and log_const = c(nu) as in R.

// Everything in a = log(p f1(x)) and b = log((1 - p) f2(x)) that does not
// depend on x. It is computed once per call: the logs of p and sigma would
// otherwise be recomputed N times, which made a first version slower than R.
struct Comp {
  double k1, k2;   // log p + c(nu1) - log sigma1 and log(1 - p) + c(nu2) - log sigma2
  double h1, h2;   // (nu + 1) / 2
  double s1, s2;   // 1 / (nu sigma^2)
  double mu1, mu2;

  Comp(const Rcpp::NumericVector& par, const Rcpp::NumericVector& nu, const Rcpp::NumericVector& log_const)
  {
    k1 = std::log(par[0]) + log_const[0] - std::log(par[3]);
    k2 = std::log1p(-par[0]) + log_const[1] - std::log(par[4]);
    h1 = (nu[0] + 1.0) / 2.0;
    h2 = (nu[1] + 1.0) / 2.0;
    s1 = 1.0 / (nu[0] * par[3] * par[3]);
    s2 = 1.0 / (nu[1] * par[4] * par[4]);
    mu1 = par[1];
    mu2 = par[2];
  }

  void ab(double x, double& a, double& b) const
  {
    const double r1 = x - mu1, r2 = x - mu2;
    a = k1 - h1 * std::log1p(r1 * r1 * s1);
    b = k2 - h2 * std::log1p(r2 * r2 * s2);
  }
};

// log-likelihood, log(e^a + e^b) with the log-sum-exp trick
// [[Rcpp::export]]
double loglik_cpp(Rcpp::NumericVector x, Rcpp::NumericVector par, Rcpp::NumericVector nu,
                  Rcpp::NumericVector log_const)
{
  const Comp comp(par, nu, log_const);
  double s = 0.0, a, b;
  for (int i = 0; i < x.size(); ++i) {
    comp.ab(x[i], a, b);
    const double m = std::max(a, b);
    s += m + std::log1p(std::exp(-std::abs(a - b)));
  }
  return s;
}

// E-step, gamma = 1 / (1 + e^(b - a))
// [[Rcpp::export]]
Rcpp::NumericVector e_step_cpp(Rcpp::NumericVector x, Rcpp::NumericVector par, Rcpp::NumericVector nu,
                               Rcpp::NumericVector log_const)
{
  const Comp comp(par, nu, log_const);
  Rcpp::NumericVector gam(x.size());
  double a, b;
  for (int i = 0; i < x.size(); ++i) {
    comp.ab(x[i], a, b);
    gam[i] = 1.0 / (1.0 + std::exp(b - a));
  }
  return gam;
}

// Q = sum gam a + (1 - gam) b
// [[Rcpp::export]]
double Q_cpp(Rcpp::NumericVector x, Rcpp::NumericVector par, Rcpp::NumericVector gam,
             Rcpp::NumericVector nu, Rcpp::NumericVector log_const)
{
  const Comp comp(par, nu, log_const);
  double s = 0.0, a, b;
  for (int i = 0; i < x.size(); ++i) {
    comp.ab(x[i], a, b);
    s += gam[i] * a + (1.0 - gam[i]) * b;
  }
  return s;
}

// Gradient of Q, the same formula as grad_Q.tmix() in R
// [[Rcpp::export]]
Rcpp::NumericVector grad_Q_cpp(Rcpp::NumericVector x, Rcpp::NumericVector par, Rcpp::NumericVector gam,
                               Rcpp::NumericVector nu)
{
  const double s1 = par[3], s2 = par[4];
  const double is1 = 1.0 / (s1 * s1), is2 = 1.0 / (s2 * s2);
  double g_sum = 0.0, d_mu1 = 0.0, d_mu2 = 0.0, d_s1 = 0.0, d_s2 = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    const double r1 = x[i] - par[1];
    const double r2 = x[i] - par[2];
    const double w1 = (nu[0] + 1.0) / (nu[0] + r1 * r1 * is1);
    const double w2 = (nu[1] + 1.0) / (nu[1] + r2 * r2 * is2);
    g_sum += gam[i];
    d_mu1 += gam[i] * w1 * r1;
    d_mu2 += (1.0 - gam[i]) * w2 * r2;
    d_s1 += gam[i] * (w1 * r1 * r1 * is1 - 1.0);
    d_s2 += (1.0 - gam[i]) * (w2 * r2 * r2 * is2 - 1.0);
  }
  const double N = x.size();
  return Rcpp::NumericVector::create(g_sum / par[0] - (N - g_sum) / (1.0 - par[0]),
                                     d_mu1 * is1, d_mu2 * is2, d_s1 / s1, d_s2 / s2);
}
