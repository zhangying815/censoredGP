// Monte Carlo SOV estimator for P(Y_j > C_j for all j), Y ~ N(mu, Sigma).
// Compile in R with:
  //   Rcpp::sourceCpp("mc_sov_censored_nd.cpp")

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

namespace {
  
  std::vector<double> threshold_vector(const NumericVector& C, int r) {
    std::vector<double> Cvec(r);
    if (C.size() == 1) {
      for (int j = 0; j < r; ++j) Cvec[j] = C[0];
    } else if (C.size() == r) {
      for (int j = 0; j < r; ++j) Cvec[j] = C[j];
    } else {
      stop("C must have length 1 or length(mu).");
    }
    return Cvec;
  }
  
  NumericMatrix chol_lower_cpp(const NumericMatrix& Sigma) {
    int r = Sigma.nrow();
    NumericMatrix L(r, r);
    
    for (int i = 0; i < r; ++i) {
      for (int j = 0; j <= i; ++j) {
        double s = Sigma(i, j);
        for (int k = 0; k < j; ++k) {
          s -= L(i, k) * L(j, k);
        }
        
        if (i == j) {
          if (!R_finite(s) || s <= 0.0) {
            stop("Sigma is not positive definite.");
          }
          L(i, j) = std::sqrt(s);
        } else {
          L(i, j) = s / L(j, j);
        }
      }
    }
    
    return L;
  }
  
  double log_mean_exp(const std::vector<double>& log_w) {
    double max_log_w = R_NegInf;
    int n = static_cast<int>(log_w.size());
    
    for (int i = 0; i < n; ++i) {
      if (log_w[i] > max_log_w) max_log_w = log_w[i];
    }
    if (!R_finite(max_log_w)) return R_NegInf;
    
    double sum_scaled = 0.0;
    for (int i = 0; i < n; ++i) {
      sum_scaled += std::exp(log_w[i] - max_log_w);
    }
    
    return max_log_w + std::log(sum_scaled / static_cast<double>(n));
  }
  
} // namespace

// [[Rcpp::export]]
double mc_sov_censored_nd_cpp(NumericVector mu,
                              NumericMatrix Sigma,
                              NumericVector C,
                              int M = 10000,
                              double eps = 1e-300) {
  int r = mu.size();
  if (r < 1) stop("mu must have positive length.");
  if (Sigma.nrow() != r || Sigma.ncol() != r) {
    stop("Sigma dimensions must match length(mu).");
  }
  if (M <= 0) stop("M must be a positive integer.");
  
  std::vector<double> Cvec = threshold_vector(C, r);
  NumericMatrix L = chol_lower_cpp(Sigma);
  std::vector<double> log_w(M, 0.0);
  std::vector<double> z(r, 0.0);
  
  RNGScope scope;
  
  for (int m = 0; m < M; ++m) {
    std::fill(z.begin(), z.end(), 0.0);
    bool active = true;
    
    for (int j = 0; j < r; ++j) {
      double shift = 0.0;
      for (int k = 0; k < j; ++k) {
        shift += L(j, k) * z[k];
      }
      
      double alpha = (Cvec[j] - mu[j] - shift) / L(j, j);
      double tail_prob = R::pnorm(alpha, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/0);
      
      if (!R_finite(tail_prob) || tail_prob <= 0.0) {
        log_w[m] = R_NegInf;
        active = false;
        break;
      }
      
      log_w[m] += std::log(std::max(tail_prob, eps));
      
      // Equivalent to z = qnorm(Phi(alpha) + U * (1 - Phi(alpha))).
      // We work with the upper-tail probability for numerical stability.
      double u = unif_rand();
      double survival_q = std::max(u * tail_prob, eps);
      z[j] = R::qnorm(survival_q, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/0);
    }
    
    if (!active) continue;
  }
  
  double lme = log_mean_exp(log_w);
  if (!R_finite(lme)) return 0.0;
  return std::exp(lme);
}

// [[Rcpp::export]]
List mc_sov_censored_nd_cpp_details(NumericVector mu,
                                    NumericMatrix Sigma,
                                    NumericVector C,
                                    int M = 10000,
                                    double eps = 1e-300) {
  int r = mu.size();
  if (r < 1) stop("mu must have positive length.");
  if (Sigma.nrow() != r || Sigma.ncol() != r) {
    stop("Sigma dimensions must match length(mu).");
  }
  if (M <= 0) stop("M must be a positive integer.");
  
  std::vector<double> Cvec = threshold_vector(C, r);
  NumericMatrix L = chol_lower_cpp(Sigma);
  std::vector<double> log_w(M, 0.0);
  std::vector<double> z(r, 0.0);
  
  RNGScope scope;
  
  for (int m = 0; m < M; ++m) {
    std::fill(z.begin(), z.end(), 0.0);
    bool active = true;
    
    for (int j = 0; j < r; ++j) {
      double shift = 0.0;
      for (int k = 0; k < j; ++k) {
        shift += L(j, k) * z[k];
      }
      
      double alpha = (Cvec[j] - mu[j] - shift) / L(j, j);
      double tail_prob = R::pnorm(alpha, 0.0, 1.0, 0, 0);
      
      if (!R_finite(tail_prob) || tail_prob <= 0.0) {
        log_w[m] = R_NegInf;
        active = false;
        break;
      }
      
      log_w[m] += std::log(std::max(tail_prob, eps));
      
      double u = unif_rand();
      double survival_q = std::max(u * tail_prob, eps);
      z[j] = R::qnorm(survival_q, 0.0, 1.0, 0, 0);
    }
    
    if (!active) continue;
  }
  
  double lme = log_mean_exp(log_w);
  double prob = R_finite(lme) ? std::exp(lme) : 0.0;
  
  NumericVector weights(M);
  for (int i = 0; i < M; ++i) {
    weights[i] = R_finite(log_w[i]) ? std::exp(log_w[i]) : 0.0;
  }
  
  double mean_w = mean(weights);
  double ss = 0.0;
  for (int i = 0; i < M; ++i) {
    double d = weights[i] - mean_w;
    ss += d * d;
  }
  double se = M > 1 ? std::sqrt(ss / (M - 1)) / std::sqrt(static_cast<double>(M)) : NA_REAL;
  
  return List::create(
    _["prob"] = prob,
    _["se"] = se,
    _["weights"] = weights
  );
}
