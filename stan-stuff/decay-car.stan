data {
  int<lower = 0> N; // number of non-isolated houses in the village
  int<lower = 1> K; //number of cols in the design matrix
  int<lower = 0, upper = 1> Y[N]; //observations for each house
  matrix[N, K] X; // dummy variable covariates per house, plus intercept
  matrix<lower = 0>[N, N] A; // adjacency matrix
  
  real tau_a;
  real tau_b;
  real kmin; // max for shape param k
  real kmax;
  real x0min; // max distance for cutoff param x0
  real x0max;
}

transformed data {
  vector[N] zeros;
  zeros = rep_vector(0, N);
}

parameters {
  vector[N] phi; // spatial effects
  real<lower=0> tau; // precision of spatial effects
  real<lower=0, upper=.99> alpha; // strength of spatial correlation
  real<lower=x0min, upper=x0max> x0;
  real<lower=kmin, upper=kmax> k;
  vector[K] beta; // slopes/risks of each cov in GLM
}

model {
  {
    matrix [N, N] Q; // precision
    matrix [N, N] C; // transformed A
    matrix [N, N] D; // degree of C
    C = 1 - inv(1 + exp(-k * (A - x0)));
    D = rep_matrix(0, N, N);
    for (i in 1:N) {
      C[i, i] = 0;
      D[i, i] = sum(C[i]);
    }
    Q = tau * (D - alpha * C);
    phi ~ multi_normal_prec(zeros, Q);
  }
  tau ~ gamma(tau_a, tau_b);
  beta ~ cauchy(0, 2.5);
  Y ~ bernoulli_logit(X * beta + phi);
}

generated quantities {
  vector[N] Z_hat; // est. risk without random effect
  int Y_hat[N]; // est. for replicated Y

  Z_hat = X * beta;
  Y_hat = bernoulli_rng(inv_logit(Z_hat + phi));
}
