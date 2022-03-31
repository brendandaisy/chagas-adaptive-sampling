data {
  int<lower = 0> N; // number of houses in the village
  int<lower = 0> K; // number of cols in the design matrix, NO intercept
  int<lower = 0, upper = 1> Y[N]; // observations for each house
  matrix<lower = 0>[N, K] X; // neighborhood covariates per house, NO intercept

  real<lower = 0, upper = 1> D[N, N];
}

parameters {
  vector[K] beta; // contrasts of each risk factor
  real alpha;
  real gamma;
}

model {
  matrix[N, N] A;
  matrix[N, K] Z;
  beta ~ cauchy(0, 2.5);
  alpha ~ cauchy(0, 2.5);
  gamma ~ beta(1, 3);
  for (i in 1:N) {
    for (j in 1:N) {
      A[i, j] = D[i, j] < h;
    }
  }
  Z = A' * X;
  Y ~ bernoulli_logit(alpha + Z * beta);
}

/* generated quantities { */
/*   vector[N] r = inv_logit(Z * gamma); */
/*   int y_pred[N] = bernoulli_rng(r); */
/* } */

