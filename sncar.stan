data {
  /* regression data */
  int<lower = 0> N; // number of non-isolated houses in the village
  int<lower = 0> K; // number of cols in the design matrix
  int<lower = 0, upper = 1> Y[N]; //observations for each house
  matrix[N, K] X; // dummy variable covariates per house, plus intercept
  
  /* neighborhood data */
  matrix<lower = 0>[N, N] A; // distance matrix, scaled appropriately
  real d_max;
  
  real ep;
  real alpha;
}

parameters {
  vector[K] beta; // contrasts of each risk factor
  vector[N] phi; // spatial effects
  real delta;
  real rho;
  real<lower=1 + ep, upper=d_max> gap;
}

transformed parameters {
  
}

model {
  Y ~ bernoulli_logit(X * beta + phi);
}

