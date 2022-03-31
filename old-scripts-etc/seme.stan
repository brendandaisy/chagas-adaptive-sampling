data {
  /* regression data */
  int<lower = 0> N; // number of non-isolated houses in the village
  int<lower = 0> K; // number of cols in the design matrix
  int<lower = 0, upper = 1> Y[N]; //observations for each house
  matrix[N, K] X; // dummy variable covariates per house, plus intercept
  
  /* neighborhood data */
  matrix<lower = 0>[N, N] A; // distance matrix, scaled appropriately
}

parameters {
  vector[K] beta; // contrasts of each risk factor
  vector[N] phi; // spatial effects
  vector<lower=0, upper=1>[N] C[N];
  
  real tau;
  real a;
  real<lower=0> b;
  real<lower=0> c;
}

model {
  Y ~ bernoulli_logit(X * beta + phi);
}
