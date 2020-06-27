data {
  int<lower = 0> N; //number of houses in the village
  int<lower = 0> K; //number of cols in the design matrix (= num covariates x levels per covariate), plus intercept
  int<lower = 0, upper = 1> Y[N]; //observations for each house
  matrix[N, K] X; // dummy variable covariates per house, plus intercept
}

parameters {
  vector[K] beta; // contrasts of each risk factor
}

model {
  beta ~ cauchy(0, 2.5);
  Y ~ bernoulli_logit(X * beta);
}

generated quantities {
  vector[N] r = inv_logit(X * beta);
  int y_pred[N] = bernoulli_rng(r);
}

