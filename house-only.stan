/* 
TODO: 
make choice of beta priors passable as a function/argument
*/

data {
  int<lower = 0> N; //number of houses in the village
  int<lower = 0> K; //number of cols in the design matrix
  int<lower = 0, upper = 1> Y[N]; //observations for each house
  matrix[N, K] X; // dummy variable covariates per house, plus intercept
}

parameters {
  vector[K] beta; // slopes/risks of each cov in GLM
}

model {
  /* target += cauchy_lpdf(beta | 0, 2.5); */
  /* target += bernoulli_lpmf(Y | inv_logit(X * beta)); */
  beta ~ cauchy(0, 2.5);
  Y ~ bernoulli_logit(X * beta);
}

generated quantities {
  vector[N] r = inv_logit(X * beta); // prob each house is infested
  int y_pred[N] = bernoulli_rng(r);
}

