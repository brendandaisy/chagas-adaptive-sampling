data {
  int<lower=0> E; // number of edges in the village
  int<lower = 0> K; // number of cols in the design matrix
  int<lower=0, upper=1> Y[E]; // observations for each edge's target
  matrix[E, K] X; // dummy variable covariates per edge's source
}

parameters {
  vector[K] beta; // slopes for each covariate in GLM
}

model {
  beta ~ cauchy(0, 2.5);
  Y ~ bernoulli_logit(X * beta);
}

/* model { */
/*   target += cauchy_lpdf(beta | 0, 2.5); */
/*   target += bernoulli_lpmf(Y | inv_logit(X * beta)); */
/* } */
