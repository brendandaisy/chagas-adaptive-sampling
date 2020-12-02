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

// transformed parameters {
//   matrix<lower=0, upper=1>[N, N] P;
//   P = a + b * exp(-c * A);
// }

model {
  a ~ cauchy(0, 2.5);
  b ~ cauchy(0, 2.5);
  c ~ cauchy(0, 2.5);
  tau ~ inv_gamma(.1, .1);
  {
    matrix[N, N] D;
    
    D = rep_matrix(0, N, N);
    for (i in 1:N) {
      for (j in 1:)
      target += bernoulli_logit_lpmf(C[i] | a + b * exp(-c * A[i, 1]))
      C[i, i] = 0;
      D[i, i] = sum(C[i]);
    }
    
  phi ~ multi_normal_prec(rep_row_vector(0, N), tau * (D - C));
  }
  beta ~ cauchy(0, 2.5);
  Y ~ bernoulli_logit(X * beta + phi);
}
