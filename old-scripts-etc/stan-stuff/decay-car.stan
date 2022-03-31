functions {
  matrix lcar(real x0, real k, real alpha, real tau, matrix A, int N) {
    matrix [N, N] C; // transformed A
    matrix [N, N] D; // degree of C
    D = rep_matrix(0, N, N);
    C = 1 - inv(1 + exp(-k * (A - x0)));
    D = rep_matrix(0, N, N);
    for (i in 1:N) {
      C[i, i] = 0;
      D[i, i] = sum(C[i]);
    }
    return tau * (D - alpha * C);
  }
  
  matrix ecar(real x0, real k, real alpha, real tau, matrix A, int N) {
    matrix [N, N] C; // transformed A
    matrix [N, N] D; // degree of C
    D = rep_matrix(0, N, N);
    for (i in 1:N) {
      for (j in i:N) {
	if (i == j)
	  C[i, j] = 0;
	else
	  C[i, j] = A[i, j] <= (x0 + 1) ? 1 : (A[i, j] - x0)^(-k);
	C[j, i] = C[i, j];
      }
      D[i, i] = sum(C[i]);
    }
    return tau * (D - alpha * C);
  }
}

data {
  int<lower = 0> N; // number of non-isolated houses in the village
  int<lower = 1> K; //number of cols in the design matrix
  int<lower = 0, upper = 1> Y[N]; //observations for each house
  matrix[N, K] X; // dummy variable covariates per house, plus intercept
  matrix<lower = 0>[N, N] A; // adjacency matrix

  int dc;
  real k;
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
  vector[K] beta; // slopes/risks of each cov in GLM
}

model {
  if (dc == 1)
    phi ~ multi_normal_prec(zeros, lcar(x0, k, alpha, tau, A, N));
  else if (dc == 2)
    phi ~ multi_normal_prec(zeros, ecar(x0, k, alpha, tau, A, N));
  tau ~ gamma(1, .00005);
  beta ~ cauchy(0, 2.5);
  Y ~ bernoulli_logit(X * beta + phi);
}

generated quantities {
  vector[N] Z_hat; // est. risk without random effect
  int Y_hat[N]; // est. for replicated Y

  Z_hat = X * beta;
  Y_hat = bernoulli_rng(inv_logit(Z_hat + phi));
}
