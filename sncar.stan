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
  real psi;
  matrix<lower=0, upper=N>[N, N] C;
  matrix<lower=0, upper=N>[N, N] D;
  D = rep_matrix(0, N, N);
  psi = -log(alpha) / log(1 + gap);
  for (i in 1:N-1) {
    for (j in i+1:N) { // work on upper tri
      if (A[i, j] <= 1) C[i, j] = 1;
      else if (A[i, j] < gap + 1) C[i, j] = pow(A[i, j], -psi);
      else C[i, j] = 0;
      C[j, i] = C[i, j];
    }
    C[i, i] = 0;
    D[i, i] = sum(C[i]);
  }
}

model {
  delta ~ inv_gamma(.1, .1);
  rho ~ beta(2, 1);
  phi ~ multi_normal_prec(rep_row_vector(0, N), delta * (D - rho * C));
  beta ~ cauchy(0, 2.5);
  Y ~ bernoulli_logit(X * beta + phi);
}

