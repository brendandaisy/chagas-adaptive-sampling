functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  * From https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}

data {
  int<lower = 0> N; // number of non-isolated houses in the village
  int N_iso; // number of isolated houses
  int<lower = 0> K; // number of cols in the design matrix
  int<lower = 0, upper = 1> Y[N + N_iso]; //observations for each house
  matrix[N + N_iso, K] X; // dummy variable covariates per house, plus intercept

  /* neighborhood data */
  matrix<lower = 0, upper = 1>[N, N] W; // adjacency matrix
  int W_n; // number of adjacent region pairs
}

transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[N] lambda;     // eigenvalues of invsqrtD * W * invsqrtD
  
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:N) D_sparse[i] = sum(W[i]);
  { // find the eigenvalues for Jin et al. simplification
    vector[N] invsqrtD;  
    for (i in 1:N) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}

parameters {
  vector[K] beta; // contrasts of each risk factor
  vector[N] phi; // spatial effects
  real<lower=0> tau; // precision of spatial effects
  real<lower = 0, upper = 1> alpha; // strength of spatial correlation
}

model {
  phi ~ sparse_car(tau, alpha, W_sparse, D_sparse, lambda, N, W_n);
  beta ~ cauchy(0, 2.5);
  tau ~ gamma(2, 2);
  Y[1:N] ~ bernoulli_logit(X[1:N] * beta + phi);
  if (N_iso > 0)
    Y[N+1:(N+N_iso)] ~ bernoulli_logit(X[N+1:(N+N_iso)] * beta);
}

/* generated quantities { */
/*   vector[N] r = inv_logit(X * beta); */
/*   int y_pred[N] = bernoulli_rng(r); */
/* } */

