functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  * From https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param C_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param C_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*C*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] C_sparse, vector D_sparse, vector lambda, int n, int C_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_C; // phi' * C
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_C = rep_row_vector(0, n);
      for (i in 1:C_n) {
        phit_C[C_sparse[i, 1]] = phit_C[C_sparse[i, 1]] + phi[C_sparse[i, 2]];
        phit_C[C_sparse[i, 2]] = phit_C[C_sparse[i, 2]] + phi[C_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_C * phi)));
  }
}

data {
  int<lower = 1> N;
  int<lower = 1> K; //number of cols in the design matrix
  int<lower = 0, upper = 1> Y[N]; //observations for each house
  matrix[N, K] X; // dummy variable covariates per house, plus intercept
  matrix<lower = 0, upper = 1>[N, N] C; // adjacency matrix
  int C_n;                // number of edges

  real tau_a;
  real tau_b;
}

transformed data {
  int C_sparse[C_n, 2];   // adjacency pairs
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[N] lambda;       // eigenvalues of invsqrtD * C * invsqrtD
  
  { // generate sparse representation for C
  int counter;
  counter = 1;
  // loop over upper triangular part of C to identify neighbor pairs
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (C[i, j] == 1) {
          C_sparse[counter, 1] = i;
          C_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:N) D_sparse[i] = sum(C[i]);
  {
    vector[N] invsqrtD;  
    for (i in 1:N) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(C, diag_matrix(invsqrtD)));
  }
}

parameters {
  vector[N] phi;
  real<lower = 0> tau;
  real<lower = 0, upper = 1> alpha;
  vector[K] beta; // slopes/risks of each cov in GLM
}

model {
  phi ~ sparse_car(tau, alpha, C_sparse, D_sparse, lambda, N, C_n);
  tau ~ gamma(tau_a, tau_b);
  beta ~ cauchy(0, 2.5);
  Y ~ bernoulli_logit(X * beta + phi);
}

generated quantities {
  vector[N] Z_hat; // est. risk without random effect
  int Y_hat[N]; // est. for replicated Y

  Z_hat = X * beta;
  Y_hat = bernoulli_rng(inv_logit(Z_hat + phi));
}
