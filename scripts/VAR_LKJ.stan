////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> K; // number of predictors
  int<lower=0> T; // number of time points
  //int<lower=0> N; // number of 
  array[T] vector[K] Y; // responses
  // Priors
  matrix[K,K] prior_Beta_loc; // locations for priors on Beta matrix
  matrix[K,K] prior_Beta_scale; // scales for priors on Beta matrix
  matrix[K,K] prior_Rho_loc; // locations for priors on partial correlations
  matrix[K,K] prior_Rho_scale; // scales for priors on partial correlations
}
////////////////////////////////////////////////////////////////////////////////
parameters {
  // Temporal
  matrix[K,K] Beta_raw; //
  //real mu_Beta;
  //real<lower=0> sigma_Beta;
  
  // Contemporaneous
  cholesky_factor_corr[K] L_Theta;
  vector<lower=0>[K] sigma_theta;
}
////////////////////////////////////////////////////////////////////////////////
transformed parameters{
  // Non-centered parameterization for Beta matrix
  matrix[K,K] Beta = Beta_raw .* prior_Beta_scale + prior_Beta_loc;
  //matrix[K,K] Beta = Beta_raw * sigma_Beta + mu_Beta;
  
  // Precision matrix
  // BS: So this is just the Cholesky decomposition right?
  // and the pre-multiplication scales it?
  matrix[K,K] Sigma = diag_pre_multiply(sigma_theta, L_Theta) * 
                      diag_pre_multiply(sigma_theta, L_Theta)'; 
  // Partial correlation matrix
  matrix[K,K] Rho;
  {
    matrix[K,K] Theta = inverse_spd(Sigma); 
    for(i in 1:K){
      for(j in 1:K){
        if(i != j){
          Rho[i,j] = -Theta[i,j] / sqrt(Theta[i,i] * Theta[j,j]);
        }else{
          Rho[i,j] = 0;
        } // end else
      } // end j
    } // end i
  }
}
////////////////////////////////////////////////////////////////////////////////
model {
  // Priors
  target+=   std_normal_lpdf(to_vector(Beta_raw));    // prior on Beta
  //target+= student_t_lpdf(mu_Beta | 3,0,2);
  //target+= student_t_lpdf(sigma_Beta | 3,0,2);
  target+=   lkj_corr_cholesky_lpdf(L_Theta | 1);  // prior on Cholesky factor
  target+=   student_t_lpdf(sigma_theta | 3,0,2);   // prior on sigma_theta
  // Priors on partial correlations
  for(i in 1:K){
    for(j in 1:K){
      if(i < j){
        // BS: I think this maps beta to the unit interval?
        target+= beta_proportion_lpdf(
          Rho[i,j] / 2 + 0.5 | prior_Rho_loc[i,j], prior_Rho_scale[i,j]);
        }
      }
    }
  {
    matrix[K, K] Sigma_chol = diag_pre_multiply(sigma_theta, L_Theta);
    for(t in 2:T){
      // BS: What about intercept?
      vector[K] mu = Beta * Y[t-1,];
      target += multi_normal_cholesky_lpdf(Y[t,] | mu, Sigma_chol);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
generated quantities{
  vector[T-1] log_lik;
  {
    matrix[K, K] Sigma_chol = diag_pre_multiply(sigma_theta, L_Theta);
    for(t in 2:T){
      // BS: What about intercept?
      vector[K] mu = Beta * Y[t-1,];
      log_lik[t-1] = multi_normal_cholesky_lpdf(Y[t, ] | mu, Sigma_chol);
    }
  }
}
