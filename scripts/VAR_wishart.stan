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
  // matrix[K,K] prior_Rho_loc; // locations for priors on partial correlations
  // matrix[K,K] prior_Rho_scale; // scales for priors on partial correlations
}
////////////////////////////////////////////////////////////////////////////////
transformed data{
  matrix[K,K] I = diag_matrix(rep_vector(1, K));
}
////////////////////////////////////////////////////////////////////////////////
parameters {
  // Temporal
  matrix[K,K] Beta_raw; //
  //real mu_Beta;
  //real<lower=0> sigma_Beta;
  
  // Contemporaneous
  cov_matrix[K] Theta;
}
////////////////////////////////////////////////////////////////////////////////
transformed parameters{
  // Non-centered parameterization for Beta matrix
  matrix[K,K] Beta = Beta_raw .* prior_Beta_scale + prior_Beta_loc;
  //matrix[K,K] Beta = Beta_raw * sigma_Beta + mu_Beta;

  matrix[K,K] Sigma = inverse_spd(Theta);
  
  // Partial correlation matrix
  matrix[K,K] Rho;
  {
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
  target+=   inv_wishart_lpdf(Theta | 3 + K - 1, I);  // prior on precision matrix
  {
    for(t in 2:T){
      // BS: What about intercept?
      vector[K] mu = Beta * Y[t-1,];
       target += multi_normal_lpdf(Y[t,] | mu, Sigma);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
generated quantities{
  vector[T-1] log_lik;
  {
    for(t in 2:T){
      // BS: What about intercept?
      vector[K] mu = Beta * Y[t-1,];
      log_lik[t-1] = multi_normal_lpdf(Y[t, ] | mu, Sigma);
    }
  }
}
