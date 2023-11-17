////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> K; // number of predictors
  int<lower=0> T; // number of time points
  //int<lower=0> N; // number of 
  array[T] vector[K] Y; // responses
  // Priors
  matrix[K,K] prior_beta_loc; // locations for priors on Beta matrix
  matrix[K,K] prior_beta_scale; // scales for priors on Beta matrix
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
  matrix[K,K] Beta = Beta_raw .* prior_beta_scale + prior_beta_loc;
  //matrix[K,K] Beta = Beta_raw * sigma_Beta + mu_Beta;
  
  // Precision matrix
  matrix[K,K] Theta = inverse_spd(
    diag_pre_multiply(sigma_theta, L_Theta) * 
    diag_pre_multiply(sigma_theta, L_Theta)'
    ); 
    // Partial correlation matrix
    matrix[K,K] Rho;
    for(i in 1:K){
      for(j in 1:K){
        if(i != j){
          Rho[i,j] = -Theta[i,j] / sqrt(Theta[i,i] * Theta[j,j]);
        }else{
          Rho[i,j] = 0;
        }
      }
    }
}
////////////////////////////////////////////////////////////////////////////////
model {
  // Priors
  to_vector(Beta_raw) ~ std_normal();
  // mu_Beta          ~ student_t(3,0,2);
  // sigma_Beta       ~ student_t(3,0,2);
  L_Theta             ~ lkj_corr_cholesky(1);
  sigma_theta         ~ student_t(3,0,2);
  // Priors on partial correlations
  for(i in 1:K){
    for(j in 1:K){
      if(i < j){
        target+= beta_proportion_lpdf(
          Rho[i,j] / 2 + 0.5 | prior_Rho_loc[i,j], prior_Rho_scale[i,j]);
        }
      }
    }
  {
    matrix[K, K] Sigma = diag_pre_multiply(sigma_theta, L_Theta);
    for(t in 2:T){
      // BS: What about intercept?
      vector[K] mu = Beta * Y[t-1,];
      Y[t,] ~ multi_normal_cholesky(mu, Sigma);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
generated quantities{}
