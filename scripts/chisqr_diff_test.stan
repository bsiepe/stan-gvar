////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> P; // number of parameters
  int<lower=0> N; // number of iterations
  matrix[N,P] Y_1; // responses
  matrix[N,P] Y_2; // responses
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
parameters {
  // Temporal
  matrix<lower=0>[P,2] mu;
  vector<lower=0>[2] mu_mu;
  vector<lower=0>[2] sigma_mu;
  
  matrix<lower=0>[P,2] sigma;
  vector<lower=0>[2] mu_sigma;
  vector<lower=0>[2] sigma_sigma;
  
  // BS: difference also needs to be specified?
  real diff;
  
  
}
////////////////////////////////////////////////////////////////////////////////
model {
  // Priors
  // BS: Shouldn't we loop here? It's an array above
  // BS: Not sure if we also need to index the mu_sigma and sigma_sigma
  for (p in 1:P) {
    mu[p] ~ normal(mu_mu, sigma_mu);
    sigma[p] ~ normal(mu_sigma, sigma_sigma);
  }  
  mu_mu       ~ student_t(3,0,1);
  sigma_mu    ~ student_t(3,0,1);
  
  diff        ~ std_normal();
  
  sigma[,1]   ~ normal(mu_sigma[1], sigma_sigma[1]);
  sigma[,2]   ~ normal(mu_sigma[1], sigma_sigma[1]); //BS: Should be sigma_sigma[2]?
  mu_sigma    ~ student_t(3,0,1);
  sigma_sigma ~ student_t(3,0,1);
   
  for(p in 1:P){
    Y_1[,p]  ~ normal(mu[p,1], sigma[p,1]);
    Y_2[,p] ~ normal(mu[p,2] + diff, sigma[p,2]);
  }

}
////////////////////////////////////////////////////////////////////////////////
generated quantities{
  real diff_sqr = diff^2;
}
