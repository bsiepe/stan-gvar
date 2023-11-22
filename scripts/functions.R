#------------------------------------------------------------------------------>
# Functions to extract the log_lik for the gVAR
#------------------------------------------------------------------------------>
# Helper function to transform draws matrix into a list of matrices for Beta or Sigma
draws2list <- function(draws_matrix) {
  iterations_list <-
    lapply(
      X = 1:nrow(draws_matrix),
      FUN = function(X) {
        matrix(draws_matrix[X, ], ncol = sqrt(ncol(draws_matrix)), byrow = FALSE)
      }
    )
  return(iterations_list)
}


# Function to compute the log_lik for the gVAR
log_lik_gVAR <- function(Y, draws_beta, draws_sigma, n_cores = 1) {
  
  # # save chain ids
  # chain_ids <- draws_beta %>%
  #   as_draws_df() %>%
  #   select(.chain) %>% unlist()
  
  # prepare matrices from draws
  Beta <- draws2list(draws_beta)
  Sigma <- draws2list(draws_sigma)
  # number of iterations
  n_iter <- length(Beta)
  n_t <- nrow(Y)
  # parallelization
  if (n_cores > 1) {
    future::plan("multisession", workers = n_cores)
  } else{
    future::plan("sequential")
  }
  future::plan("multisession", workers = n_cores)
  # progress bar
  p <- progressr::progressor(along = 1:n_iter)
  # loop over iteraions
  progressr::with_progress(
    log_lik_list <- furrr::future_map(
      .x = 1:n_iter,
      .f = function(n) {
        # loop over time points
        log_lik_row <-
          lapply(2:n_t, function(t) {
            return(mvtnorm::dmvnorm(
              x = Y[t, ],
              mean = Beta[[n]] %*% Y[t - 1, ],
              sigma = Sigma[[n]],
              log = TRUE
            ))
          }) # end timepoint
        log_lik_row <- do.call(cbind, log_lik_row)
        return(log_lik_row)
        p()
      }
      # end iterations
      # end progress bar
    ))
  # end paralleliztion
  
  if (n_cores > 1) {
    future::plan("sequential")
  }
  log_lik_list <- do.call(rbind, log_lik_list)
  
  log_lik_mat <- posterior::as_draws_matrix(log_lik_list)
  return(log_lik_mat)
}

#------------------------------------------------------------------------------>
fit_gVAR_stan <-
  function(data,
           priors,
           backend = "rstan",
           iter_sampling = 500,
           iter_warmup = 500,
           n_chains = 4) {
    # Specify Priors
    prior_Rho_loc <- priors[["prior_Rho_loc"]]
    prior_Rho_scale <- priors[["prior_Rho_scale"]]
    prior_Beta_loc <- priors[["prior_Beta_loc"]]
    prior_Beta_scale <- priors[["prior_Beta_scale"]]
    
    Y <- data %>% apply(., 2, scale)
    K <- ncol(data)
    n_t <- nrow(data)
    
    stan_data <- list(
      K = K,
      "T" = n_t,
      Y = as.matrix(Y),
      prior_Rho_loc = prior_Rho_loc,
      prior_Rho_scale = prior_Rho_scale,
      prior_Beta_loc = prior_Beta_loc,
      prior_Beta_scale = prior_Beta_scale
    )
    # Choose model to fit
    model_name <- "VAR_lkj"
    
    if (backend == "rstan") {
      # Compile model
      stan_model <-
        rstan::stan_model(file = here("scripts", paste0(model_name, ".stan")))
      # Run sampler
      stan_fit <- rstan::sampling(
        object = stan_model,
        data = stan_data,
        pars = c("Beta_raw"),
        include = FALSE,
        seed = 2023,
        chains = n_chains,
        cores = n_chains,
        iter = iter_sampling + iter_warmup,
        warmup = iter_warmup,
        refresh = 500,
        thin = 1,
        init = .1,
        control = list(adapt_delta = .8)
      )
    } else{
      # Compile model
      stan_model <-
        cmdstanr::cmdstan_model(stan_file = here("scripts", paste0(model_name, ".stan")),
                                pedantic = TRUE)
      # Run sampler
      stan_fit <- stan_model$sample(
        data = stan_data,
        seed = 35032,
        chains = n_chains,
        parallel_chains = n_chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        refresh = 500,
        thin = 1,
        adapt_delta = .8,
        init = .1
      )
    }
    return(stan_fit)
  }

loo_gVAR <- function(stan_fit, data, n_cores = 1) {
  c <- class(stan_fit)
  if (attr(c, "package") == "rstan") {
    log_lik <-
      log_lik_gVAR(
        Y = data %>% apply(., 2, scale),
        draws_beta = posterior::as_draws_matrix(rstan::extract(
          stan_fit, pars = "Beta", permuted = FALSE
        )),
        draws_sigma = posterior::as_draws_matrix(rstan::extract(
          stan_fit, pars = "Sigma", permuted = FALSE
        )),
        n_cores = n_cores
      )
    chain_ids <-
      rstan::extract(stan_fit, pars = "Beta", permuted = FALSE) %>%
      posterior::as_draws_df() %>%
      dplyr::select(.chain) %>%
      unlist()
  } else{
    log_lik <-
      log_lik_gVAR(
        Y = data %>% apply(., 2, scale),
        draws_beta = posterior::as_draws_matrix(stan_fit$draws("Beta")),
        draws_sigma = posterior::as_draws_matrix(stan_fit$draws("Sigma")),
        n_cores = n_cores
      )
    chain_ids <- stan_fit$draws("Beta") %>%
      posterior::as_draws_df() %>%
      dplyr::select(.chain) %>%
      unlist()
  }
  loo <- loo::loo(log_lik, r_eff = relative_eff(log_lik, chain_ids))
  return(loo)
}




# -------------------------------------------------------------------------
# Helper functions for model evaluation -----------------------------------
# -------------------------------------------------------------------------
# Convert Stan fit to array -----------------------------------------------
# TODO should maybe be flexible to incorporate something else besides Sigma?
# i.e. also theta (precision matrix?)

stan_fit_convert <- 
  function(stan_fit){
   # check fitting backend
   c <- class(stan_fit)
    
   if (attr(c, "package") == "rstan") {  
    draws_beta <- posterior::as_draws_matrix(rstan::extract(
       stan_fit, pars = "Beta", permuted = FALSE
     ))
    draws_sigma <- posterior::as_draws_matrix(rstan::extract(
       stan_fit, pars = "Sigma", permuted = FALSE
     ))
    draws_rho <- posterior::as_draws_matrix(rstan::extract(
       stan_fit, pars = "Rho", permuted = FALSE
     ))
  } 
   else{
    draws_beta <- posterior::as_draws_matrix(stan_fit$draws("Beta"))
    draws_sigma <- posterior::as_draws_matrix(stan_fit$draws("Sigma"))
    draws_rho <- posterior::as_draws_matrix(stan_fit$draws("Rho"))
  }
  # Convert to array of p x p matrices
  nvar <- sqrt(ncol(draws_beta)) 
   
  # Beta 
  split_beta <- split(draws_beta, seq(nrow(draws_beta)))
  beta_l <- lapply(split_beta, function(x) {
    matrix(x, nrow = nvar, ncol = nvar, byrow = TRUE)
  })
  beta_array <- array(unlist(beta_l), dim = c(nvar, nvar, nrow(draws_beta)))
  
  # Sigma
  split_sigma <- split(draws_sigma, seq(nrow(draws_sigma)))
  sigma_l <- lapply(split_sigma, function(x) {
    matrix(x, nrow = nvar, ncol = nvar, byrow = TRUE)
  })
  sigma_array <- array(unlist(sigma_l), dim = c(nvar, nvar, nrow(draws_sigma)))
  
  # Rho
  split_rho <- split(draws_rho, seq(nrow(draws_rho)))
  rho_l <- lapply(split_rho, function(x) {
    matrix(x, nrow = nvar, ncol = nvar, byrow = TRUE)
  })
  rho_array <- array(unlist(rho_l), dim = c(nvar, nvar, nrow(draws_rho)))
  
  # Return
  return(list(beta = beta_array, sigma = sigma_array, rho = rho_array))
  
}



# Compare fit to DGP ------------------------------------------------------
array_compare_dgp <- function(post_samples, 
                              dgp = NULL,
                              plot = TRUE,
                              dgp_pcor_name = "pcor") {

  # Compute mean for each array element across the third dimension 
  # of post_samples
  post_samples_mean <- lapply(post_samples, function(x) {
    apply(x, c(1, 2), mean)
  })
  post_samples_median <- lapply(post_samples, function(x) {
    apply(x, c(1, 2), median)
  })
  
  # Compare median of beta to DGP
  beta_diff <- post_samples_median$beta - dgp$beta
  rho_diff <- post_samples_median$rho - dgp[[dgp_pcor_name]]
  
  
  result <- list(beta_diff = beta_diff, rho_diff = rho_diff)
  
  if(isTRUE(plot)){
    # Plot both difference matrixes using cor.plot
    par(mfrow = c(1, 2))
    psych::cor.plot(beta_diff, main = "Beta difference")
    psych::cor.plot(rho_diff, main = "Rho difference", upper = FALSE)
  }

  return(result)
}


   
   
   
   