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
# Function to fit the gVAR Model in Stan
#------------------------------------------------------------------------------>
fit_gVAR_stan <-
  function(data,
           priors = NULL,
           backend = "rstan",
           method = "sampling",
           cov_prior = "LKJ", # c("LKJ", "IW")
           iter_sampling = 500,
           iter_warmup = 500,
           n_chains = 4,
           ...) {
    
    
    Y <- data %>% apply(., 2, scale)
    K <- ncol(data)
    n_t <- nrow(data)
    
    # Specify Priors
    if (is.null(priors)){
      prior_Rho_loc <- matrix(.5, nrow = K, ncol = K)
      prior_Rho_scale <- matrix(.4, nrow = K, ncol = K)
      prior_Beta_loc <- matrix(0, nrow = K, ncol = K)
      prior_Beta_scale <- matrix(.5, nrow = K, ncol = K)
    }else{
      prior_Rho_loc <- priors[["prior_Rho_loc"]]
      prior_Rho_scale <- priors[["prior_Rho_scale"]]
      prior_Beta_loc <- priors[["prior_Beta_loc"]]
      prior_Beta_scale <- priors[["prior_Beta_scale"]]
    }
    
    # Stan Data
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
    if (cov_prior == "LKJ") {
      model_name <- "VAR_lkj"
    }
    if (cov_prior == "IW") {
      model_name <- "VAR_wishart"
    }
    
    
    if (backend == "rstan") {
      # Compile model
      stan_model <-
        rstan::stan_model(file = here("scripts", paste0(model_name, ".stan")))
      
      if (method == "sampling") {
        # Run sampler
        stan_fit <- rstan::sampling(
          object = stan_model,
          data = stan_data,
          #pars = c("Beta_raw"),
          #include = FALSE,
          chains = n_chains,
          cores = n_chains,
          iter = iter_sampling + iter_warmup,
          warmup = iter_warmup,
          refresh = 500,
          thin = 1,
          init = .1,
          control = list(adapt_delta = .8),
          ...
        )
      }
      if (method == "variational") {
        stan_fit <- rstan::vb(
          object = stan_model,
          data = stan_data,
          #pars = c("Beta_raw"),
          #include = FALSE,
          init = .1,
          tol_rel_obj = .001,
          output_samples = iter_sampling * n_chains,
          ...
        )
      }
    } else{
      # Compile model
      stan_model <-
        cmdstanr::cmdstan_model(stan_file = here("scripts", paste0(model_name, ".stan")),
                                pedantic = TRUE)
      if (method == "sampling") {
        # Run sampler
        stan_fit <- stan_model$sample(
          data = stan_data,
          chains = n_chains,
          parallel_chains = n_chains,
          iter_warmup = iter_warmup,
          iter_sampling = iter_sampling,
          refresh = 500,
          thin = 1,
          adapt_delta = .8,
          init = .1,
          ...
        )
      }
      if (method == "variational") {
        stan_fit <- stan_model$variational(
          data = stan_data,
          tol_rel_obj = .001,
          init = .1,
          output_samples = iter_sampling * n_chains,
          ...
        )
      }
    }
    return(stan_fit)
  }

#------------------------------------------------------------------------------>
# Function to Compute the LOO-CV for Independent Stan Model and Data
#------------------------------------------------------------------------------>
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


   
   

# -------------------------------------------------------------------------
# Mplus -------------------------------------------------------------------
# -------------------------------------------------------------------------

# Create VAR Syntax in Mplus ----------------------------------------------
mplus_var_syntax <- function(data, 
                             varnames = colnames(data), 
                             laggedvars = NULL, 
                             cores = 1, 
                             iterations = 4000, 
                             model = NULL, 
                             datafile = "mplus_data.dat",
                             syntaxfile = "mplus_model.inp",
                             samplesfile = "mplus_samples.dat") {
  # browser()
  
  # Default variable names (V1, V2, ...)
  if (is.null(varnames)) {
    varnames <- paste0("V", 1:ncol(data))
  }
  
  # Default lagged variables (V1(1), V2(1), ...)
  if (is.null(laggedvars)) {
    laggedvars <- paste0(varnames, "(1)")
  }
  
  # Default model if not specified
  if (is.null(model)) {
    model <- paste(paste(varnames, " ON ", paste(varnames, "&1", sep = "", collapse = " "), ";", sep = ""), collapse = "\n")
  }
  
  # # Write data to Mplus input format using MplusAutomation
  suppressMessages(MplusAutomation::prepareMplusData(data,
                                                     filename = datafile,
                                                     quiet = TRUE))
  
  # Fill in the template
  mplus_template <- "
DATA:	    FILE = @@DATAFILE@@;    ! data file (should be in same folder)

VARIABLE:	NAMES = @@VARNAMES@@;          ! providing names to the variables 
          USEVARIABLES = @@USEVARIABLES@@;   ! select variables for the analysis
	        LAGGED = @@LAGGEDVARS@@;  ! creating first-order
                                    ! lagged observed variables                                    
            MISSING = *;            ! missing value code

ANALYSIS:	ESTIMATOR = BAYES;      ! set estimator (must be Bayes for DSEM) 
	        PROCESSORS = @@CORES@@;         ! using 2 processors
	        BITERATIONS = @@ITERATIONS@@;   ! choose number of iterations;
                                    ! minimum is now 2000; will be more if 
                                    ! the convergence criterion indicates
                                    ! convergence was not reached

MODEL:	    @@MODEL@@

OUTPUT:	    TECH1 TECH8;            ! asking additional output
SAVEDATA: BPARAMETERS IS @@SAMPLESFILE@@; ! saving posterior samples
"
  mplus_syntax <- mplus_template
  
  # Substitute placeholders with actual values
  mplus_syntax <- gsub("@@DATAFILE@@", datafile, mplus_syntax)
  mplus_syntax <- gsub("@@VARNAMES@@", paste(varnames, collapse = " "), mplus_syntax)
  mplus_syntax <- gsub("@@USEVARIABLES@@", paste(varnames, collapse = " "), mplus_syntax)
  mplus_syntax <- gsub("@@LAGGEDVARS@@", paste(laggedvars, collapse = " "), mplus_syntax)
  mplus_syntax <- gsub("@@CORES@@", cores, mplus_syntax)
  mplus_syntax <- gsub("@@ITERATIONS@@", iterations, mplus_syntax)
  mplus_syntax <- gsub("@@MODEL@@", model, mplus_syntax)
  mplus_syntax <- gsub("@@SAMPLESFILE@@", samplesfile, mplus_syntax)
  
  # Write syntax to file
  invisible(writeLines(mplus_syntax, syntaxfile))
  
  return(list(syntax = mplus_syntax))
}


# Convert Mplus samples to array ------------------------------------------
convert_mplus_samples <- function(mplus_output_file) {
  # browser()
  # Extract posterior samples
  mplus_res <- MplusAutomation::readModels(mplus_output_file)
  mplus_samples <- as.data.frame(do.call(rbind, mplus_res$bparameters$valid_draw))
  
  # Extract columns
  samples_beta <- mplus_samples %>% 
    dplyr::select(contains(".ON"))
  samples_sigma <- mplus_samples %>% 
    dplyr::select(contains("with"))
  selected_columns <- names(mplus_samples) %>% 
    str_detect("\\d+_V\\d$")
  samples_diag <- mplus_samples %>% 
    dplyr::select(which(selected_columns))
  nvar <- sum(selected_columns)  # Number of variables
  
  # Convert posterior samples to matrix format
  split_beta <- split(samples_beta, 1:nrow(samples_beta))
  beta_samples <- lapply(split_beta, function(x) {
    matrix(x, nrow = nvar, ncol = nvar, byrow = TRUE)
  })
  # convert this list to an array
  beta_samples <- array(unlist(beta_samples), 
                        dim = c(nvar, nvar, length(beta_samples)))
  split_sigma <- split(samples_sigma, 1:nrow(samples_sigma))
  
  
  # Function to extract covariance values and create covariance matrix
  create_cov_mat <- function(dat, diag_values) {
    cov_matrix <- matrix(NA, nrow = nvar, ncol = nvar)
    diag(cov_matrix) <- 0
    for (i in 1:ncol(dat)) {
      col_name <- names(dat)[i]
      values <- as.vector(strsplit(sub(".*_", "", col_name), "\\."))
      
      # Extract row and column indices from the column name
      row_index <- as.numeric(gsub("V", "", values[[1]][1]))
      col_index <- as.numeric(gsub("V", "", values[[1]][3]))
      
      # Assign the value to the corresponding position in the covariance matrix
      cov_matrix[row_index, col_index] <- as.numeric(dat[, i])
      cov_matrix[col_index, row_index] <- as.numeric(dat[, i])  # Covariance matrix is symmetric
    }
    diag(cov_matrix) <- as.numeric(diag_values)
    return(cov_matrix)
  }
  
  # this is a very ugly mix of a loop and lapply, sorry...
  sigma_samples <- simplify2array(lapply(seq_along(split_sigma), function(i) {
    create_cov_mat(split_sigma[[i]], samples_diag[i,])
  }))
  
  
  # Convert to partial correlations
  pcor_samples <- array(NA, dim = dim(sigma_samples))
  for(i in 1:dim(sigma_samples)[3]) {
    # take inverse of covariance matrix (i.e. precision)
    # and convert to partial correlation matrix
    tmp <- -stats::cov2cor(solve(sigma_samples[,,i]))
    diag(tmp) <- 0
    pcor_samples[, , i] <- tmp
  }
  
  return(list(beta = beta_samples, sigma = sigma_samples, rho = pcor_samples))
}

   