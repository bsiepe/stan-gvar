#------------------------------------------------------------------------------>
# Functions to extract the log_lik for the gVAR
#------------------------------------------------------------------------------>
# Helper function to transform draws matrix into a list of matrices for Beta or Sigma
draws2list <- function(draws_matrix) {
  iterations_list <-
    lapply(
      X = 1:nrow(draws_matrix),
      FUN = function(X) {
        matrix(draws_matrix[X,], ncol = sqrt(ncol(draws_matrix)), byrow = FALSE)
      }
    )
  return(iterations_list)
}


# Function to compute the log_lik for the gVAR
log_lik_gVAR <- function(Y, draws_beta, draws_sigma) {
  # prepare matrices from draws
  Beta <- draws2list(draws_beta)
  Sigma <- draws2list(draws_sigma)
  # number of iterations
  n_iter <- length(Beta)
  n_t <- nrow(Y)
  # init log_lik vector
  log_lik <- matrix(NA, nrow = n_iter, ncol = n_t - 1)
  # loop over iteraions
  for(n in 1:n_iter){
    # loop over time points
    for (t in 2:n_t) {
      log_lik[n, t - 1] <-
        mvtnorm::dmvnorm(
          x = Y[t,],
          mean = Beta[[n]] %*% Y[t - 1,],
          sigma = Sigma[[n]],
          log = TRUE
        )
      }
  }
  log_lik <- posterior::as_draws_matrix(log_lik)
  return(log_lik)
}

#------------------------------------------------------------------------------>