#' Compute ME-LS Estimates and Cluster-Robust Standard Errors (iid)
#'
#' This function calculates iid classical measurement error (ME) corrected
#' least squares (LS) estimates along with (up to two-way) cluster-robust
#' standard errors.
#'
#' @param X Matrix of covariates
#' @param Y Vector or matrix of outcomes
#' @param Sigma Covariance matrix for ME (for non-MARS first-steps, which have
#'   no ME, entries will be zero)
#' @param cluster1 Vector specifying primary clustering variable (default NULL)
#' @param cluster2 Optional vector specifying secondary clustering variable for
#'   two-way clustering
#' @param eigfix Logical indicating whether to apply eigenvalue correction
#'   (default FALSE)
#' @param tol Tolerance level for eigenvalue correction (default 1e-8)
#' @param level Level for the confidence interval (default: 0.95)
#'
#' @return A list containing:
#' \itemize{
#'   \item coefficients: Point estimates
#'   \item vcov: Variance-covariance matrix
#'   \item robust_se: Cluster-robust standard errors
#'   \item p_values: Two-sided p-values
#'   \item ci_lower: Lower bounds of confidence intervals
#'   \item ci_upper: Upper bounds of confidence intervals
#' }
#'
#' @export
me_ls_iid <- function(
  X, Y, Sigma, cluster1 = NULL,
  cluster2 = NULL, eigfix = FALSE, tol = 1e-8,
  level = 0.95
) {
  # Get dimensions
  K <- ncol(X)
  n <- nrow(X)

  # Set default cluster1 if NULL
  if (is.null(cluster1)) {
    cluster1 <- seq_len(n)
  }

  # Calculate point estimate
  beta_hat <- solve((t(X) %*% X) - (n * Sigma)) %*% (t(X) %*% Y)
  
  # Calculate bread matrix
  bread <- solve((t(X) %*% X) - (n * Sigma))
  
  # Helper function to create meat matrix
  create_meat <- function(cluster, X, Y, Sigma, beta_hat) {
    meat <- matrix(0, nrow = K, ncol = K)
    clusters <- unique(cluster)
    pb <- txtProgressBar(min = 0, max = length(clusters), initial = 0, style = 3)
    step <- 0
    for (g in clusters) {
      indices <- which(cluster == g)
      X_g <- X[indices, , drop = FALSE]
      Y_g <- Y[indices, , drop = FALSE]
      moment <- matrix(0, nrow = ncol(X_g), ncol = 1)
      for (i in seq_len(nrow(X_g))) {
        moment <- moment +
          t(X_g[i, , drop = FALSE]) %*%
          (Y_g[i, ] - X_g[i, , drop = FALSE] %*% beta_hat) +
          (Sigma %*% beta_hat)
      }
      meat <- meat + (moment %*% t(moment))
      step <- step + 1
      setTxtProgressBar(pb, step)
    }
    close(pb)
    meat
  }
  
  # Calculate variance-covariance matrix
  if (is.null(cluster2)) {
    # One-way clustering
    vcov_matrix <- bread %*% create_meat(cluster1, X, Y, Sigma, beta_hat) %*% bread
  } else {
    # Two-way clustering
    meat1 <- create_meat(cluster1, X, Y, Sigma, beta_hat)
    meat2 <- create_meat(cluster2, X, Y, Sigma, beta_hat)
    cluster12 <- interaction(cluster1, cluster2, drop = TRUE)
    meat12 <- create_meat(cluster12, X, Y, Sigma, beta_hat)
    vcov_matrix <- bread %*% (meat1 + meat2 - meat12) %*% bread
  }
  
  # Apply eigenvalue correction if requested
  if (eigfix) {
    eig <- eigen(vcov_matrix)
    corrected_values <- pmax(eig$values, tol)
    vcov_matrix <- eig$vectors %*% diag(corrected_values) %*% t(eig$vectors)
  }
  
  # Calculate standard errors and confidence intervals
  robust_se <- sqrt(diag(vcov_matrix))
  p_values <- 2 * (1 - pnorm(abs(beta_hat / robust_se)))
  ci_width <- qnorm((1 + level) / 2) * robust_se
  ci_upper <- beta_hat + ci_width
  ci_lower <- beta_hat - ci_width
  
  list(
    coefficients = beta_hat,
    vcov = vcov_matrix,
    robust_se = robust_se,
    p_values = p_values,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}

#' Compute ME-LS Estimates and Cluster-Robust Standard Errors (inid)
#'
#' This function calculates independent (but not identically distributed)
#' classical measurement error (ME) corrected least squares (LS) estimates
#' along with (up to two-way) cluster-robust standard errors.
#'
#' @param X Matrix of covariates
#' @param Y Vector or matrix of outcomes
#' @param Sigma List of covariance matrices for ME, one per observation
#'   (for non-MARS first-steps, which have no ME, entries will be zero)
#' @param cluster1 Vector specifying primary clustering variable (default NULL)
#' @param cluster2 Optional vector specifying secondary clustering variable for
#'   two-way clustering
#' @param eigfix Logical indicating whether to apply eigenvalue correction
#'   (default FALSE)
#' @param tol Tolerance level for eigenvalue correction (default 1e-8)
#' @param level Level for the confidence interval (default: 0.95)
#'
#' @return A list containing:
#' \itemize{
#'   \item coefficients: Point estimates
#'   \item vcov: Variance-covariance matrix
#'   \item robust_se: Cluster-robust standard errors
#'   \item p_values: Two-sided p-values
#'   \item ci_lower: Lower bounds of confidence intervals
#'   \item ci_upper: Upper bounds of confidence intervals
#' }
#'
#' @export
me_ls_inid <- function(
  X, Y, Sigma, cluster1 = NULL, cluster2 = NULL,
  eigfix = FALSE, tol = 1e-8, level = 0.95
) {
  # Get dimensions
  K <- ncol(X)
  n <- nrow(X)
  
  # Set default cluster1 if NULL
  if (is.null(cluster1)) {
    cluster1 <- seq_len(n)
  }
  
  # Calculate average Sigma matrix
  Sigma_bar <- Reduce("+", Sigma) / length(Sigma)

  # Calculate point estimate
  beta_hat <- solve((t(X) %*% X) - (n * Sigma_bar)) %*% (t(X) %*% Y)
  
  # Calculate bread matrix
  bread <- solve((t(X) %*% X) - (n * Sigma_bar))
  
  # Helper function to create meat matrix
  create_meat <- function(cluster, X, Y, Sigma, beta_hat) {
    meat <- matrix(0, nrow = K, ncol = K)
    clusters <- unique(cluster)
    pb <- txtProgressBar(min = 0, max = length(clusters), initial = 0, style = 3)
    step <- 0
    for (g in clusters) {
      indices <- which(cluster == g)
      X_g <- X[indices, , drop = FALSE]
      Y_g <- Y[indices, , drop = FALSE]
      Sigma_g <- Sigma[indices]
      moment <- matrix(0, nrow = ncol(X_g), ncol = 1)
      for (i in seq_len(nrow(X_g))) {
        moment <- moment +
          t(X_g[i, , drop = FALSE]) %*%
          (Y_g[i, ] - X_g[i, , drop = FALSE] %*% beta_hat) +
          (Sigma_g[[i]] %*% beta_hat)
      }
      meat <- meat + (moment %*% t(moment))
      step <- step + 1
      setTxtProgressBar(pb, step)
    }
    close(pb)
    meat
  }
  
  # Calculate variance-covariance matrix
  if (is.null(cluster2)) {
    # One-way clustering
    vcov_matrix <- bread %*% create_meat(cluster1, X, Y, Sigma, beta_hat) %*% bread
  } else {
    # Two-way clustering
    meat1 <- create_meat(cluster1, X, Y, Sigma, beta_hat)
    meat2 <- create_meat(cluster2, X, Y, Sigma, beta_hat)
    cluster12 <- interaction(cluster1, cluster2, drop = TRUE)
    meat12 <- create_meat(cluster12, X, Y, Sigma, beta_hat)
    vcov_matrix <- bread %*% (meat1 + meat2 - meat12) %*% bread
  }
  
  # Apply eigenvalue correction if requested
  if (eigfix) {
    eig <- eigen(vcov_matrix)
    corrected_values <- pmax(eig$values, tol)
    vcov_matrix <- eig$vectors %*% diag(corrected_values) %*% t(eig$vectors)
  }
  
  # Calculate standard errors and confidence intervals
  robust_se <- sqrt(diag(vcov_matrix))
  p_values <- 2 * (1 - pnorm(abs(beta_hat / robust_se)))
  ci_width <- qnorm((1 + level) / 2) * robust_se
  ci_upper <- beta_hat + ci_width
  ci_lower <- beta_hat - ci_width
  
  list(
    coefficients = beta_hat,
    vcov = vcov_matrix,
    robust_se = robust_se,
    p_values = p_values,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}