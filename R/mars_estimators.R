#' MARS Mean Estimator
#'
#' Computes MARS mean estimates and confidence intervals.
#'
#' @param hat_mu Vector of imputed values
#' @param pi Vector of annotation scores
#' @param A Vector of annotation indicators
#' @param M Vector of ground truth values (with non-zero/non-missing entries indicated by A)
#' @param scale Scaling factor (default: 1)
#' @param level Level for the confidence interval (default: 0.95)
#' 
#' @return A list containing:
#' \itemize{
#'   \item est: Point estimate
#'   \item se: Standard error
#'   \item ci_lower: Lower bound of confidence interval
#'   \item ci_upper: Upper bound of confidence interval
#' }
#'
#' @examples
#' # Example usage:
#' # result <- mars_mean_estimator(
#' #   hat_mu = imputed_values_vector,
#' #   pi = annotation_scores_vector,
#' #   A = annotation_indicators_vector,
#' #   M = ground_truth_values_vector,
#' #   scale = 1,
#' #   level = 0.95
#' # )
#'
#' @export
mars_mean_estimator <- function(
  M, A, hat_mu, pi, scale = 1, level = 0.95
) {
  # make sure the vectors are of the same length
  stopifnot(length(hat_mu) == length(pi) && length(hat_mu) == length(A))
  
  # sample size
  I <- length(hat_mu)
  
  # compute the influence function
  pseudo <- hat_mu + (A / pi) * (M - hat_mu)
  
  # estimate (scaled)
  hat_theta <- mean(pseudo) * scale
 
  # compute avar
  psi <- pseudo - hat_theta
  avar <- mean(psi^2) * (scale^2)

  # compute se
  se <- sqrt(avar / I)
  
  # CI width
  ci_width <- qnorm((1 + level) / 2) * se
  
  # output
  list(
    est = hat_theta,
    se   = se,
    ci_lower = hat_theta - ci_width,
    ci_upper = hat_theta + ci_width
  )
}

#' MARS Linear Regression Estimator
#'
#' Computes MARS linear regression coefficient estimates and confidence intervals.
#' 
#' @param C_resid Vector of residualized covariates
#' @param A Vector of annotation indicators
#' @param M Vector of ground truth outcome values (with non-zero/non-missing entries indicated by A)
#' @param hat_mu Vector of imputed values
#' @param pi Vector of annotation scores
#' @param level Level for the confidence interval (default: 0.95)
#'
#' @return A list containing:
#' \itemize{
#'   \item est: Point estimate
#'   \item se: Standard error
#'   \item ci_lower: Lower bound of confidence interval
#'   \item ci_upper: Upper bound of confidence interval
#' }
#' 
#' @examples
#' # Example usage:
#' # result <- mars_linreg_estimator(
#' #   C_resid = residualized_covariates_vector,
#' #   A = annotation_indicators_vector,
#' #   M = ground_truth_outcomes_vector,
#' #   hat_mu = imputed_values_vector,
#' #   pi = annotation_scores_vector,
#' #   level = 0.95
#' # )
#'
#' @export
mars_linreg_estimator <- function(
  C_resid, A, M, hat_mu, pi, level = 0.95
) {
  # make sure all the vectors are of the same length
  stopifnot(length(C_resid) == length(A) && length(C_resid) == length(M) && length(C_resid) == length(hat_mu) && length(C_resid) == length(pi))
  
  # sample size
  I <- length(C_resid)

  # assert that C_resid is mean zero
  stopifnot(mean(C_resid) == 0)

  # compute the pseudo-outcome
  pseudo <- hat_mu + (A / pi) * (M - hat_mu)

  # compute numerator and denominator estimates
  hat_theta_num <- mean(pseudo * C_resid)
  hat_theta_denom <- mean(C_resid^2)

  # compute the estimate
  hat_theta <- hat_theta_num / hat_theta_denom
  
  # compute the avar
  psi <- hat_theta_denom^(-1) * C_resid * (pseudo - C_resid * hat_theta)
  avar <- mean(psi^2)

  # compute se
  se <- sqrt(avar / I)

  # CI width
  ci_width <- qnorm((1 + level) / 2) * se

  # output
  list(
    est = hat_theta,
    se = se,
    ci_lower = hat_theta - ci_width,
    ci_upper = hat_theta + ci_width
  )
}

#' MARS Linear IV Estimator
#'
#' Computes MARS linear IV (constant treatment effect) estimates and confidence intervals.
#' 
#' @param Z_resid Vector of residualized instrument values
#' @param Y Vector of outcome values
#' @param M Vector of ground truth treatment values (with non-zero/non-missing entries indicated by A)
#' @param A Vector of annotation indicators
#' @param hat_mu Vector of imputed values
#' @param pi Vector of annotation scores
#' @param level Level for the confidence interval (default: 0.95)
#' 
#' @return A list containing:
#' \itemize{
#'   \item est: Point estimate
#'   \item se: Standard error
#'   \item ci_lower: Lower bound of confidence interval
#'   \item ci_upper: Upper bound of confidence interval
#' }
#' 
#' @examples
#' # Example usage:
#' # result <- mars_liv_estimator(
#' #   Z_resid = residualized_instrument_vector,
#' #   Y = outcome_values_vector,
#' #   M = ground_truth_outcomes_vector,
#' #   A = annotation_indicators_vector,
#' #   hat_mu = imputed_values_vector,
#' #   pi = annotation_scores_vector,
#' #   level = 0.95
#' # )  
#'
#' @export
mars_liv_estimator <- function(
  Z_resid, Y, M, A, hat_mu, pi, level = 0.95
) {
  # make sure all the vectors are of the same length
  stopifnot(length(Z_resid) == length(Y) && length(Z_resid) == length(M) && length(Z_resid) == length(A) && length(Z_resid) == length(hat_mu) && length(Z_resid) == length(pi))
  
  # sample size
  I <- length(Z_resid)

  # compute the pseudo-outcome
  pseudo <- hat_mu + (A / pi) * (M - hat_mu)
  
  # compute numerator and denominator estimates
  hat_theta_num <- mean(Z_resid * Y)
  hat_theta_denom <- mean(Z_resid * pseudo)

  # compute the estimate
  hat_theta <- hat_theta_num / hat_theta_denom

  # compute the avar
  psi <- hat_theta_denom^(-1) * Z_resid * (Y - pseudo * hat_theta)
  avar <- mean(psi^2)

  # compute se
  se <- sqrt(avar / I)

  # CI width
  ci_width <- qnorm((1 + level) / 2) * se

  # output
  list(
    est = hat_theta,
    se = se,
    ci_lower = hat_theta - ci_width,
    ci_upper = hat_theta + ci_width
  )
}

#' MARS DiD Estimator
#'
#' Computes MARS DiD estimates (Callaway and Sant'Anna, 2021) and confidence intervals.
#' 
#' @param M Vector of ground truth outcome values
#' @param C Vector of indicators for never-treated units
#' @param G Vector of staggered treatment indicators
#' @param A Vector of annotation indicators
#' @param muG Vector of imputed outcomes for (staggered) treated units
#' @param muC Vector of imputed outcomes for never-treated units
#' @param pi Vector of annotation scores
#' @param probG Probability of being in the treated cohort
#' @param probC Probability of being never-treated
#' @param level Level for the confidence interval (default: 0.95)
#' 
#' @return A list containing:
#' \itemize{
#'   \item est: Point estimate
#'   \item se: Standard error
#'   \item ci_lower: Lower bound of confidence interval
#'   \item ci_upper: Upper bound of confidence interval
#' }
#' 
#' @examples
#' # Example usage:
#' # result <- mars_did_estimator(
#' #   M = ground_truth_outcomes_vector,
#' #   C = never_treated_indicator_vector,
#' #   G = staggered_treatment_indicator_vector,
#' #   A = annotation_indicators_vector,
#' #   muG = imputed_outcomes_treated_vector,
#' #   muC = imputed_outcomes_control_vector,
#' #   pi = annotation_scores_vector,
#' #   probG = 0.3,
#' #   probC = 0.7,
#' #   level = 0.95
#' # )
#' 
#' @export
mars_did_estimator <- function(
  M, C, G, A, muG, muC, pi, probG, probC, level = 0.95
) {
  # make sure all the vectors are of the same length
  stopifnot(length(M) == length(C) && length(M) == length(G) && length(M) == length(A) && length(M) == length(muG) && length(M) == length(muC) && length(M) == length(pi))
  
  # validate probability parameters
  stopifnot(probG > 0, probC > 0)

  # sample size
  I <- length(M)

  # compute the pseudo-outcomes
  pseudoC <- muC + (A / pi) * (M - muC)
  pseudoG <- muG + (A / pi) * (M - muG)

  # compute estimates
  hat_theta_G <- mean((G / probG) * pseudoG)
  hat_theta_C <- mean((C / probC) * pseudoC)
  hat_theta <- hat_theta_G - hat_theta_C

  # compute the avar
  psi <- (G / probG) * pseudoG - (C / probC) * pseudoC - hat_theta
  avar <- mean(psi^2)

  # compute se
  se <- sqrt(avar / I)

  # CI width
  ci_width <- qnorm((1 + level) / 2) * se

  # output
  list(
    est = hat_theta,
    se = se,
    ci_lower = hat_theta - ci_width,
    ci_upper = hat_theta + ci_width
  )
}

#' MARS Local RDD Estimator
#'
#' Computes MARS local RDD estimates and confidence intervals.
#' 
#' @param M Vector of ground truth outcome values
#' @param R Vector of running variable values
#' @param D Vector of treatment indicators
#' @param A Vector of annotation indicators
#' @param mu1 Vector of imputed outcomes for treated units
#' @param mu0 Vector of imputed outcomes for control units
#' @param pi Vector of annotation scores
#' @param B Vector of window bounds
#' @param prob1 Probability of being treated and within the window
#' @param prob0 Probability of being control and within the window
#' @param level Level for the confidence interval (default: 0.95)
#' 
#' @return A list containing:
#' \itemize{
#'   \item est: Point estimate
#'   \item se: Standard error
#'   \item ci_lower: Lower bound of confidence interval
#'   \item ci_upper: Upper bound of confidence interval
#' }
#' 
#' @examples
#' # Example usage:
#' # result <- mars_lrdd_estimator(
#' #   M = ground_truth_outcomes_vector,
#' #   R = running_variable_vector,
#' #   D = treatment_indicator_vector,
#' #   A = annotation_indicators_vector,
#' #   mu1 = imputed_outcomes_treated_vector,
#' #   mu0 = imputed_outcomes_control_vector,
#' #   pi = annotation_scores_vector,
#' #   B = c(lower_bound, upper_bound),
#' #   prob1 = 0.5,
#' #   prob0 = 0.5,
#' #   level = 0.95
#' # )
#' 
#' @export
mars_lrdd_estimator <- function(
  M, R, D, A, mu1, mu0, pi, B, prob1, prob0, level = 0.95
) {
  # make sure all the vectors are of the same length
  stopifnot(length(M) == length(R) && length(M) == length(D) && length(M) == length(A) && length(M) == length(mu1) && length(M) == length(mu0) && length(M) == length(pi))

  # validate probability parameters
  stopifnot(prob1 > 0, prob0 > 0)

  # sample size
  I <- length(M)

  # indicator for running variable in window B
  R_in_B <- (R >= B[1]) & (R <= B[2])
  D_is_1 <- ifelse(D == 1, 1, 0)
  D_is_0 <- ifelse(D == 0, 1, 0)

  # compute pseudo-outcomes
  pseudo1 <- mu1 + (A / pi) * (M - mu1)
  pseudo0 <- mu0 + (A / pi) * (M - mu0)

  # compute estimates
  hat_theta_1 <- mean(((R_in_B * D_is_1) / prob1) * pseudo1)
  hat_theta_0 <- mean(((R_in_B * D_is_0) / prob0) * pseudo0)
  hat_theta <- hat_theta_1 - hat_theta_0

  # compute avar
  psi <- ((R_in_B * D_is_1) / prob1) * pseudo1 - ((R_in_B * D_is_0) / prob0) * pseudo0 - hat_theta
  avar <- mean(psi^2)

  # compute se
  se <- sqrt(avar / I)

  # CI width
  ci_width <- qnorm((1 + level) / 2) * se

  # output
  list(
    est = hat_theta,
    se = se,
    ci_lower = hat_theta - ci_width,
    ci_upper = hat_theta + ci_width
  )
}