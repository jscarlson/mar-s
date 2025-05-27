#' Residualize with a Linear Projection
#'
#' Computes residuals from a linear projection of the target variable on the feature matrix.
#'
#' @param target_vec Vector of the target variable to be residualized
#' @param featur_mat Matrix of the features for residualization
#'
#' @return Vector of residuals from the regression
#'
#' @keywords internal
residualize <- function(target_vec, featur_mat) {
  # add a constant term to the feature matrix if it is not already present
  if (!any(apply(featur_mat, 2, function(x) all(x == 1)))) {
    featur_mat <- cbind(1, featur_mat)
  }
  # fit the linear model
  fit <- lm(target_vec ~ featur_mat)

  # return the residuals
  as.vector(fit$residuals)
}
