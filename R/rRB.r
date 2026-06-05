#' Random Generation from the Redistributed-Beta Distribution
#'
#' Generates random observations from the Redistributed-Beta (RB)
#' distribution.
#'
#' The RB distribution is defined through the stochastic representation
#' \deqn{T = X / U,}
#' where \eqn{X \sim Beta(\alpha,\beta)}, \eqn{U \sim Uniform(1,2)},
#' and \eqn{X} and \eqn{U} are independent.
#'
#' @param n Number of observations.
#' @param alpha Positive shape parameter.
#' @param beta Positive shape parameter.
#'
#' @return A numeric vector of random observations in \eqn{(0,1)}.
#'
#' @examples
#' set.seed(123)
#' x <- rRB(100, alpha = 2, beta = 5)
#' hist(x, probability = TRUE)
#' curve(dRB(x, alpha = 2, beta = 5), add = TRUE)
#'
#' @export
rRB <- function(n, alpha, beta) {
  
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 0) {
    stop("n must be a non-negative integer.", call. = FALSE)
  }
  
  n <- as.integer(n)
  
  if (length(alpha) != 1L || length(beta) != 1L) {
    stop("alpha and beta must be scalars.", call. = FALSE)
  }
  
  if (!is.numeric(alpha) || !is.numeric(beta) ||
      !is.finite(alpha) || !is.finite(beta) ||
      alpha <= 0 || beta <= 0) {
    stop("alpha and beta must be positive finite scalars.", call. = FALSE)
  }
  
  X <- stats::rbeta(n, shape1 = alpha, shape2 = beta)
  U <- stats::runif(n, min = 1, max = 2)
  
  X / U
}