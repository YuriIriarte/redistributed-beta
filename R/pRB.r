#' Distribution Function of the Redistributed-Beta Distribution
#'
#' Computes the cumulative distribution function (CDF) of the
#' Redistributed-Beta (RB) distribution.
#'
#' The RB distribution is defined through the stochastic representation
#' \deqn{T = X / U,}
#' where \eqn{X \sim Beta(\alpha,\beta)}, \eqn{U \sim Uniform(1,2)},
#' and \eqn{X} and \eqn{U} are independent.
#'
#' @param x Numeric vector of quantiles.
#' @param alpha Positive shape parameter.
#' @param beta Positive shape parameter.
#' @param lower.tail Logical. If \code{TRUE}, probabilities are \eqn{P(T \le x)};
#' otherwise, \eqn{P(T > x)}.
#' @param log.p Logical. If \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return A numeric vector of probabilities.
#'
#' @details
#' The cumulative distribution function is computed using a closed-form
#' expression involving Beta distribution functions. The implementation
#' is numerically stable and handles boundary values appropriately.
#'
#' @examples
#' pRB(0.5, alpha = 2, beta = 5)
#' pRB(c(0.25, 0.5, 0.75), alpha = 2, beta = 5)
#' pRB(0.5, alpha = 2, beta = 5, lower.tail = FALSE)
#'
#' @export
pRB <- function(x, alpha, beta,
                lower.tail = TRUE,
                log.p = FALSE) {
  
  if (length(alpha) != 1L || length(beta) != 1L) {
    stop("alpha and beta must be scalars.", call. = FALSE)
  }
  
  if (!is.numeric(alpha) || !is.numeric(beta) ||
      !is.finite(alpha) || !is.finite(beta) ||
      alpha <= 0 || beta <= 0) {
    stop("alpha and beta must be positive finite scalars.", call. = FALSE)
  }
  
  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }
  
  out <- rep(NA_real_, length(x))
  
  out[is.finite(x) & x <= 0] <- 0
  out[is.finite(x) & x >= 1] <- 1
  
  rho <- exp(lbeta(alpha + 1, beta) - lbeta(alpha, beta))
  
  idx1 <- is.finite(x) & x > 0 & x <= 0.5
  
  if (any(idx1)) {
    
    xx <- x[idx1]
    
    F1_2x <- stats::pbeta(2 * xx, shape1 = alpha,     shape2 = beta)
    F1_x  <- stats::pbeta(xx,     shape1 = alpha,     shape2 = beta)
    
    F2_2x <- stats::pbeta(2 * xx, shape1 = alpha + 1, shape2 = beta)
    F2_x  <- stats::pbeta(xx,     shape1 = alpha + 1, shape2 = beta)
    
    val <- 2 * F1_2x - F1_x -
      rho * xx^(-1) * (F2_2x - F2_x)
    
    out[idx1] <- val
  }
  
  idx2 <- is.finite(x) & x > 0.5 & x < 1
  
  if (any(idx2)) {
    
    xx <- x[idx2]
    
    F1_x <- stats::pbeta(xx, shape1 = alpha,     shape2 = beta)
    F2_x <- stats::pbeta(xx, shape1 = alpha + 1, shape2 = beta)
    
    val <- 2 - F1_x -
      rho * xx^(-1) * (1 - F2_x)
    
    out[idx2] <- val
  }
  
  ok <- is.finite(out)
  out[ok] <- pmin(pmax(out[ok], 0), 1)
  
  if (!lower.tail) {
    out[ok] <- 1 - out[ok]
  }
  
  if (log.p) {
    log_out <- rep(NA_real_, length(out))
    log_out[ok & out == 0] <- -Inf
    log_out[ok & out > 0] <- log(out[ok & out > 0])
    return(log_out)
  }
  
  out
}