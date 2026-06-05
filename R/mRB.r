#' Raw Moments of the Redistributed-Beta Distribution
#'
#' Computes the raw moment of order \eqn{r} for the
#' Redistributed-Beta (RB) distribution.
#'
#' The r-th raw moment is given by
#'
#' \deqn{
#' \mathbb{E}(T^r)
#' =
#' \frac{B(\alpha + r, \beta)}{B(\alpha, \beta)}
#' \times
#' \tau_r,
#' }
#'
#' where
#'
#' \deqn{
#' \tau_r =
#' \left\{
#' \begin{array}{ll}
#' \log(2), & r = 1, \\
#' \dfrac{1 - 2^{1-r}}{r - 1}, & r \neq 1.
#' \end{array}
#' \right.
#' }
#'
#' @param r Non-negative order or vector of orders of the moment.
#' @param alpha Positive shape parameter.
#' @param beta Positive shape parameter.
#'
#' @return A numeric vector with the raw moments of order \eqn{r}.
#'
#' @examples
#' mRB(1, alpha = 2, beta = 5)
#' mRB(1:4, alpha = 2, beta = 5)
#'
#' @export
mRB <- function(r, alpha, beta) {
  
  if (!is.numeric(r) || any(!is.finite(r)) || any(r < 0)) {
    stop("r must be non-negative and finite.", call. = FALSE)
  }
  
  if (length(alpha) != 1L || length(beta) != 1L) {
    stop("alpha and beta must be scalars.", call. = FALSE)
  }
  
  if (!is.numeric(alpha) || !is.numeric(beta) ||
      !is.finite(alpha) || !is.finite(beta) ||
      alpha <= 0 || beta <= 0) {
    stop("alpha and beta must be positive finite scalars.", call. = FALSE)
  }
  
  beta_ratio <- exp(lbeta(alpha + r, beta) - lbeta(alpha, beta))
  
  tau_r <- ifelse(
    abs(r - 1) < 1e-10,
    log(2),
    (1 - 2^(1 - r)) / (r - 1)
  )
  
  tau_r * beta_ratio
}


#' Mean of the Redistributed-Beta Distribution
#'
#' Computes the theoretical mean of the RB distribution.
#'
#' @inheritParams mRB
#'
#' @return A numeric value.
#'
#' @examples
#' meanRB(alpha = 2, beta = 5)
#'
#' @export
meanRB <- function(alpha, beta) {
  mRB(1, alpha = alpha, beta = beta)
}


#' Variance of the Redistributed-Beta Distribution
#'
#' Computes the theoretical variance of the RB distribution.
#'
#' @inheritParams mRB
#'
#' @return A numeric value.
#'
#' @examples
#' varRB(alpha = 2, beta = 5)
#'
#' @export
varRB <- function(alpha, beta) {
  m1 <- mRB(1, alpha = alpha, beta = beta)
  m2 <- mRB(2, alpha = alpha, beta = beta)
  
  m2 - m1^2
}


#' Coefficient of Variation of the Redistributed-Beta Distribution
#'
#' Computes the theoretical coefficient of variation of the RB distribution.
#'
#' @inheritParams mRB
#'
#' @return A numeric value.
#'
#' @examples
#' cvRB(alpha = 2, beta = 5)
#'
#' @export
cvRB <- function(alpha, beta) {
  mu <- meanRB(alpha = alpha, beta = beta)
  sig2 <- varRB(alpha = alpha, beta = beta)
  
  sqrt(sig2) / mu
}


#' Skewness of the Redistributed-Beta Distribution
#'
#' Computes the theoretical skewness coefficient of the RB distribution.
#'
#' @inheritParams mRB
#'
#' @return A numeric value.
#'
#' @examples
#' skewRB(alpha = 2, beta = 5)
#'
#' @export
skewRB <- function(alpha, beta) {
  m1 <- mRB(1, alpha = alpha, beta = beta)
  m2 <- mRB(2, alpha = alpha, beta = beta)
  m3 <- mRB(3, alpha = alpha, beta = beta)
  
  mu3 <- m3 - 3 * m1 * m2 + 2 * m1^3
  sig2 <- m2 - m1^2
  
  mu3 / sig2^(3 / 2)
}


#' Kurtosis of the Redistributed-Beta Distribution
#'
#' Computes the theoretical kurtosis or excess kurtosis of the RB distribution.
#'
#' @inheritParams mRB
#' @param excess Logical. If \code{TRUE}, the excess kurtosis is returned;
#' otherwise, the ordinary kurtosis is returned.
#'
#' @return A numeric value.
#'
#' @examples
#' kurtRB(alpha = 2, beta = 5)
#' kurtRB(alpha = 2, beta = 5, excess = FALSE)
#'
#' @export
kurtRB <- function(alpha, beta, excess = TRUE) {
  m1 <- mRB(1, alpha = alpha, beta = beta)
  m2 <- mRB(2, alpha = alpha, beta = beta)
  m3 <- mRB(3, alpha = alpha, beta = beta)
  m4 <- mRB(4, alpha = alpha, beta = beta)
  
  mu4 <- m4 - 4 * m1 * m3 + 6 * m1^2 * m2 - 3 * m1^4
  sig2 <- m2 - m1^2
  
  kurt <- mu4 / sig2^2
  
  if (excess) {
    kurt <- kurt - 3
  }
  
  kurt
}