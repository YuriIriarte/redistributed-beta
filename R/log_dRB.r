log1mexp <- function(a) {
  out <- rep(NaN, length(a))
  
  ok <- is.finite(a) & (a <= 0)
  
  if (any(ok)) {
    a_ok <- a[ok]
    idx <- a_ok > -log(2)
    
    tmp <- numeric(length(a_ok))
    tmp[idx]  <- log(-expm1(a_ok[idx]))
    tmp[!idx] <- log1p(-exp(a_ok[!idx]))
    
    out[ok] <- tmp
  }
  
  out[is.finite(a) & abs(a) < 1e-15] <- -Inf
  
  out
}


log_dRB_core <- function(x, alpha, beta) {
  
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
  
  out <- rep(-Inf, length(x))
  inside <- is.finite(x) & x > 0 & x < 1
  
  if (!any(inside)) {
    return(out)
  }
  
  xx <- x[inside]
  
  log_const <- lbeta(alpha + 1, beta) - lbeta(alpha, beta)
  log_kernel <- -2 * log(xx)
  
  log_prob <- rep(-Inf, length(xx))
  
  idx1 <- xx <= 0.5
  idx2 <- !idx1
  
  if (any(idx1)) {
    
    z <- xx[idx1]
    z2 <- 2 * z
    
    logF_2z <- stats::pbeta(
      z2,
      shape1 = alpha + 1,
      shape2 = beta,
      log.p = TRUE
    )
    
    logF_z <- stats::pbeta(
      z,
      shape1 = alpha + 1,
      shape2 = beta,
      log.p = TRUE
    )
    
    delta <- logF_z - logF_2z
    delta[is.finite(delta) & delta > 0 & delta < 1e-12] <- 0
    
    lp <- rep(-Inf, length(delta))
    ok <- is.finite(logF_2z) & is.finite(delta) & delta <= 0
    
    if (any(ok)) {
      lp[ok] <- logF_2z[ok] + log1mexp(delta[ok])
    }
    
    log_prob[idx1] <- lp
  }
  
  if (any(idx2)) {
    
    z <- xx[idx2]
    
    log_prob[idx2] <- stats::pbeta(
      z,
      shape1 = alpha + 1,
      shape2 = beta,
      lower.tail = FALSE,
      log.p = TRUE
    )
  }
  
  ans <- log_const + log_kernel + log_prob
  ans[!is.finite(ans)] <- -Inf
  
  out[inside] <- ans
  out
}