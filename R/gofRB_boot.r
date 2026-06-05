#' Bootstrap Goodness-of-Fit Tests for the RB Distribution
#'
#' Performs bootstrap goodness-of-fit tests for the 
#' Redistributed-Beta (RB) distribution.
#'
#' The function fits the RB distribution by maximum likelihood, computes the
#' Anderson-Darling and Cramer-von Mises statistics, and obtains bootstrap
#' p-values using a parametric bootstrap procedure.
#'
#' @param x Numeric vector of observations in \eqn{(0,1)}.
#' @param B Number of bootstrap replicates.
#' @param seed Optional seed for reproducibility.
#' @param eps Small positive constant used to avoid evaluating probabilities
#' exactly at 0 or 1.
#' @param verbose Logical. If \code{TRUE}, progress information is printed.
#' @param suppress_warnings Logical. If \code{TRUE}, warnings generated during
#' fitting and bootstrap refitting are suppressed. The default is \code{TRUE}.
#' @param method Optimization method used by \code{fitRB_mle()}.
#' @param multistart Logical. If \code{TRUE}, multiple starting values are used
#' in the maximum likelihood estimation.
#'
#' @return An object of class \code{"gofRB_boot"}, which is a list containing:
#' \itemize{
#'   \item \code{model}: fitted model name.
#'   \item \code{n}: sample size.
#'   \item \code{B}: requested number of bootstrap replicates.
#'   \item \code{B_valid}: number of successful bootstrap replicates.
#'   \item \code{seed}: random seed used.
#'   \item \code{method}: optimization method used.
#'   \item \code{multistart}: whether multiple starting values were used.
#'   \item \code{par_hat}: maximum likelihood estimates.
#'   \item \code{logLik}: maximized log-likelihood.
#'   \item \code{AIC}: Akaike information criterion.
#'   \item \code{BIC}: Bayesian information criterion.
#'   \item \code{AD}: Anderson-Darling statistic and bootstrap p-value.
#'   \item \code{CvM}: Cramer-von Mises statistic and bootstrap p-value.
#'   \item \code{bootstrap}: bootstrap statistics and convergence indicators.
#'   \item \code{fit}: fitted object returned by \code{fitRB_mle()}.
#'   \item \code{call}: matched function call.
#' }
#'
#' @details
#' For each bootstrap replicate, a sample is generated from the fitted RB
#' distribution, the model is refitted, and the goodness-of-fit statistics are
#' recomputed. Bootstrap p-values are computed using the proportion of bootstrap
#' statistics greater than or equal to the observed statistic, with a standard
#' finite-sample correction.
#'
#' @examples
#' set.seed(123)
#' x <- rRB(n = 200, alpha = 2, beta = 5)
#'
#' gof <- gofRB_boot(x, B = 100, seed = 123)
#' summary(gof)
#'
#' @export
gofRB_boot <- function(x,
                       B = 999,
                       seed = 2026,
                       eps = 1e-10,
                       verbose = TRUE,
                       suppress_warnings = TRUE,
                       method = "L-BFGS-B",
                       multistart = FALSE) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    stop("x must have length >= 2.", call. = FALSE)
  }
  
  if (any(x <= 0 | x >= 1)) {
    stop("All observations must lie in (0,1).", call. = FALSE)
  }
  
  if (!is.numeric(B) || length(B) != 1L || B < 1) {
    stop("B must be a positive integer.", call. = FALSE)
  }
  
  B <- as.integer(B)
  n <- length(x)
  
  # ----------------------------------------------------------
  # Statistics
  # ----------------------------------------------------------
  
  AD_stat <- function(x, par) {
    x <- sort(x)
    n <- length(x)
    
    u <- pRB(x, alpha = par["alpha"], beta = par["beta"])
    u <- pmin(pmax(u, eps), 1 - eps)
    
    -n - mean((2 * seq_len(n) - 1) *
                (log(u) + log(1 - rev(u))))
  }
  
  CvM_stat <- function(x, par) {
    x <- sort(x)
    n <- length(x)
    
    u <- pRB(x, alpha = par["alpha"], beta = par["beta"])
    u <- pmin(pmax(u, eps), 1 - eps)
    
    1 / (12 * n) + sum((u - (2 * seq_len(n) - 1) / (2 * n))^2)
  }
  
  # ----------------------------------------------------------
  # Fit
  # ----------------------------------------------------------
  
  fit_main <- function(z) {
    if (isTRUE(suppress_warnings)) {
      suppressWarnings(
        fitRB_mle(
          z,
          methods = method,
          multistart = multistart
        )
      )
    } else {
      fitRB_mle(
        z,
        methods = method,
        multistart = multistart
      )
    }
  }
  
  fit <- fit_main(x)
  
  if (is.null(fit$par) ||
      any(!is.finite(fit$par)) ||
      (!is.null(fit$convergence) && fit$convergence != 0)) {
    stop("The RB model could not be fitted to the original data.", call. = FALSE)
  }
  
  par_hat <- fit$par
  
  AD_obs  <- AD_stat(x, par_hat)
  CvM_obs <- CvM_stat(x, par_hat)
  
  AD_boot   <- rep(NA_real_, B)
  CvM_boot  <- rep(NA_real_, B)
  conv_boot <- rep(FALSE, B)
  
  set.seed(seed)
  
  # ----------------------------------------------------------
  # Bootstrap
  # ----------------------------------------------------------
  
  for (b in seq_len(B)) {
    
    xb <- rRB(
      n = n,
      alpha = par_hat["alpha"],
      beta  = par_hat["beta"]
    )
    
    fit_b <- try(fit_main(xb), silent = TRUE)
    
    if (!inherits(fit_b, "try-error") &&
        !is.null(fit_b$par) &&
        all(is.finite(fit_b$par)) &&
        (is.null(fit_b$convergence) || fit_b$convergence == 0)) {
      
      par_b <- fit_b$par
      
      AD_boot[b]  <- AD_stat(xb, par_b)
      CvM_boot[b] <- CvM_stat(xb, par_b)
      conv_boot[b] <- TRUE
    }
    
    if (isTRUE(verbose) && (b %% 100 == 0L || b == B)) {
      cat("Bootstrap replication:", b, "of", B, "\n")
    }
  }
  
  AD_boot_valid  <- AD_boot[conv_boot]
  CvM_boot_valid <- CvM_boot[conv_boot]
  
  B_valid <- sum(conv_boot)
  
  if (B_valid == 0L) {
    p_AD  <- NA_real_
    p_CvM <- NA_real_
  } else {
    p_AD  <- (1 + sum(AD_boot_valid  >= AD_obs))  / (B_valid + 1)
    p_CvM <- (1 + sum(CvM_boot_valid >= CvM_obs)) / (B_valid + 1)
  }
  
  out <- list(
    model = "RB",
    n = n,
    B = B,
    B_valid = B_valid,
    seed = seed,
    method = method,
    multistart = multistart,
    par_hat = par_hat,
    logLik = fit$logLik,
    AIC = fit$AIC,
    BIC = fit$BIC,
    AD = list(
      statistic = AD_obs,
      p_value = p_AD
    ),
    CvM = list(
      statistic = CvM_obs,
      p_value = p_CvM
    ),
    bootstrap = data.frame(
      AD = AD_boot,
      CvM = CvM_boot,
      convergence = conv_boot
    ),
    fit = fit,
    call = match.call()
  )
  
  class(out) <- c("gofRB_boot", "list")
  out
}