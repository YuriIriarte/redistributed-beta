#' Maximum Likelihood Estimation for the Redistributed-Beta Distribution
#'
#' Estimates the parameters of the Redistributed-Beta (RB)
#' distribution by maximum likelihood.
#'
#' The log-likelihood is evaluated using a numerically stable log-density
#' implementation. The optimization is performed on the logarithmic scale
#' to ensure positivity of the parameters. Moment-based estimates are used
#' as initial values whenever available, and additional deterministic starting
#' values can be used when \code{multistart = TRUE}.
#'
#' @param x Numeric vector of observations in \eqn{(0,1)}.
#' @param start Optional named numeric vector with initial values for
#' \code{alpha} and \code{beta}.
#' @param starts Optional list of named numeric vectors with additional starting
#' values.
#' @param methods Character vector with optimization methods passed to
#' \code{\link[stats]{optim}}.
#' @param multistart Logical. If \code{TRUE}, a deterministic grid of starting
#' values is used.
#' @param lower_alpha,upper_alpha Lower and upper bounds for the parameter
#' \code{alpha}.
#' @param lower_beta,upper_beta Lower and upper bounds for the parameter
#' \code{beta}.
#' @param control List of control parameters passed to \code{\link[stats]{optim}}.
#' @param show_warnings Logical. If \code{TRUE}, optimization warnings are displayed.
#'
#' @return An object of class \code{"fitRB_mle"}, which is a list containing:
#' \itemize{
#'   \item \code{par}: maximum likelihood estimates.
#'   \item \code{logLik}: maximized log-likelihood.
#'   \item \code{se}: approximate standard errors based on the observed Hessian
#'   and the delta method.
#'   \item \code{convergence}: convergence code returned by \code{\link[stats]{optim}}.
#'   \item \code{counts}: number of function and gradient evaluations.
#'   \item \code{n_iter}: total number of reported optimizer evaluations.
#'   \item \code{hessian_ok}: logical indicator of whether the Hessian was finite.
#'   \item \code{vcov_ok}: logical indicator of whether the covariance matrix was obtained.
#'   \item \code{se_ok}: logical indicator of whether the standard errors were finite and positive.
#'   \item \code{hessian}: observed Hessian matrix on the log-parameter scale.
#'   \item \code{method}: optimization method associated with the best fit.
#'   \item \code{start}: starting values associated with the best fit.
#'   \item \code{message}: optimizer message.
#'   \item \code{warnings}: warnings generated during optimization.
#'   \item \code{success}: logical indicator of successful convergence.
#'   \item \code{n}: sample size.
#'   \item \code{k}: number of parameters.
#'   \item \code{AIC}: Akaike information criterion.
#'   \item \code{AICc}: corrected Akaike information criterion.
#'   \item \code{BIC}: Bayesian information criterion.
#'   \item \code{all_fits}: list with all optimization attempts.
#'   \item \code{x}: analyzed sample.
#'   \item \code{call}: matched function call.
#' }
#'
#' @details
#' The Redistributed-Beta distribution is defined through a redistribution
#' mechanism applied to a Beta\eqn{(\alpha,\beta)} baseline distribution.
#' The resulting model preserves support on \eqn{(0,1)} while allowing
#' flexible distributional shapes.
#'
#' @examples
#' set.seed(123)
#' x <- rRB(n = 200, alpha = 2, beta = 5)
#'
#' fit_mle <- fitRB_mle(x)
#' summary(fit_mle)
#'
#' @export
fitRB_mle <- function(x,
                      start = NULL,
                      starts = NULL,
                      methods = c("L-BFGS-B", "BFGS", "Nelder-Mead"),
                      multistart = TRUE,
                      lower_alpha = 1e-6,
                      upper_alpha = 1e6,
                      lower_beta = 1e-6,
                      upper_beta = 1e6,
                      control = list(maxit = 1000),
                      show_warnings = FALSE) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    stop("x must have length >= 2.", call. = FALSE)
  }
  if (any(x <= 0 | x >= 1)) {
    stop("x must lie strictly in (0,1).", call. = FALSE)
  }
  
  n <- length(x)
  start_list <- list()
  
  mom_fit <- try(fitRB_mom(x), silent = TRUE)
  
  if (!inherits(mom_fit, "try-error") &&
      isTRUE(mom_fit$success) &&
      !is.null(mom_fit$par)) {
    start_list[[length(start_list) + 1L]] <- mom_fit$par
  }
  
  tau1 <- log(2)
  
  x_corr <- x / tau1
  x_corr <- pmin(pmax(x_corr, 1e-6), 1 - 1e-6)
  
  m_corr <- mean(x_corr)
  v_corr <- stats::var(x_corr)
  
  if (is.finite(m_corr) && is.finite(v_corr) &&
      m_corr > 0 && m_corr < 1 && v_corr > 0) {
    
    phi_corr <- m_corr * (1 - m_corr) / v_corr - 1
    
    if (is.finite(phi_corr) && phi_corr > 0) {
      start_list[[length(start_list) + 1L]] <- c(
        alpha = m_corr * phi_corr,
        beta  = (1 - m_corr) * phi_corr
      )
    }
  }
  
  if (!is.null(start)) {
    start_list[[length(start_list) + 1L]] <- start
  }
  
  if (!is.null(starts)) {
    if (!is.list(starts)) {
      stop("starts must be a list of named numeric vectors.", call. = FALSE)
    }
    start_list <- c(start_list, starts)
  }
  
  neutral_starts <- list(
    c(alpha = 1, beta = 1),
    c(alpha = 2, beta = 2),
    c(alpha = 0.5, beta = 0.5),
    c(alpha = 2, beta = 5),
    c(alpha = 5, beta = 2)
  )
  
  start_list <- c(start_list, neutral_starts)
  
  if (isTRUE(multistart)) {
    grid_alpha <- c(0.5, 1, 2, 5, 10)
    grid_beta  <- c(0.5, 1, 2, 5, 10)
    
    grid_starts <- expand.grid(alpha = grid_alpha, beta = grid_beta)
    grid_starts <- split(grid_starts, seq_len(nrow(grid_starts)))
    
    grid_starts <- lapply(grid_starts, function(z) {
      c(alpha = z$alpha, beta = z$beta)
    })
    
    start_list <- c(start_list, grid_starts)
  }
  
  clean_start <- function(s) {
    
    if (is.null(s) || is.null(names(s)) ||
        !all(c("alpha", "beta") %in% names(s))) {
      return(NULL)
    }
    
    s <- as.numeric(s[c("alpha", "beta")])
    names(s) <- c("alpha", "beta")
    
    if (any(!is.finite(s)) || any(s <= 0)) {
      return(NULL)
    }
    
    s["alpha"] <- min(max(s["alpha"], lower_alpha), upper_alpha)
    s["beta"]  <- min(max(s["beta"], lower_beta), upper_beta)
    
    s
  }
  
  start_list <- lapply(start_list, clean_start)
  start_list <- start_list[!vapply(start_list, is.null, logical(1))]
  
  if (length(start_list) == 0L) {
    stop("No valid starting values available.", call. = FALSE)
  }
  
  keyfun <- function(s) paste(round(log(s), 6), collapse = "_")
  keys <- vapply(start_list, keyfun, character(1))
  start_list <- start_list[!duplicated(keys)]
  
  nll_eta <- function(eta) {
    
    alpha <- exp(eta[1])
    beta  <- exp(eta[2])
    
    ll <- sum(dRB(
      x,
      alpha = alpha,
      beta  = beta,
      log   = TRUE
    ))
    
    if (!is.finite(ll)) {
      return(1e12)
    }
    
    -ll
  }
  
  fit_one <- function(st, method) {
    
    eta0 <- log(c(st["alpha"], st["beta"]))
    warn_list <- character(0)
    
    res <- tryCatch(
      {
        withCallingHandlers(
          {
            if (method == "L-BFGS-B") {
              stats::optim(
                par = eta0,
                fn = nll_eta,
                method = method,
                lower = log(c(lower_alpha, lower_beta)),
                upper = log(c(upper_alpha, upper_beta)),
                hessian = TRUE,
                control = control
              )
            } else {
              stats::optim(
                par = eta0,
                fn = nll_eta,
                method = method,
                hessian = TRUE,
                control = control
              )
            }
          },
          warning = function(w) {
            warn_list <<- c(warn_list, conditionMessage(w))
            if (!isTRUE(show_warnings)) {
              invokeRestart("muffleWarning")
            }
          }
        )
      },
      error = function(e) e
    )
    
    if (inherits(res, "error")) {
      return(list(
        par = c(alpha = NA_real_, beta = NA_real_),
        logLik = -Inf,
        se = c(alpha = NA_real_, beta = NA_real_),
        convergence = NA_integer_,
        counts = c("function" = NA_integer_, "gradient" = NA_integer_),
        n_iter = NA_real_,
        hessian_ok = FALSE,
        vcov_ok = FALSE,
        se_ok = FALSE,
        hessian = NULL,
        method = method,
        start = st,
        message = conditionMessage(res),
        warnings = unique(warn_list),
        success = FALSE,
        AIC = Inf,
        AICc = Inf,
        BIC = Inf
      ))
    }
    
    par_hat <- exp(res$par)
    names(par_hat) <- c("alpha", "beta")
    
    logLik <- -res$value
    
    counts <- res$counts
    
    n_iter <- if (!is.null(counts)) {
      sum(counts, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    hessian_ok <- FALSE
    vcov_ok <- FALSE
    se_ok <- FALSE
    
    se <- c(alpha = NA_real_, beta = NA_real_)
    
    H <- try(as.matrix(res$hessian), silent = TRUE)
    
    if (!inherits(H, "try-error") &&
        all(dim(H) == c(2, 2)) &&
        all(is.finite(H))) {
      
      hessian_ok <- TRUE
      
      V <- try(solve(H), silent = TRUE)
      
      if (!inherits(V, "try-error") &&
          all(is.finite(V))) {
        
        vcov_ok <- TRUE
        
        se_eta <- sqrt(pmax(diag(V), 0))
        
        se <- par_hat * se_eta
        names(se) <- c("alpha", "beta")
        
        se_ok <- all(is.finite(se)) && all(se > 0)
      }
    }
    
    k <- 2L
    AIC  <- -2 * logLik + 2 * k
    BIC  <- -2 * logLik + log(n) * k
    AICc <- AIC + (2 * k * (k + 1)) / max(n - k - 1, 1)
    
    list(
      par = par_hat,
      logLik = logLik,
      se = se,
      convergence = res$convergence,
      method = method,
      counts = counts,
      n_iter = n_iter,
      hessian_ok = hessian_ok,
      vcov_ok = vcov_ok,
      se_ok = se_ok,
      hessian = H,
      start = st,
      message = res$message,
      warnings = unique(warn_list),
      success = is.finite(logLik) &&
        all(is.finite(par_hat)) &&
        all(par_hat > 0) &&
        res$convergence == 0,
      AIC = AIC,
      AICc = AICc,
      BIC = BIC
    )
  }
  
  all_fits <- list()
  id <- 1L
  
  for (st in start_list) {
    for (meth in methods) {
      all_fits[[id]] <- fit_one(st, meth)
      id <- id + 1L
    }
  }
  
  valid <- vapply(all_fits, function(z) is.finite(z$logLik), logical(1))
  
  if (!any(valid)) {
    stop("All optimization attempts failed.", call. = FALSE)
  }
  
  success_vec <- vapply(all_fits, function(z) isTRUE(z$success), logical(1))
  loglik_vec  <- vapply(all_fits, function(z) z$logLik, numeric(1))
  
  ord <- order(success_vec, loglik_vec, decreasing = TRUE)
  all_fits <- all_fits[ord]
  best <- all_fits[[1]]
  
  all_warnings <- unique(unlist(lapply(all_fits, function(z) z$warnings)))
  
  out <- list(
    par = best$par,
    logLik = best$logLik,
    se = best$se,
    convergence = best$convergence,
    counts = best$counts,
    n_iter = best$n_iter,
    hessian_ok = best$hessian_ok,
    vcov_ok = best$vcov_ok,
    se_ok = best$se_ok,
    hessian = best$hessian,
    method = best$method,
    start = best$start,
    message = best$message,
    warnings = all_warnings,
    success = best$success,
    n = n,
    k = 2L,
    AIC = best$AIC,
    AICc = best$AICc,
    BIC = best$BIC,
    all_fits = all_fits,
    x = x,
    call = match.call()
  )
  
  class(out) <- c("fitRB_mle", "list")
  out
}