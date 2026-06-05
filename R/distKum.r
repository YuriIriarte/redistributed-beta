# ============================================
# Kumaraswamy distribution
# ============================================

dKum <- function(x, a, b, log = FALSE) {
  x <- as.numeric(x)
  
  if (length(a) != 1 || length(b) != 1) {
    stop("'a' and 'b' must be scalars.", call. = FALSE)
  }
  
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    stop("'a' and 'b' must be positive finite numbers.", call. = FALSE)
  }
  
  logdens <- rep(-Inf, length(x))
  idx <- is.finite(x) & x > 0 & x < 1
  
  if (any(idx)) {
    logdens[idx] <- log(a) + log(b) +
      (a - 1) * log(x[idx]) +
      (b - 1) * log1p(-x[idx]^a)
  }
  
  if (log) return(logdens)
  
  dens <- exp(logdens)
  dens[!is.finite(dens)] <- 0
  dens
}

pKum <- function(q, a, b) {
  q <- as.numeric(q)
  
  if (length(a) != 1 || length(b) != 1) {
    stop("'a' and 'b' must be scalars.", call. = FALSE)
  }
  
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    stop("'a' and 'b' must be positive finite numbers.", call. = FALSE)
  }
  
  out <- rep(NA_real_, length(q))
  
  out[q <= 0] <- 0
  out[q >= 1] <- 1
  
  idx <- is.finite(q) & q > 0 & q < 1
  
  if (any(idx)) {
    out[idx] <- 1 - (1 - q[idx]^a)^b
  }
  
  out
}


# ============================================
# Maximum likelihood estimation
# ============================================

fitKum_mle <- function(x,
                       start = NULL,
                       starts = NULL,
                       methods = c("L-BFGS-B", "BFGS", "Nelder-Mead"),
                       multistart = TRUE,
                       lower_a = 1e-8,
                       upper_a = 1e8,
                       lower_b = 1e-8,
                       upper_b = 1e8,
                       control = list()) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    stop("x must have length >= 2.", call. = FALSE)
  }
  
  if (any(x <= 0 | x >= 1)) {
    stop("All observations must lie in (0,1).", call. = FALSE)
  }
  
  n <- length(x)
  
  start_list <- list()
  
  if (!is.null(start)) {
    start_list[[length(start_list) + 1L]] <- start
  }
  
  if (isTRUE(multistart)) {
    
    grid_a <- c(0.3, 0.5, 1, 2, 5, 10)
    grid_b <- c(0.3, 0.5, 1, 2, 5, 10)
    
    grid_starts <- expand.grid(a = grid_a, b = grid_b)
    grid_starts <- split(grid_starts, seq_len(nrow(grid_starts)))
    
    grid_starts <- lapply(grid_starts, function(z) {
      c(a = z$a, b = z$b)
    })
    
    start_list <- c(start_list, grid_starts)
    
  } else if (length(start_list) == 0L) {
    start_list <- list(c(a = 1, b = 1))
  }
  
  if (!is.null(starts)) {
    if (!is.list(starts)) {
      stop("starts must be a list of named numeric vectors.", call. = FALSE)
    }
    start_list <- c(start_list, starts)
  }
  
  clean_start <- function(s) {
    
    if (is.null(s) || is.null(names(s)) ||
        !all(c("a", "b") %in% names(s))) {
      return(NULL)
    }
    
    s <- as.numeric(s[c("a", "b")])
    names(s) <- c("a", "b")
    
    if (any(!is.finite(s)) || any(s <= 0)) {
      return(NULL)
    }
    
    s["a"] <- min(max(s["a"], lower_a), upper_a)
    s["b"] <- min(max(s["b"], lower_b), upper_b)
    
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
    
    a <- exp(eta[1])
    b <- exp(eta[2])
    
    ll <- sum(dKum(x, a = a, b = b, log = TRUE))
    
    if (!is.finite(ll)) {
      return(1e12)
    }
    
    -ll
  }
  
  fit_one <- function(st, method) {
    
    eta0 <- log(c(st["a"], st["b"]))
    
    res <- tryCatch(
      {
        if (method == "L-BFGS-B") {
          stats::optim(
            par = eta0,
            fn = nll_eta,
            method = method,
            lower = log(c(lower_a, lower_b)),
            upper = log(c(upper_a, upper_b)),
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
      error = function(e) e
    )
    
    if (inherits(res, "error")) {
      return(list(
        par = c(a = NA_real_, b = NA_real_),
        logLik = -Inf,
        se = c(a = NA_real_, b = NA_real_),
        convergence = NA_integer_,
        method = method,
        start = st,
        message = conditionMessage(res),
        success = FALSE,
        AIC = Inf,
        BIC = Inf
      ))
    }
    
    par_hat <- exp(res$par)
    names(par_hat) <- c("a", "b")
    
    logLik <- -res$value
    
    se <- c(a = NA_real_, b = NA_real_)
    
    H <- try(as.matrix(res$hessian), silent = TRUE)
    
    if (!inherits(H, "try-error") &&
        all(dim(H) == c(2, 2)) &&
        all(is.finite(H))) {
      
      V <- try(solve(H), silent = TRUE)
      
      if (!inherits(V, "try-error") &&
          all(is.finite(V))) {
        se_eta <- sqrt(pmax(diag(V), 0))
        se <- par_hat * se_eta
        names(se) <- c("a", "b")
      }
    }
    
    list(
      par = par_hat,
      logLik = logLik,
      se = se,
      convergence = res$convergence,
      method = method,
      start = st,
      message = res$message,
      success = is.finite(logLik) &&
        all(is.finite(par_hat)) &&
        all(par_hat > 0) &&
        res$convergence == 0,
      AIC = -2 * logLik + 2 * 2,
      BIC = -2 * logLik + log(n) * 2
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

  ord <- order(
    success_vec,
    loglik_vec,
    decreasing = TRUE
  )

  all_fits <- all_fits[ord]
  best <- all_fits[[1]]
  
  out <- list(
    par = best$par,
    logLik = best$logLik,
    se = best$se,
    convergence = best$convergence,
    method = best$method,
    start = best$start,
    message = best$message,
    success = best$success,
    n = n,
    k = 2L,
    AIC = best$AIC,
    BIC = best$BIC,
    all_fits = all_fits,
    call = match.call()
  )
  
  class(out) <- c("fitKum", "list")
  out
}