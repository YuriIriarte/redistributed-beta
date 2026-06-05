skewBeta <- function(shape1, shape2) {
  
  if (length(shape1) != 1L || length(shape2) != 1L) {
    stop("shape1 and shape2 must be scalars.", call. = FALSE)
  }
  
  if (!is.numeric(shape1) || !is.numeric(shape2) ||
      !is.finite(shape1) || !is.finite(shape2) ||
      shape1 <= 0 || shape2 <= 0) {
    stop("shape1 and shape2 must be positive finite scalars.", call. = FALSE)
  }
  
  a <- shape1
  b <- shape2
  
  2 * (b - a) * sqrt(a + b + 1) /
    ((a + b + 2) * sqrt(a * b))
}


kurtBeta <- function(shape1, shape2, excess = TRUE) {
  
  if (length(shape1) != 1L || length(shape2) != 1L) {
    stop("shape1 and shape2 must be scalars.", call. = FALSE)
  }
  
  if (!is.numeric(shape1) || !is.numeric(shape2) ||
      !is.finite(shape1) || !is.finite(shape2) ||
      shape1 <= 0 || shape2 <= 0) {
    stop("shape1 and shape2 must be positive finite scalars.", call. = FALSE)
  }
  
  if (length(excess) != 1L || !is.logical(excess)) {
    stop("excess must be a logical scalar.", call. = FALSE)
  }
  
  a <- shape1
  b <- shape2
  
  excess_kurt <- 6 * (
    (a - b)^2 * (a + b + 1) - a * b * (a + b + 2)
  ) / (a * b * (a + b + 2) * (a + b + 3))
  
  if (excess) {
    return(excess_kurt)
  }
  
  excess_kurt + 3
}

fitBeta_mle <- function(x,
                        start = NULL,
                        starts = NULL,
                        methods = c("L-BFGS-B", "BFGS", "Nelder-Mead"),
                        multistart = TRUE,
                        lower_shape1 = 1e-8,
                        upper_shape1 = 1e8,
                        lower_shape2 = 1e-8,
                        upper_shape2 = 1e8,
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
  
  # ----------------------------------------------------------
  # Starting values
  # ----------------------------------------------------------
  start_list <- list()
  
  if (!is.null(start)) {
    start_list[[length(start_list) + 1L]] <- start
  }
  
  # MOM start
  m <- mean(x)
  v <- stats::var(x)
  
  if (is.finite(m) && is.finite(v) && v > 0 && v < m * (1 - m)) {
    tmp <- m * (1 - m) / v - 1
    mom_start <- c(
      shape1 = m * tmp,
      shape2 = (1 - m) * tmp
    )
    
    if (all(is.finite(mom_start)) && all(mom_start > 0)) {
      start_list[[length(start_list) + 1L]] <- mom_start
    }
  }
  
  # Deterministic grid
  if (isTRUE(multistart)) {
    grid_shape1 <- c(0.3, 0.5, 1, 2, 5, 10)
    grid_shape2 <- c(0.3, 0.5, 1, 2, 5, 10)
    
    grid_starts <- expand.grid(
      shape1 = grid_shape1,
      shape2 = grid_shape2
    )
    
    grid_starts <- split(grid_starts, seq_len(nrow(grid_starts)))
    
    grid_starts <- lapply(grid_starts, function(z) {
      c(shape1 = z$shape1, shape2 = z$shape2)
    })
    
    start_list <- c(start_list, grid_starts)
    
  } else if (length(start_list) == 0L) {
    start_list <- list(c(shape1 = 1, shape2 = 1))
  }
  
  if (!is.null(starts)) {
    if (!is.list(starts)) {
      stop("starts must be a list of named numeric vectors.", call. = FALSE)
    }
    start_list <- c(start_list, starts)
  }
  
  clean_start <- function(s) {
    
    if (is.null(s)) {
      return(NULL)
    }
    
    if (is.null(names(s)) || 
        !all(c("shape1", "shape2") %in% names(s))) {
      return(NULL)
    }
    
    s <- as.numeric(s[c("shape1", "shape2")])
    names(s) <- c("shape1", "shape2")
    
    if (any(!is.finite(s)) || any(s <= 0)) {
      return(NULL)
    }
    
    s["shape1"] <- min(max(s["shape1"], lower_shape1), upper_shape1)
    s["shape2"] <- min(max(s["shape2"], lower_shape2), upper_shape2)
    
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
  
  # ----------------------------------------------------------
  # Negative log-likelihood in log-scale
  # ----------------------------------------------------------
  nll_eta <- function(eta) {
    
    shape1 <- exp(eta[1])
    shape2 <- exp(eta[2])
    
    ll <- sum(stats::dbeta(
      x,
      shape1 = shape1,
      shape2 = shape2,
      log = TRUE
    ))
    
    if (!is.finite(ll)) {
      return(1e12)
    }
    
    -ll
  }
  
  # ----------------------------------------------------------
  # Fit one optimization
  # ----------------------------------------------------------
  fit_one <- function(st, method) {
    
    eta0 <- log(c(st["shape1"], st["shape2"]))
    
    res <- tryCatch(
      {
        if (method == "L-BFGS-B") {
          stats::optim(
            par = eta0,
            fn = nll_eta,
            method = method,
            lower = log(c(lower_shape1, lower_shape2)),
            upper = log(c(upper_shape1, upper_shape2)),
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
        par = c(shape1 = NA_real_, shape2 = NA_real_),
        logLik = -Inf,
        se = c(shape1 = NA_real_, shape2 = NA_real_),
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
    names(par_hat) <- c("shape1", "shape2")
    
    logLik <- -res$value
    
    se <- c(shape1 = NA_real_, shape2 = NA_real_)
    
    H <- try(as.matrix(res$hessian), silent = TRUE)
    
    if (!inherits(H, "try-error") &&
        all(dim(H) == c(2, 2)) &&
        all(is.finite(H))) {
      
      V <- try(solve(H), silent = TRUE)
      
      if (!inherits(V, "try-error") &&
          all(is.finite(V))) {
        
        se_eta <- sqrt(pmax(diag(V), 0))
        se <- par_hat * se_eta
        names(se) <- c("shape1", "shape2")
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
  
  # ----------------------------------------------------------
  # Run all fits
  # ----------------------------------------------------------
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
  
  class(out) <- c("fitBeta", "list")
  out
}