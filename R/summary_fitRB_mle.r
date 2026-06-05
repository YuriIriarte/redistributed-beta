#' @export
summary.fitRB_mle <- function(object, digits = 4, ...) {
  
  par_tab <- data.frame(
    Parameter = names(object$par),
    Estimate  = as.numeric(object$par),
    SE        = as.numeric(object$se),
    row.names = NULL
  )
  
  fit_tab <- data.frame(
    logLik = object$logLik,
    AIC    = object$AIC,
    AICc   = object$AICc,
    BIC    = object$BIC,
    Method = object$method,
    stringsAsFactors = FALSE
  )
  
  out <- list(
    call = object$call,
    coefficients = par_tab,
    fit_statistics = fit_tab,
    convergence = object$convergence,
    success = object$success,
    n = object$n,
    k = object$k,
    digits = digits
  )
  
  class(out) <- c("summary.fitRB_mle", "list")
  out
}


#' @export
print.summary.fitRB_mle <- function(x, digits = x$digits, ...) {
  
  cat("\n")
  cat("Maximum Likelihood Estimation for the Redistributed-Beta Distribution\n")
  cat(strrep("-", 70), "\n\n", sep = "")
  
  cat("Call:\n")
  print(x$call)
  
  cat("\nSample size:", x$n, "\n")
  cat("Number of parameters:", x$k, "\n")
  
  cat("\nParameter estimates:\n")
  coef_tab <- x$coefficients
  coef_tab$Estimate <- round(coef_tab$Estimate, digits)
  coef_tab$SE <- round(coef_tab$SE, digits)
  print(coef_tab, row.names = FALSE)
  
  cat("\nFit statistics:\n")
  fit_tab <- x$fit_statistics
  fit_tab$logLik <- round(fit_tab$logLik, digits)
  fit_tab$AIC    <- round(fit_tab$AIC, digits)
  fit_tab$AICc   <- round(fit_tab$AICc, digits)
  fit_tab$BIC    <- round(fit_tab$BIC, digits)
  print(fit_tab, row.names = FALSE)
  
  cat("\nConvergence status:\n")
  cat("Convergence code:", x$convergence, "\n")
  cat("Successful estimation:", ifelse(isTRUE(x$success), "yes", "no"), "\n")
  
  invisible(x)
}