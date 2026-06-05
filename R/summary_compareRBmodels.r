#' @export
summary.compareRBmodels <- function(object, digits = 4, ...) {
  
  out <- list(
    call = object$call,
    n = object$n,
    comparison = object$comparison,
    estimates = object$estimates,
    best_AIC = object$best_AIC,
    best_AICc = object$best_AICc,
    best_BIC = object$best_BIC,
    digits = digits
  )
  
  class(out) <- c("summary.compareRBmodels", "list")
  out
}


#' @export
print.summary.compareRBmodels <- function(x, digits = x$digits, ...) {
  
  cat("\n")
  cat("Comparison of Fitted Models for Unit Data\n")
  cat(strrep("-", 42), "\n\n", sep = "")
  
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  
  cat("Sample size:", x$n, "\n\n")
  
  # ----------------------------------------------------------
  # Comparison table
  # ----------------------------------------------------------
  
  cat("Model comparison:\n")
  
  comp <- x$comparison
  
  if ("success" %in% names(comp)) {
    comp$success <- NULL
  }
  
  comp$logLik <- round(comp$logLik, digits)
  comp$AIC    <- round(comp$AIC, digits)
  
  if ("AICc" %in% names(comp)) {
    comp$AICc <- round(comp$AICc, digits)
  }
  
  comp$BIC <- round(comp$BIC, digits)
  
  print(comp, row.names = FALSE)
  
  # ----------------------------------------------------------
  # Parameter estimates
  # ----------------------------------------------------------
  
  cat("\nParameter estimates:\n")
  
  est <- x$estimates
  est$estimate <- round(est$estimate, digits)
  est$se       <- round(est$se, digits)
  
  print(est, row.names = FALSE)
  
  # ----------------------------------------------------------
  # Best models
  # ----------------------------------------------------------
  
  cat("\nModel selection summary:\n")
  cat("Best model according to AIC:  ", x$best_AIC, "\n", sep = "")
  
  if (!is.null(x$best_AICc)) {
    cat("Best model according to AICc: ", x$best_AICc, "\n", sep = "")
  }
  
  cat("Best model according to BIC:  ", x$best_BIC, "\n", sep = "")
  
  invisible(x)
}