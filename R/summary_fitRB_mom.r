#' @export
summary.fitRB_mom <- function(object, digits = 4, ...) {
  
  par_tab <- data.frame(
    Parameter = names(object$par),
    Estimate = as.numeric(object$par),
    row.names = NULL
  )
  
  mom_tab <- data.frame(
    Quantity = names(object$moments),
    Value = as.numeric(object$moments),
    row.names = NULL
  )
  
  out <- list(
    call = object$call,
    coefficients = par_tab,
    moments = mom_tab,
    success = object$success,
    digits = digits
  )
  
  class(out) <- c("summary.fitRB_mom", "list")
  out
}


#' @export
print.summary.fitRB_mom <- function(x, digits = x$digits, ...) {
  
  cat("\n")
  cat("Moment Estimation for the Redistributed-Beta Distribution\n")
  cat(strrep("-", 58), "\n\n", sep = "")
  
  cat("Call:\n")
  print(x$call)
  
  cat("\nParameter estimates:\n")
  coef_tab <- x$coefficients
  coef_tab$Estimate <- round(coef_tab$Estimate, digits)
  print(coef_tab, row.names = FALSE)
  
  cat("\nMoment-based quantities:\n")
  mom_tab <- x$moments
  mom_tab$Value <- round(mom_tab$Value, digits)
  print(mom_tab, row.names = FALSE)
  
  cat("\nConvergence status:\n")
  cat("Successful estimation:", ifelse(isTRUE(x$success), "yes", "no"), "\n")
  
  invisible(x)
}