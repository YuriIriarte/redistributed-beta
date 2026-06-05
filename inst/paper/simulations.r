# ============================================================
# Monte Carlo study for the Redistributed-Beta (RB) distribution
# ============================================================

library(RB)
library(dplyr)
library(tidyr)
library(patchwork)

# ------------------------------------------------------------
# One Monte Carlo replicate
# ------------------------------------------------------------

mc_one_rb <- function(n, alpha_true, beta_true) {
  
  x <- rRB(n = n, alpha = alpha_true, beta = beta_true)
  
  # ----------------------------------------------------------
  # Moment estimator
  # ----------------------------------------------------------
  
  fit_mom <- try(fitRB_mom(x), silent = TRUE)
  
  if (!inherits(fit_mom, "try-error") &&
      isTRUE(fit_mom$success) &&
      !is.null(fit_mom$par) &&
      all(is.finite(fit_mom$par))) {
    
    alpha_mom <- unname(fit_mom$par["alpha"])
    beta_mom  <- unname(fit_mom$par["beta"])
    success_mom <- TRUE
    
  } else {
    
    alpha_mom <- NA_real_
    beta_mom  <- NA_real_
    success_mom <- FALSE
  }
  
  # ----------------------------------------------------------
  # Maximum likelihood estimator
  # ----------------------------------------------------------
  
  fit_mle <- try(
    fitRB_mle(
      x,
      methods = "L-BFGS-B",
      multistart = FALSE
    ),
    silent = TRUE
  )
  
  conv_code  <- NA_integer_
  n_iter     <- NA_real_
  hessian_ok <- FALSE
  vcov_ok    <- FALSE
  se_ok      <- FALSE
  ci_ok      <- FALSE
  
  if (inherits(fit_mle, "try-error") ||
      is.null(fit_mle$par) ||
      any(!is.finite(fit_mle$par)) ||
      is.null(fit_mle$convergence) ||
      fit_mle$convergence != 0) {
    
    return(data.frame(
      alpha_true = alpha_true,
      beta_true = beta_true,
      n = n,
      alpha_mle = NA_real_,
      beta_mle = NA_real_,
      alpha_mom = alpha_mom,
      beta_mom = beta_mom,
      se_alpha = NA_real_,
      se_beta = NA_real_,
      cp_alpha = NA_real_,
      cp_beta = NA_real_,
      success_mle = FALSE,
      success_mom = success_mom,
      conv_code = conv_code,
      n_iter = n_iter,
      hessian_ok = hessian_ok,
      vcov_ok = vcov_ok,
      se_ok = se_ok,
      ci_ok = ci_ok
    ))
  }
  
  conv_code  <- fit_mle$convergence
  n_iter     <- fit_mle$n_iter
  hessian_ok <- fit_mle$hessian_ok
  vcov_ok    <- fit_mle$vcov_ok
  se_ok      <- fit_mle$se_ok
  
  alpha_mle <- unname(fit_mle$par["alpha"])
  beta_mle  <- unname(fit_mle$par["beta"])
  
  se_alpha <- if (!is.null(fit_mle$se["alpha"])) {
    unname(fit_mle$se["alpha"])
  } else {
    NA_real_
  }
  
  se_beta <- if (!is.null(fit_mle$se["beta"])) {
    unname(fit_mle$se["beta"])
  } else {
    NA_real_
  }
  
  cp_alpha <- NA_real_
  cp_beta  <- NA_real_
  
  if (is.finite(se_alpha) && se_alpha > 0) {
    li_alpha <- alpha_mle - 1.96 * se_alpha
    ls_alpha <- alpha_mle + 1.96 * se_alpha
    cp_alpha <- as.numeric(alpha_true >= li_alpha && alpha_true <= ls_alpha)
  }
  
  if (is.finite(se_beta) && se_beta > 0) {
    li_beta <- beta_mle - 1.96 * se_beta
    ls_beta <- beta_mle + 1.96 * se_beta
    cp_beta <- as.numeric(beta_true >= li_beta && beta_true <= ls_beta)
  }
  
  ci_ok <- is.finite(cp_alpha) && is.finite(cp_beta)
  
  data.frame(
    alpha_true = alpha_true,
    beta_true = beta_true,
    n = n,
    alpha_mle = alpha_mle,
    beta_mle = beta_mle,
    alpha_mom = alpha_mom,
    beta_mom = beta_mom,
    se_alpha = se_alpha,
    se_beta = se_beta,
    cp_alpha = cp_alpha,
    cp_beta = cp_beta,
    success_mle = TRUE,
    success_mom = success_mom,
    conv_code = conv_code,
    n_iter = n_iter,
    hessian_ok = hessian_ok,
    vcov_ok = vcov_ok,
    se_ok = se_ok,
    ci_ok = ci_ok
  )
}

# ------------------------------------------------------------
# Full Monte Carlo study
# ------------------------------------------------------------

mc_rb_study <- function(R = 1000,
                        n_values = c(50, 100, 200, 300, 500),
                        alpha_values = c(0.5, 2, 4),
                        beta_values = c(1, 3, 6),
                        seed = 123,
                        verbose = TRUE) {
  
  set.seed(seed)
  
  scenarios <- expand.grid(
    alpha_true = alpha_values,
    beta_true = beta_values,
    n = n_values
  )
  
  results <- vector("list", nrow(scenarios) * R)
  id <- 1L
  
  for (i in seq_len(nrow(scenarios))) {
    
    alpha_true <- scenarios$alpha_true[i]
    beta_true  <- scenarios$beta_true[i]
    n          <- scenarios$n[i]
    
    if (isTRUE(verbose)) {
      message(
        "Scenario ", i, "/", nrow(scenarios),
        ": alpha = ", alpha_true,
        ", beta = ", beta_true,
        ", n = ", n
      )
    }
    
    for (r in seq_len(R)) {
      
      aux <- mc_one_rb(
        n = n,
        alpha_true = alpha_true,
        beta_true = beta_true
      )
      
      aux$replicate <- r
      results[[id]] <- aux
      id <- id + 1L
    }
  }
  
  bind_rows(results) %>%
    select(
      replicate,
      alpha_true,
      beta_true,
      n,
      alpha_mle,
      beta_mle,
      alpha_mom,
      beta_mom,
      se_alpha,
      se_beta,
      cp_alpha,
      cp_beta,
      success_mle,
      success_mom,
      conv_code,
      n_iter,
      hessian_ok,
      vcov_ok,
      se_ok,
      ci_ok
    )
}

# ------------------------------------------------------------
# Statistical summary:
# Bias, RMSE, CP
# ------------------------------------------------------------

summarize_rb_mc <- function(mc_results) {
  
  mc_results %>%
    group_by(alpha_true, beta_true, n) %>%
    summarise(
      Bias_alpha_MLE = mean(alpha_mle - alpha_true, na.rm = TRUE),
      Bias_beta_MLE  = mean(beta_mle  - beta_true,  na.rm = TRUE),
      RMSE_alpha_MLE = sqrt(mean((alpha_mle - alpha_true)^2, na.rm = TRUE)),
      RMSE_beta_MLE  = sqrt(mean((beta_mle  - beta_true)^2,  na.rm = TRUE)),
      
      Bias_alpha_MoM = mean(alpha_mom - alpha_true, na.rm = TRUE),
      Bias_beta_MoM  = mean(beta_mom  - beta_true,  na.rm = TRUE),
      RMSE_alpha_MoM = sqrt(mean((alpha_mom - alpha_true)^2, na.rm = TRUE)),
      RMSE_beta_MoM  = sqrt(mean((beta_mom  - beta_true)^2,  na.rm = TRUE)),
      
      CP_alpha_MLE = mean(cp_alpha, na.rm = TRUE),
      CP_beta_MLE  = mean(cp_beta,  na.rm = TRUE),
      
      SuccessRate_MLE = mean(success_mle, na.rm = TRUE),
      SuccessRate_MoM = mean(success_mom, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# Computational stability summary
# ------------------------------------------------------------

summarize_rb_comp <- function(mc_results) {
  
  mc_results %>%
    group_by(alpha_true, beta_true, n) %>%
    summarise(
      ConvRate_MLE    = mean(success_mle, na.rm = TRUE),
      HessianRate_MLE = mean(hessian_ok, na.rm = TRUE),
      VcovRate_MLE    = mean(vcov_ok, na.rm = TRUE),
      SERate_MLE      = mean(se_ok, na.rm = TRUE),
      CIRate_MLE      = mean(ci_ok, na.rm = TRUE),
      AvgIter_MLE     = mean(n_iter, na.rm = TRUE),
      SuccessRate_MoM = mean(success_mom, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# Run the Monte Carlo study
# ------------------------------------------------------------

rb_mc_res <- mc_rb_study(
  R = 1000,
  n_values = c(50, 100, 200, 300, 500),
  alpha_values = c(0.5, 2, 4),
  beta_values = c(1, 3, 6),
  seed = 123,
  verbose = TRUE
)

# ------------------------------------------------------------
# Monte Carlo summaries
# ------------------------------------------------------------

rb_mc_summary_full <- summarize_rb_mc(rb_mc_res)
rb_mc_comp_full    <- summarize_rb_comp(rb_mc_res)

rb_mc_summary_alpha05 <- rb_mc_summary_full %>%
  filter(alpha_true == 0.5) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

rb_mc_summary_alpha2 <- rb_mc_summary_full %>%
  filter(alpha_true == 2) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

rb_mc_summary_alpha4 <- rb_mc_summary_full %>%
  filter(alpha_true == 4) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

rb_mc_comp_alpha05 <- rb_mc_comp_full %>%
  filter(alpha_true == 0.5) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

rb_mc_comp_alpha2 <- rb_mc_comp_full %>%
  filter(alpha_true == 2) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

rb_mc_comp_alpha4 <- rb_mc_comp_full %>%
  filter(alpha_true == 4) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

as.data.frame(rb_mc_summary_alpha05)
as.data.frame(rb_mc_summary_alpha2)
as.data.frame(rb_mc_summary_alpha4)

as.data.frame(rb_mc_comp_alpha05)
as.data.frame(rb_mc_comp_alpha2)
as.data.frame(rb_mc_comp_alpha4)