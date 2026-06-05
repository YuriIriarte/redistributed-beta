# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("Ecdat")
# remotes::install_github("YuriIriarte/redistributed-beta")

library(ggplot2)
library(dplyr)
library(tidyr)
library(RB)
library(Ecdat)

# ------------------------------------------------------------
# Budget share for other good expenditure
# ------------------------------------------------------------

x <- BudgetUK[, 6]

# ------------------------------------------------------------
# Fit competing models
# ------------------------------------------------------------

cmp <- compareRBmodels(x)

summary(cmp)

# Comparison of Fitted Models for Unit Data
# ------------------------------------------
#
# Call:
# compareRBmodels(x = x)
#
# Sample size: 1519
#
# Model comparison:
#  Model npar   logLik       AIC      AICc       BIC convergence
#     RB    2 1389.898 -2775.796 -2775.788 -2765.144           0
#      B    2 1373.147 -2742.293 -2742.285 -2731.642           0
#      K    2 1322.564 -2641.128 -2641.121 -2630.477           0
#
# Parameter estimates:
#  Model parameter estimate     se
#     RB     alpha   5.0973 0.2288
#     RB      beta   8.8818 0.4300
#      B    shape1   4.4476 0.1560
#      B    shape2  13.1181 0.4775
#      K         a   2.4957 0.0503
#      K         b  22.1273 1.3596
#
# Model selection summary:
# Best model according to AIC:  RB
# Best model according to AICc: RB
# Best model according to BIC:  RB

gof_fit <- gofRB_boot(
  x = x,
  B = 1000,
  method = "L-BFGS-B",
  multistart = FALSE,
  seed = 2026
)

summary(gof_fit)

# Bootstrap Goodness-of-Fit for the Redistributed-Beta Distribution
# -------------------------------------------------------------------
#
# Call:
# gofRB_boot(x = x, B = 1000, seed = 2026, method = "L-BFGS-B",
#     multistart = FALSE)
#
# Sample size: 1519
# Bootstrap replicates: 1000
# Valid bootstrap replicates: 1000
#
# Parameter estimates:
#  Parameter Estimate
#      alpha   5.0973
#       beta   8.8818
#
# Fit statistics:
#    logLik       AIC       BIC
#  1389.898 -2775.796 -2765.144
#
# Goodness-of-fit statistics:
#              Test Statistic p_value
#  Anderson-Darling    0.4829  0.2488
#  Cramer-von Mises    0.0591  0.4066

# ------------------------------------------------------------
# Fitted histogram
# ------------------------------------------------------------

dens_grid <- seq(0.01, 0.75, length.out = 600)

dens_df <- data.frame(
  x = dens_grid,
  B = dbeta(dens_grid, shape1 = 4.4476, shape2 = 13.1181),
  K = RB:::dKum(dens_grid, a = 2.4957, b = 22.1273),
  RB = dRB(dens_grid, alpha = 5.0973, beta = 8.8818)
) %>%
  pivot_longer(
    cols = c(B, K, RB),
    names_to = "Distribution",
    values_to = "Density"
  ) %>%
  mutate(
    Distribution = factor(
      Distribution,
      levels = c("RB", "B", "K")
    )
  )

ggplot(data.frame(x = x), aes(x = x)) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(0, 0.7, by = 0.04),
    fill = "white",
    color = "black",
    alpha = 0.65
  ) +
  geom_line(
    data = dens_df,
    aes(
      x = x,
      y = Density,
      linetype = Distribution,
      linewidth = Distribution
    ),
    color = "black"
  ) +
  scale_linetype_manual(
    values = c(
      RB = "solid",
      B = "dashed",
      K = "dotdash"
    )
  ) +
  scale_linewidth_manual(
    values = c(
      RB = 1.1,
      B = 0.9,
      K = 0.9
    )
  ) +
  labs(
    x = "Budget share",
    y = "Density function",
    linetype = "Distribution",
    linewidth = "Distribution"
  ) +
  guides(
    linewidth = "none",
    linetype = guide_legend(
      override.aes = list(linewidth = c(1.1, 0.9, 0.9))
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    
    legend.position = c(0.74, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      color = "grey70",
      linewidth = 0.3
    ),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.6, "cm"),
    legend.key.height = unit(0.45, "cm")
  )

# ------------------------------------------------------------
# Cumulative distributions
# ------------------------------------------------------------

alpha_rb <- 5.0973
beta_rb  <- 8.8818

grid_x <- seq(0.001, 0.8, length.out = 600)

rb_df <- data.frame(
  x = grid_x,
  CDF = pRB(grid_x, alpha = alpha_rb, beta = beta_rb)
)

ggplot() +
  stat_ecdf(
    data = data.frame(x = x),
    aes(x = x, linetype = "Empirical"),
    geom = "step",
    linewidth = 1,
    color = "black"
  ) +
  geom_line(
    data = rb_df,
    aes(x = x, y = CDF, linetype = "RB"),
    linewidth = 1,
    color = "black"
  ) +
  scale_linetype_manual(
    values = c(
      Empirical = "solid",
      RB = "longdash"
    )
  ) +
  labs(
    x = "Budget share",
    y = "Cumulative function",
    linetype = "Distribution"
  ) +
  coord_cartesian(xlim = c(0.01, 0.8), ylim = c(0, 1)) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      color = "grey70",
      linewidth = 0.3
    ),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.6, "cm"),
    legend.key.height = unit(0.45, "cm")
  )