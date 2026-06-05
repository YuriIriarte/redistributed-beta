# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("nasapower")
# remotes::install_github("YuriIriarte/redistributed-beta")

library(ggplot2)
library(nasapower)
library(dplyr)
library(tidyr)
library(RB)

# ------------------------------------------------------------
# Weekly relative humidity data for Calama, Chile
# Source: NASA POWER
# ------------------------------------------------------------

get_weekly_RH_Calama <- function(
    lon = -68.93,
    lat = -22.46,
    start_date = "2010-01-01",
    end_date = "2025-12-31",
    file = "humedadsemanal.csv"
) {
  
  months_es <- c(
    "enero", "febrero", "marzo", "abril",
    "mayo", "junio", "julio", "agosto",
    "septiembre", "octubre", "noviembre", "diciembre"
  )
  
  rh_data <- get_power(
    community = "AG",
    lonlat = c(lon, lat),
    pars = "RH2M",
    dates = c(start_date, end_date),
    temporal_api = "daily"
  )
  
  data_meses <- rh_data %>%
    mutate(
      semana = case_when(
        DD <= 7  ~ 1L,
        DD <= 14 ~ 2L,
        DD <= 21 ~ 3L,
        TRUE     ~ 4L
      ),
      mes = factor(MM, levels = 1:12, labels = months_es)
    ) %>%
    group_by(YEAR, mes, semana) %>%
    summarise(
      RH_promedio = mean(RH2M, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(YEAR, mes, semana) %>%
    group_by(mes) %>%
    mutate(id = row_number()) %>%
    ungroup() %>%
    select(id, mes, RH_promedio) %>%
    pivot_wider(
      names_from = mes,
      values_from = RH_promedio
    ) %>%
    arrange(id)
  
  data_meses
}

data_meses <- get_weekly_RH_Calama()

# ------------------------------------------------------------
# Fit competing models
# ------------------------------------------------------------

x <- data_meses$marzo / 100

cmp <- compareRBmodels(x)

summary(cmp)

# Comparison of Fitted Models for Unit Data
# ------------------------------------------
#
# Call:
# compareRBmodels(x = x)
#
# Sample size: 64
#
# Model comparison:
#  Model npar  logLik       AIC      AICc       BIC convergence
#     RB    2 66.8683 -129.7366 -129.5399 -125.4188           0
#      B    2 64.5587 -125.1174 -124.9206 -120.7996           0
#      K    2 62.2575 -120.5149 -120.3182 -116.1971           0
#
# Parameter estimates:
#  Model parameter estimate      se
#     RB     alpha  20.0397  7.3745
#     RB      beta  23.3766  8.9340
#      B    shape1   8.6440  1.5017
#      B    shape2  18.0168  3.1778
#      K         a   3.8019  0.3754
#      K         b  48.6629 17.9120
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
# Sample size: 64
# Bootstrap replicates: 1000
# Valid bootstrap replicates: 1000
#
# Parameter estimates:
#  Parameter Estimate
#      alpha  20.0397
#       beta  23.3766
#
# Fit statistics:
#   logLik       AIC       BIC
#  66.8683 -129.7366 -125.4188
#
# Goodness-of-fit statistics:
#              Test Statistic p_value
#  Anderson-Darling    0.7429  0.1419
#  Cramer-von Mises    0.1196  0.1738

# ------------------------------------------------------------
# Fitted histogram
# ------------------------------------------------------------

dens_grid <- seq(0.1, 0.65, length.out = 600)

dens_df <- data.frame(
  x = dens_grid,
  B = dbeta(dens_grid, shape1 = 8.6440, shape2 = 18.0168),
  K = RB:::dKum(dens_grid, a = 3.8019, b = 48.6629),
  RB = dRB(dens_grid, alpha = 20.0397, beta = 23.3766)
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
    breaks = seq(0.15, 0.56, by = 0.05),
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
    x = "Relative humidity",
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
# Empirical CDF vs fitted RB CDF
# ------------------------------------------------------------

alpha_rb <- 20.0397
beta_rb  <- 23.3766

grid_x <- seq(0.01, 0.75, length.out = 600)

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
    x = "Relative humidity",
    y = "Cumulative function",
    linetype = "Distribution"
  ) +
  coord_cartesian(xlim = c(0.01, 0.75), ylim = c(0, 1)) +
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