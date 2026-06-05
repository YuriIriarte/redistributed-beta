# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(grid)
library(RB)

# ------------------------------------------------------------
# Global settings
# ------------------------------------------------------------

x_values <- seq(0.001, 0.999, length.out = 1000)

theme_rb <- theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill = "grey85", color = NA),
    strip.text = element_text(face = "bold", size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    legend.key.width = unit(1.6, "cm"),
    legend.key.height = unit(0.45, "cm")
  )

# ------------------------------------------------------------
# Scenarios (6 panels)
# ------------------------------------------------------------

scenarios <- data.frame(
  alpha = c(0.5, 1, 2, 1, 3, 4),
  beta  = c(0.5, 1, 1, 2, 4, 3),
  label = c(
    "alpha==0.5*','~beta==0.5",
    "alpha==1*','~beta==1",
    "alpha==2*','~beta==1",
    "alpha==1*','~beta==2",
    "alpha==3*','~beta==4",
    "alpha==4*','~beta==3"
  )
)

# ------------------------------------------------------------
# Data construction
# ------------------------------------------------------------

dens_data <- bind_rows(
  lapply(seq_len(nrow(scenarios)), function(i) {
    
    a <- scenarios$alpha[i]
    b <- scenarios$beta[i]
    lab <- scenarios$label[i]
    
    data.frame(
      t = rep(x_values, 2),
      density = c(
        dRB(x_values, alpha = a, beta = b),
        dbeta(x_values, shape1 = a, shape2 = b)
      ),
      Model = rep(c("RB", "Beta baseline"), each = length(x_values)),
      Panel = lab
    )
  })
) %>%
  mutate(
    Model = factor(Model, levels = c("RB", "Beta baseline")),
    Panel = factor(Panel, levels = scenarios$label)
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig_density <- ggplot(
  dens_data,
  aes(x = t, y = density, color = Model, linetype = Model)
) +
  geom_line(linewidth = 0.9) +
  facet_wrap(
    ~ Panel,
    scales = "free_y",
    ncol = 3,
    labeller = label_parsed
  ) +
  scale_color_manual(
    values = c("black", "black")
  ) +
  scale_linetype_manual(
    values = c("solid", "dashed")
  ) +
  labs(
    x = expression(t),
    y = "Density function"
  ) +
  theme_rb

fig_density