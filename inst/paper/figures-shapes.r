# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------
# install.packages("plot3D")
# remotes::install_github("YuriIriarte/redistributed-beta")

library(RB)
library(plot3D)

# ------------------------------------------------------------
# Figure: Skewness 3D
# ------------------------------------------------------------

alpha_grid <- seq(0.05, 50, length.out = 20)
beta_grid  <- seq(0.05, 50, length.out = 20)

skew_fun <- Vectorize(function(alpha, beta) {
  skewRB(alpha = alpha, beta = beta)
})

z_skew <- outer(alpha_grid, beta_grid, skew_fun)

persp3D(
  x = alpha_grid, y = beta_grid, z = z_skew,
  xlab = "", ylab = "", zlab = "",
  theta = 130, phi = 20, d = 200,
  zlim = c(-10, 10),
  bty = "f",
  ticktype = "detailed",
  facets = FALSE,
  colkey = FALSE,
  col = "black",
  cex.axis = 0.7
)

skew_beta_fun <- Vectorize(function(alpha, beta) {
  RB:::skewBeta(shape1 = alpha, shape2 = beta)
})

z_skew_beta <- outer(alpha_grid, beta_grid, skew_beta_fun)

persp3D(
  x = beta_grid, y = alpha_grid, z = z_skew_beta,
  add = TRUE,
  facets = FALSE,
  colkey = FALSE,
  col = gray.colors(100),
  lwd = 0.7
)

text3D(
  x = 45, y = -15, z = -6.5,
  labels = "Skewness", srt = 90, add = TRUE,
  col = "black", cex = 1.2
)

text3D(
  x = 45, y = 10, z = -15,
  labels = expression(beta), add = TRUE,
  col = "black", cex = 1.2
)

text3D(
  x = 10, y = 42, z = -15,
  labels = expression(alpha), add = TRUE,
  col = "black", cex = 1.2
)

# ------------------------------------------------------------
# Figure: Kurtosis 3D
# ------------------------------------------------------------

alpha_grid <- seq(1, 50, length.out = 20)
beta_grid  <- seq(1, 50, length.out = 20)

kurt_fun <- Vectorize(function(alpha, beta) {
  kurtRB(alpha = alpha, beta = beta)
})

z_kurt <- outer(alpha_grid, beta_grid, kurt_fun)

persp3D(
  x = alpha_grid, y = beta_grid, z = z_kurt,
  xlab = "", ylab = "", zlab = "",
  theta = 150, phi = 20, d = 200,
  zlim = c(-1, 8.5),
  bty = "f",
  ticktype = "detailed",
  facets = FALSE,
  colkey = FALSE,
  col = "black",
  cex.axis = 0.7
)

kurt_beta_fun <- Vectorize(function(alpha, beta) {
  RB:::kurtBeta(shape1 = alpha, shape2 = beta, excess = FALSE)
})

z_kurt_beta <- outer(alpha_grid, beta_grid, kurt_beta_fun)

persp3D(
  x = alpha_grid, y = beta_grid, z = z_kurt_beta,
  add = TRUE,
  facets = FALSE,
  colkey = FALSE,
  col = gray.colors(100),
  lwd = 0.8
)

text3D(
  x = 46, y = -21, z = 0.5,
  labels = "Kurtosis", srt = 90, add = TRUE,
  col = "black", cex = 1.2
)

text3D(
  x = 45, y = 5, z = -4,
  labels = expression(beta), add = TRUE,
  col = "black", cex = 1.2
)

text3D(
  x = 2, y = 15, z = -5.2,
  labels = expression(alpha), add = TRUE,
  col = "black", cex = 1.2
)