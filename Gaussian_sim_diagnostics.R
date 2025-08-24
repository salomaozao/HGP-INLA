rm(list = ls())
set.seed(1232)

source("./blockNNGPrgeneric.R")
source("C:/Users/GabrielNascimento/Documents/Gabriel/IC/TESTE/HGP-INLA/utils.R")
library(tidyverse)
library(sf)
library(INLA)

diagnosticar_inferencia <- function(out, true_values, n_star, M) {
  # Nomes dos parâmetros
  parametros <- c("Intercept", "Beta", "sigma.sq", "phi", "tau.sq")

  # Calcula métricas de diagnóstico
  means <- colMeans(out)
  bias <- means - true_values
  variance <- apply(out, 2, var)
  mse <- bias^2 + variance

  # Cria dataframe com resultados
  resultados <- data.frame(
    Parâmetro = parametros,
    Valor_Verdadeiro = true_values,
    Média_Estimativa = round(means, 4),
    Viés = round(bias, 4),
    Variância = round(variance, 4),
    EQM = round(mse, 4)
  )

  # Exibe informações da amostra
  cat("DIAGNÓSTICO DA INFERÊNCIA\n")
  cat("=========================\n")
  cat("Tamanho da amostra (n):", n_star, "\n")
  cat("Número de iterações (M):", M, "\n\n")

  # Exibe tabela de resultados
  print(resultados, row.names = FALSE)
  cat("\n")

  # Exibe estatísticas descritivas
  cat("Estatísticas Descritivas das Estimativas:\n")
  print(summary(out))
}

n_total <- 1000
sample_sizes = c(100, 200, 500)

#  pass spatial parameters
B <- as.matrix(c(1, 5))
tau.sq <- 0.1
sigma.sq <- 1
phi <- 1 / 3
alpha <- 1

priors = list(a = 0, b = 2.432049)
n.partition <- 8
n.blocks <- n.partition^2
num.nb <- 2

M <- 3


loc <- cbind(runif(n_total, 0, 1), runif(n_total, 0, 1))
colnames(loc) <- c("x", "y")

sf <- st_as_sf(
  as.data.frame(loc),
  coords = c("x", "y"),
  crs = 4326
)
loc <- st_coordinates(st_centroid(sf))


# range <- log(10)^(1 / alpha) * phi

distMatrix <- st_distance(sf, which = "Hausdorff") |>
  units::set_units(value = NULL) |>
  as.matrix()
R <- exp((-1) * (distMatrix / phi)^alpha)


diag(R) <- 1

C <- sigma.sq * R

D <- chol(C)
rnorm_n.obs <- rnorm(n_total)
w <- t(matrix(rnorm_n.obs, ncol = (n_total)) %*% D)

X <- as.matrix(cbind(1, rnorm(n_total))) ## X = intercept + covariate

p <- length(B)

y <- rnorm(n_total, X %*% B + w, sqrt(tau.sq)) ## y= X beta + w(spatial) + nugget

true_values <- c(1, 5, sigma.sq, phi, tau.sq)

for (n_star in sample_sizes) {
  print(n_star)
  n = n_star
  w_new = w[1:n_star]
  X_new = X[1:n_star, ]
  y_new = y[1:n_star]
  sf_new = sf[1:n_star, ]

  HGPdata <- get_HGPdata(sf_new, y_new, X_new, n.blocks, num.nb, alpha, priors)

  y_new <- y_new[(HGPdata$order)]
  X_new <- X_new[(HGPdata$order), ]

  nb <- HGPdata$nb
  ind_obs1 <- HGPdata$ind_obs1
  indb <- HGPdata$indb
  coords.D <- HGPdata$coords.D

  data1 <- data.frame(y = y_new, x = X_new[, 2])
  data1$idx <- 1:nrow(data1)
  blockNNGP.model <- HGPdata$model

  f.blockNNGP <- y ~ 1 + x + f(idx, model = blockNNGP.model)

  out <- matrix(0, ncol = 5, nrow = M)
  for (i in 1:M) {
    set.seed(i)
    print(i)

    resf <- inla(f.blockNNGP, data = as.data.frame(data1), family = "gaussian")

    # Extraction of hyperparameters
    tau.est <- inla.emarginal(function(x) 1 / x, resf$marginals.hyperpar[[1]])
    sigmasq.est <- inla.emarginal(
      function(x) exp(-x),
      resf$marginals.hyperpar[[2]]
    )
    phi.est <- inla.emarginal(
      function(x) 1 - 1 / (1 + exp(x)),
      resf$marginals.hyperpar[[3]]
    )
    w.est <- resf$summary.random$idx$mean

    summary.theta <- c(tau.est, sigmasq.est, phi.est)
    print(c("tau.est", "sigmasq.est", "phi.est"))
    print(summary.theta)

    print(resf$summary.fixed[, 1])
    out[i, ] <- c(resf$summary.fixed[, 1], sigmasq.est, phi.est, tau.est)
  }
  print("CASO GAUSSIAN")
  diagnosticar_inferencia(out, true_values, n_star, M)
}
