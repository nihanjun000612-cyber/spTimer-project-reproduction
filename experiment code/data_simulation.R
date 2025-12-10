
# ---- R/functions/data_simulation.R ----
# Simulate spatio-temporal data using GP, AR, and GPP models
library(MASS)

# 1. Gaussian Process (GP) Simulation
data.sim.GP <- function(n = 12, time_length = 100, beta = c(5, 2, 1, 0.5),
                        phi = 3, sig2eps = 1, sig2eta = 1) {
  set.seed(123)
  coords <- matrix(runif(n * 2), ncol = 2)
  dist.mat <- as.matrix(dist(coords))
  Sigma.eta <- sig2eta * exp(-phi * dist.mat)
  eta <- mvrnorm(n = 1, mu = rep(0, n), Sigma = Sigma.eta)
  time <- 1:time_length
  X <- cbind(1, sin(2 * pi * time / time_length), cos(2 * pi * time / time_length), (time / time_length)^2)
  Y <- array(0, dim = c(n, time_length))
  for (t in 1:time_length) {
    mu <- X[t, ] %*% beta + eta
    Y[, t] <- mu + rnorm(n, sd = sqrt(sig2eps))
  }
  return(list(Y = Y, coords = coords, X = X, beta = beta, sig2eta = sig2eta))
}

# 2. Autoregressive (AR) Simulation
data.sim.AR <- function(n = 12, time_length = 100, beta = c(5, 2, 1, 0.5),
                        phi = 3, sig2eps = 1, sig2eta = 1, rho = 0.7) {
  set.seed(123)
  coords <- matrix(runif(n * 2), ncol = 2)
  dist.mat <- as.matrix(dist(coords))
  Sigma.eta <- sig2eta * exp(-phi * dist.mat)
  eta <- mvrnorm(n = 1, mu = rep(0, n), Sigma = Sigma.eta)
  time <- 1:time_length
  X <- cbind(1, sin(2 * pi * time / time_length), cos(2 * pi * time / time_length), (time / time_length)^2)
  Y <- matrix(0, nrow = n, ncol = time_length)
  Y[, 1] <- X[1, ] %*% beta + eta + rnorm(n, sd = sqrt(sig2eps))
  for (t in 2:time_length) {
    mu <- X[t, ] %*% beta + rho * Y[, t-1]
    Y[, t] <- mu + rnorm(n, sd = sqrt(sig2eps))
  }
  return(list(Y = Y, coords = coords, X = X, beta = beta, rho = rho, sig2eta = sig2eta))
}

# 3. Predictive Process (GPP) Simulation
data.sim.GPP <- function(n = 12, time_length = 100, beta = c(5, 2, 1, 0.5),
                         phi = 3, sig2eps = 1, sig2eta = 1, m = 4) {
  set.seed(123)
  coords <- matrix(runif(n * 2), ncol = 2)
  knots <- matrix(runif(m * 2), ncol = 2)
  D_nn <- as.matrix(dist(coords))
  D_mm <- as.matrix(dist(knots))
  D_nm <- as.matrix(dist(rbind(coords, knots)))[1:n, (n+1):(n+m)]
  Sigma.mm <- sig2eta * exp(-phi * D_mm)
  Sigma.nm <- sig2eta * exp(-phi * D_nm)
  Sigma.nn.approx <- Sigma.nm %*% solve(Sigma.mm) %*% t(Sigma.nm)
  eta <- mvrnorm(1, mu = rep(0, n), Sigma = Sigma.nn.approx)
  time <- 1:time_length
  X <- cbind(1, sin(2 * pi * time / time_length), cos(2 * pi * time / time_length), (time / time_length)^2)
  Y <- matrix(0, nrow = n, ncol = time_length)
  for (t in 1:time_length) {
    mu <- X[t, ] %*% beta + eta
    Y[, t] <- mu + rnorm(n, sd = sqrt(sig2eps))
  }
  return(list(Y = Y, coords = coords, X = X, beta = beta, knots = knots))
}
# 先运行这个诊断
library(spTimer)
setwd('D:/python project')
source("D:/python project/data_simulation.R")
source("D:/python project/model_fitting.R")

set.seed(123)
sim_data <- data.sim.GP(sig2eps = 1, sig2eta = 1)

cat("Data dimensions:\n")
cat("  coords:", nrow(sim_data$coords), "sites\n")
cat("  Y:", nrow(sim_data$Y), "sites x", ncol(sim_data$Y), "times\n")
cat("  Total obs:", nrow(sim_data$Y) * ncol(sim_data$Y), "\n\n")

# Fit model
model_test <- fit_model("GP", data.sim = sim_data, nItr = 100, nBurn = 50)

cat("Model$fitted structure:\n")
cat("  Class:", class(model_test$fitted), "\n")
cat("  Dimensions:", dim(model_test$fitted), "\n")
cat("  Length:", length(model_test$fitted), "\n")
cat("  First few values:", head(as.vector(model_test$fitted)), "\n\n")

# Check if it's a matrix
if (is.matrix(model_test$fitted)) {
  cat("  It's a matrix with", nrow(model_test$fitted), "rows and", ncol(model_test$fitted), "cols\n")
  cat("  Needs to be converted to vector\n")
} else if (is.vector(model_test$fitted)) {
  cat("  It's already a vector\n")
}

# Expected length
obs_Y <- as.vector(t(sim_data$Y))
cat("\nExpected length:", length(obs_Y), "\n")
cat("Actual fitted length:", length(as.vector(model_test$fitted)), "\n")