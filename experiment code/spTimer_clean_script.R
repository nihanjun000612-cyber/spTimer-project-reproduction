rm(list = ls())
library(MASS)
library(spTimer)
library(ggplot2)
library(fields)
library(dplyr)


# Gaussian Process (GP) Data Simulation

data.sim.GP <- function(n = 12, time_length = 100, 
                        beta = c(5, 2, 1, 0.5), 
                        phi = 3, sig2eps = 1, sig2eta = 1) {
  set.seed(123)
  
  coords <- matrix(runif(n * 2), ncol = 2)
  dist.mat <- as.matrix(dist(coords))
  Sigma.eta <- sig2eta * exp(-phi * dist.mat)
  eta <- mvrnorm(n = 1, mu = rep(0, n), Sigma = Sigma.eta)
  
  time <- 1:time_length
  X <- cbind(1,
             sin(2 * pi * time / time_length),
             cos(2 * pi * time / time_length),
             (time / time_length)^2)
  
  Y <- matrix(0, nrow = n, ncol = time_length)
  for (tt in 1:time_length) {
    mu <- as.numeric(X[tt, ] %*% beta + eta)
    Y[, tt] <- mu + rnorm(n, sd = sqrt(sig2eps))
  }
  
  return(list(Y = Y, coords = coords, X = X, beta = beta, sig2eta = sig2eta))
}


# Autoregressive (AR) Spatio-Temporal Data Simulation

data.sim.AR <- function(n = 12, time_length = 100, 
                        beta = c(5, 2, 1, 0.5), 
                        phi = 3, sig2eps = 1, sig2eta = 1, rho = 0.7) {
  set.seed(123)
  
  coords <- matrix(runif(n * 2), ncol = 2)
  dist.mat <- as.matrix(dist(coords))
  Sigma.eta <- sig2eta * exp(-phi * dist.mat)
  eta <- mvrnorm(n = 1, mu = rep(0, n), Sigma = Sigma.eta)
  
  time <- 1:time_length
  X <- cbind(1,
             sin(2 * pi * time / time_length),
             cos(2 * pi * time / time_length),
             (time / time_length)^2)
  
  Y <- matrix(0, nrow = n, ncol = time_length)
  Y[, 1] <- as.numeric(X[1, ] %*% beta + eta + rnorm(n, sd = sqrt(sig2eps)))
  
  for (tt in 2:time_length) {
    mu <- as.numeric(X[tt, ] %*% beta) + rho * Y[, tt - 1]
    Y[, tt] <- mu + rnorm(n, sd = sqrt(sig2eps))
  }
  
  return(list(Y = Y, coords = coords, X = X, beta = beta, rho = rho, sig2eta = sig2eta))
}

#Predictive Process (GPP) Simulation

data.sim.GPP <- function(n = 12, time_length = 100, 
                         beta = c(5, 2, 1, 0.5), 
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
  X <- cbind(1,
             sin(2 * pi * time / time_length),
             cos(2 * pi * time / time_length),
             (time / time_length)^2)
  
  Y <- matrix(0, nrow = n, ncol = time_length)
  for (tt in 1:time_length) {
    mu <- as.numeric(X[tt, ] %*% beta + eta)
    Y[, tt] <- mu + rnorm(n, sd = sqrt(sig2eps))
  }
  
  return(list(Y = Y, coords = coords, X = X, beta = beta, knots = knots))
}

# Helper: Convert to Long Format
spT.to.df <- function(data.sim) {
  coords <- data.sim$coords
  Y <- data.sim$Y
  X <- data.sim$X
  time_length <- ncol(Y)
  n <- nrow(Y)
  
  df <- data.frame()
  for (i in 1:n) {
    df <- rbind(df, data.frame(
      site = i,
      time = 1:time_length,
      Y = as.numeric(Y[i, ]),
      x1 = X[, 2],
      x2 = X[, 3],
      x3 = X[, 4],
      lon = coords[i, 1],
      lat = coords[i, 2]
    ))
  }
  return(df)
}

#  Model Fitting Wrappers

fit.GP.model <- function(data.sim, nItr = 500, nBurn = 250) {
  
  coords <- data.sim$coords
  Y <- data.sim$Y
  X <- data.sim$X
  n_time <- ncol(Y)
  n_site <- nrow(Y)
  
  
  df <- data.frame()
  for (i in 1:n_site) {
    df <- rbind(df, data.frame(
      site = i,
      time = 1:n_time,
      Y = as.numeric(Y[i, ]),
      x1 = X[, 2],
      x2 = X[, 3],
      x3 = X[, 4],
      lon = coords[i, 1],
      lat = coords[i, 2]
    ))
  }
  
  df$time <- as.numeric(df$time)
  
  
  time.data.arg <- list(1, ncol(data.sim$Y))
  
  model_out <- spT.Gibbs(
    formula = Y ~ x1 + x2 + x3,
    data = df,
    model = "GP",
    coords = ~lon + lat,
    time.data = time.data.arg,     
    nItr = nItr,
    nBurn = nBurn,
    report = 50
  )
  
  return(model_out)
}

fit.AR.model <- function(data.sim, nItr = 500, nBurn = 250) {
  df <- spT.to.df(data.sim)
  spT.Gibbs(
    formula = Y ~ x1 + x2 + x3,
    data = df,
    model = "AR",
    coords = ~lon + lat,
    time.data = list(type = "temporal", time = "time"),
    nItr = nItr, nBurn = nBurn, report = 50
  )
}

fit.GPP.model <- function(data.sim, nItr = 500, nBurn = 250) {
  df <- spT.to.df(data.sim)
  spT.Gibbs(
    formula = Y ~ x1 + x2 + x3,
    data = df,
    model = "GPP",
    coords = ~lon + lat,
    time.data = list(type = "temporal", time = "time"),
    nItr = nItr, nBurn = nBurn, report = 50
  )
}