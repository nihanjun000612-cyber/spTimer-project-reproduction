
library('devtools')

# These packages will be required to run the code in this file
library("akima")
library("coda")
library("spacetime")
library("fields")
library("forecast")
library("MASS")
library("mgcv")
library("spBayes")
library("colorspace")
library("maps")
library("MBA")
library('spTimer')
library('sp')
#Simulation study #

#Functions to simulate data #

# Function to simulate data: GP model

data.sim.GP <- function(n, r = 1, T, sig2e = 0.01, sig2eta = 0.1, phi = 0.003, beta = 5) {
  # set.seed(11)
  n <- n * n  # say, sites
  Longitude <- seq(0, 1000, by = 1000/(sqrt(n) - 1))
  Latitude <- seq(0, 1000, by = 1000/(sqrt(n) - 1))
  long.lat <- expand.grid(Longitude, Latitude)
  site <- data.frame(s.index = 1:n, Longitude = long.lat[, 1], Latitude = long.lat[, 
    2])
  d <- as.matrix(dist(site[, 2:3], method = "euclidean", diag = TRUE, upper = TRUE))
  r <- r  # year
  T <- T  # day
  N <- n * r * T
  sig2e <- sig2e
  sig2eta <- sig2eta
  phi <- phi
  D1 <- exp(-phi * d)
  beta <- beta
  Ivec <- rep(1, n)
  z <- matrix(NA, r * T, n)
  o <- matrix(NA, r * T, n)
  if (length(beta) > 1) {
    x <- rep(1, N)
    for (i in 2:length(beta)) {
      x <- cbind(x, rnorm(N))
    }
    xb <- x %*% matrix(beta)
    xb <- matrix(c(xb), r * T, n)
    for (i in 1:(r * T)) {
      o[i, ] <- xb[i, ] + mvrnorm(1, rep(0, n), sig2eta * D1)
      z[i, ] <- o[i, ] + rnorm(1, 0, sqrt(sig2e))
    }
    dat1 <- matrix(NA, n * r * T, 4 + length(beta) - 1)
    dat1[, 5:(4 + length(beta) - 1)] <- x[, -1]
    dimnames(dat1)[[2]][5:(4 + length(beta) - 1)] <- paste("x", 1:(length(beta) - 
      1), sep = "")
  } else {
    for (i in 1:(r * T)) {
      o[i, ] <- beta + mvrnorm(1, rep(0, n), sig2eta * D1)
      z[i, ] <- o[i, ] + rnorm(1, 0, sqrt(sig2e))
    }
    dat1 <- matrix(NA, n * r * T, 4)
  }
  dat1[, 1] <- sort(rep(1:n, r * T))
  dat1[, 2] <- sort(rep(1:r, T))
  dat1[, 3] <- 1:T
  dat1[, 4] <- c(z)
  dimnames(dat1)[[2]][1:4] <- c("s.index", "year", "day", "y")
  dat1 <- as.data.frame(dat1)
  dat1 <- merge(dat1, site, by = c("s.index"), all.x = TRUE)
  dat1$y_no_mis <- dat1$y
  # set.seed(11)
  dat1[sample(1:dim(dat1)[[1]], round(dim(dat1)[[1]] * 0.05)), 4] <- NA  # 5% missing values to put
  dat1 <- dat1[order(dat1$s.index, dat1$year, dat1$day), ]
  dat1
}

## Function to simulate data: AR model

data.sim.AR <- function(n, r = 1, T, sig2e = 0.01, sig2eta = 0.1, phi = 0.003, beta = 5, 
  rho = 0.2, mu = 5, sig2 = 0.5) {
  # set.seed(111)
  n <- n * n  # say, sites
  Longitude <- seq(0, 1000, by = 1000/(sqrt(n) - 1))
  Latitude <- seq(0, 1000, by = 1000/(sqrt(n) - 1))
  long.lat <- expand.grid(Longitude, Latitude)
  site <- data.frame(s.index = 1:n, Longitude = long.lat[, 1], Latitude = long.lat[, 
    2])
  d <- as.matrix(dist(site[, 2:3], method = "euclidean", diag = TRUE, upper = TRUE))
  r <- r  # year
  T <- T  # day
  N <- n * r * T
  sig2e <- sig2e
  sig2eta <- sig2eta
  phi <- phi
  D1 <- exp(-phi * d)
  beta <- beta
  rho <- rho
  mu <- mu
  sig2 <- sig2
  z <- matrix(NA, T * r, n)
  o <- matrix(NA, (T + 1) * r, n)
  if (length(beta) > 1) {
    x <- rep(1, N)
    for (i in 2:length(beta)) {
      x <- cbind(x, rnorm(N))
    }
    xb <- x %*% matrix(beta)
    xb <- matrix(c(xb), r * T, n)
    for (j in 1:r) {
      o[1 + (j - 1) * T, ] <- mvrnorm(1, rep(mu, n), sig2 * D1)
      for (i in 1:T) {
        o[(i + 1) + (j - 1) * T, ] <- rho * o[i + (j - 1) * T, ] + xb[i + 
          (j - 1) * T, ] + mvrnorm(1, rep(0, n), sig2eta * D1)
        z[i + (j - 1) * T, ] <- o[(i + 1) + (j - 1) * T, ] + rnorm(1, 0, 
          sqrt(sig2e))
      }
    }
    dat1 <- matrix(NA, n * r * T, 4 + length(beta) - 1)
    dat1[, 5:(4 + length(beta) - 1)] <- x[, -1]
    dimnames(dat1)[[2]][5:(4 + length(beta) - 1)] <- paste("x", 1:(length(beta) - 
      1), sep = "")
  } else {
    Ivec <- rep(1, n)
    for (j in 1:r) {
      o[1 + (j - 1) * T, ] <- mvrnorm(1, rep(mu, n), sig2 * D1)
      for (i in 1:T) {
        o[(i + 1) + (j - 1) * T, ] <- rho * o[i + (j - 1) * T, ] + beta * 
          Ivec + mvrnorm(1, rep(0, n), sig2eta * D1)
        z[i + (j - 1) * T, ] <- o[(i + 1) + (j - 1) * T, ] + rnorm(1, 0, 
          sqrt(sig2e))
      }
    }
    dat1 <- matrix(NA, n * r * T, 4)
  }
  dat1[, 1] <- sort(rep(1:n, r * T))
  dat1[, 2] <- sort(rep(1:r, T))
  dat1[, 3] <- 1:T
  dat1[, 4] <- c(z)
  dimnames(dat1)[[2]][1:4] <- c("s.index", "year", "day", "y")
  dat1 <- as.data.frame(dat1)
  dat1 <- merge(dat1, site, by = c("s.index"), all.x = TRUE)
  dat1$y_no_mis <- dat1$y
  # set.seed(111)
  dat1[sample(1:dim(dat1)[[1]], round(dim(dat1)[[1]] * 0.05)), 4] <- NA  # 5% missing values to put
  dat1 <- dat1[order(dat1$s.index, dat1$year, dat1$day), ]
  dat1
}

# Function to simulate data: GPP based model

data.sim.GPP <- function(n, m = 10, r = 1, T, sig2e = 0.01, sig2eta = 0.1, phi = 0.003, 
  beta = 5, rho = 0.2, mu = 5, sig2 = 0.5) {
  # set.seed(33)
  n <- n * n  # say, sites
  Longitude <- seq(0, 1000, by = 1000/(sqrt(n) - 1))
  Latitude <- seq(0, 1000, by = 1000/(sqrt(n) - 1))
  long.lat <- expand.grid(Longitude, Latitude)
  knots.coords <- spT.grid.coords(Longitude = c(990.75, 9.25), Latitude = c(990.75, 
    9.25), by = c(m, m))
  site <- data.frame(s.index = 1:n, Longitude = long.lat[, 1], Latitude = long.lat[, 
    2])
  d <- as.matrix(dist(site[, 2:3], method = "euclidean", diag = TRUE, upper = TRUE))
  d2 <- as.matrix(dist(knots.coords, method = "euclidean", diag = TRUE, upper = TRUE))
  r <- r  # year
  T <- T  # day
  N <- n * r * T
  sig2e <- sig2e
  sig2eta <- sig2eta
  phi <- phi
  D1 <- exp(-phi * d)
  beta <- beta
  rho <- rho
  mu <- mu
  sig2 <- sig2
  D2 <- exp(-phi * d2)
  m <- m * m
  dd <- as.matrix(dist(rbind(as.matrix(site[, 2:3]), knots.coords), method = "euclidean", 
    diag = TRUE, upper = TRUE))
  C <- dd[1:dim(site)[[1]], (dim(site)[[1]] + 1):(dim(site)[[1]] + dim(knots.coords)[[1]])]
  A <- exp(-phi * C) %*% solve(D2)
  z <- matrix(NA, T * r, n)
  w <- matrix(NA, (T + 1) * r, m)
  if (length(beta) > 1) {
    x <- rep(1, N)
    for (i in 2:length(beta)) {
      x <- cbind(x, rnorm(N))
    }
    xb <- x %*% matrix(beta)
    xb <- matrix(c(xb), r * T, n)
    for (j in 1:r) {
      w[1 + (j - 1) * T, ] <- mvrnorm(1, rep(0, m), sig2 * D2)
      for (i in 1:T) {
        w[(i + 1) + (j - 1) * T, ] <- rho * w[i + (j - 1) * T, ] + mvrnorm(1, 
          rep(0, m), sig2eta * D2)
      }
    }
    for (j in 1:r) {
      for (i in 1:T) {
        e <- rnorm(1, 0, sqrt(sig2e))
        z[i + (j - 1) * T, ] <- A %*% w[(i + 1) + (j - 1) * T, ] + xb[i + 
          (j - 1) * T, ] + e
      }
    }
    dat1 <- matrix(NA, n * r * T, 4 + length(beta) - 1)
    dat1[, 5:(4 + length(beta) - 1)] <- x[, -1]
    dimnames(dat1)[[2]][5:(4 + length(beta) - 1)] <- paste("x", 1:(length(beta) - 
      1), sep = "")
  } else {
    Ivec <- rep(1, n)
    for (j in 1:r) {
      w[1 + (j - 1) * T, ] <- mvrnorm(1, rep(0, m), sig2 * D2)
      for (i in 1:T) {
        w[(i + 1) + (j - 1) * T, ] <- rho * w[i + (j - 1) * T, ] + mvrnorm(1, 
          rep(0, m), sig2eta * D2)
      }
    }
    for (j in 1:r) {
      for (i in 1:T) {
        e <- rnorm(1, 0, sqrt(sig2e))
        z[i + (j - 1) * T, ] <- A %*% w[(i + 1) + (j - 1) * T, ] + beta * 
          Ivec + e
      }
    }
    dat1 <- matrix(NA, n * r * T, 4)
  }
  dat1[, 1] <- sort(rep(1:n, r * T))
  dat1[, 2] <- sort(rep(1:r, T))
  dat1[, 3] <- 1:T
  dat1[, 4] <- c(z)
  dimnames(dat1)[[2]][1:4] <- c("s.index", "year", "day", "y")
  dat1 <- as.data.frame(dat1)
  dat1 <- merge(dat1, site, by = c("s.index"), all.x = TRUE)
  dat1$y_no_mis <- dat1$y
  # set.seed(33)
  dat1[sample(1:dim(dat1)[[1]], round(dim(dat1)[[1]] * 0.05)), 4] <- NA  # 5% missing values to put
  dat1 <- dat1[order(dat1$s.index, dat1$year, dat1$day), ]
  dat1
}


#Figure 1(a) and 1(b)

# spatial domain: Figure 1(a)
n <- 12 * 12  # say, sites
Longitude <- seq(0, 1000, by = 1000/(12 - 1))
Latitude <- seq(0, 1000, by = 1000/(12 - 1))
long.lat <- expand.grid(Longitude, Latitude)
plot(long.lat, xlab = "Longitude", ylab = "Latitude", pch = 3, col = 4)

# spatial domain: Figure 1(b)
n <- 55 * 55  # say, sites
Longitude <- seq(0, 1000, by = 1000/(55 - 1))
Latitude <- seq(0, 1000, by = 1000/(55 - 1))
long.lat <- expand.grid(Longitude, Latitude)
spT.grid.coords <- function(Longitude, Latitude, by) {
  long_seq <- seq(from = min(Longitude), to = max(Longitude), length.out = by[1])
  lat_seq <- seq(from = min(Latitude), to = max(Latitude), length.out = by[2])
  coords <- expand.grid(Longitude = long_seq, Latitude = lat_seq)
  return(coords)
}
knots.coords <- spT.grid.coords(Longitude = c(990, 10), Latitude = c(990, 10), by = c(10, 10))
plot(long.lat, xlab = "Longitude", ylab = "Latitude", pch = 3, col = 4, cex = 0.6)
points(knots.coords, pch = 19, col = 2)


# Replications 

# To obtain results quickly I have changed the following numbers

replic <- 2  # number of replications of the datasets, used as 25 in the paper
nItr <- 200  # number of MCMC samples for each model, used as 5000 in the paper
nBurn <- 100  # number of burn-in from the MCMC samples, used as 1000 in the paper 

# For GP model

paraGP <- NULL
for (i in 1:replic) {
  set.seed(round(rnorm(1, mean = i, sd = 100)))
  dat <- data.sim.GP(n = 12, T = 365, beta = c(5, 2, 1, 0.5), sig2eta = runif(1, 
    0, 1))
  out <- spT.Gibbs(formula = y ~ x1 + x2 + x3, data = dat, model = "GP", coords = ~Longitude + 
    Latitude, distance.method = "euclidean", nItr = nItr, nBurn = nBurn, report = 1, 
    spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.9))
  paraGP <- rbind(paraGP, as.mcmc(out)[, 1:(dim(as.mcmc(out))[[2]] - 3)])
}

# For AR model

paraAR <- NULL
for (i in 1:replic) {
  set.seed(round(rnorm(1, mean = i, sd = 100)))
  dat <- data.sim.AR(n = 12, T = 365, beta = c(5, 2, 1, 0.5), sig2eta = runif(1, 
    0, 1))
  out <- spT.Gibbs(formula = y ~ x1 + x2 + x3, data = dat, model = "AR", coords = ~Longitude + 
    Latitude, distance.method = "euclidean", nItr = nItr, nBurn = nBurn, report = 1, 
    spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.2))
  paraAR <- rbind(paraAR, as.mcmc(out)[, 1:(dim(as.mcmc(out))[[2]] - 3)])
}

# For GPP model

spT.grid.coords <- function(Longitude, Latitude, by) {
  long_seq <- seq(from = min(Longitude), to = max(Longitude), length.out = by[1])
  lat_seq  <- seq(from = min(Latitude), to = max(Latitude), length.out = by[2])
  coords   <- expand.grid(Longitude = long_seq, Latitude = lat_seq)
  return(coords)
}

paraGPP <- NULL
for (i in 1:replic) {
  set.seed(round(rnorm(1, mean = i, sd = 100)))
  
  dat <- data.sim.GPP(
    n = 55,
    T = 365,
    beta = c(5, 2, 1, 0.5),
    sig2eta = runif(1, 0, 1)
  )
  
  knots.coords <- as.matrix(
    spT.grid.coords(
      Longitude = c(990.75, 9.25),
      Latitude = c(990.75, 9.25),
      by = c(10, 10)
    )
  )
  
  out <- spT.Gibbs(
    formula = y ~ x1 + x2 + x3,
    data = dat,
    model = "GPP",
    coords = ~Longitude + Latitude,
    knots.coords = knots.coords,
    distance.method = "euclidean",
    nItr = nItr,
    nBurn = nBurn,
    report = 1,
    tol.dist = 1,
    spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.08)
  )
  
  paraGPP <- rbind(paraGPP, as.mcmc(out)[, 1:(dim(as.mcmc(out))[[2]] - 3)])
}

# Table 2 

# summary of parameter estimates: Table 2

t(apply(paraGP, 2, quantile, prob = c(0.025, 0.5, 0.975)))
t(apply(paraAR, 2, quantile, prob = c(0.025, 0.5, 0.975)))
t(apply(paraGPP, 2, quantile, prob = c(0.025, 0.5, 0.975)))

# Figure 2(a) and 2(b) #

# Surface plot: Figure 2(a) and 2(b)
spT.grid.coords <- function(Longitude, Latitude, by) {
  long_seq <- seq(from = min(Longitude), to = max(Longitude), length.out = by[1])
  lat_seq  <- seq(from = min(Latitude),  to = max(Latitude),  length.out = by[2])
  coords   <- expand.grid(Longitude = long_seq, Latitude = lat_seq)
  as.matrix(coords)  # crucial!
}
set.seed(round(rnorm(1, mean = i, sd = 100)))
n <- 18
T <- 30
dat <- data.sim.GPP(n = n, T = T, beta = c(5, 2, 1, 0.5))

# returns a 49 × 2 numeric matrix
knots.coords <- spT.grid.coords(
  Longitude = c(990.75, 9.25),
  Latitude  = c(990.75, 9.25),
  by = c(7, 7)
)

out <- spT.Gibbs(
  formula = y ~ x1 + x2 + x3,
  data = dat,
  model = "GPP",
  coords = ~Longitude + Latitude,
  knots.coords = knots.coords,          # now a matrix
  distance.method = "euclidean",
  nItr = nItr, nBurn = nBurn, report = 1,
  spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.08)
)

# S3 class code for plot.spT and contour.spT using R package akima
plot.spT <- function(x, residuals = FALSE, surface = NULL, time = c(1), a3d = FALSE, 
  points = FALSE, title = TRUE, ...) {
  if (is.null(surface) & a3d == FALSE) {
    if (as.logical(residuals) == FALSE) {
      tmp <- as.mcmc(x)
      plot(tmp, ...)
    } else {
      plot(x$fitted[, 1], residuals(x), ylab = "Residuals", xlab = "Fitted values")
      abline(h = 0, lty = 2)
      title("Residuals vs Fitted")
      par(ask = TRUE)
      qqnorm(residuals(x))
      qqline(residuals(x), lty = 2)
    }
  } else {
    if (is.null(surface)) {
      stop("\n# Error: surface should be defined as 'Mean' or 'SD'. \n")
    }
    if (!surface %in% c("Mean", "SD")) {
      stop("\n# Error: surface only takes 'Mean' or 'SD'. \n")
    }
    z <- array(fitted(x)[, paste(surface)], dim = c(x$T * x$r, x$n))
    xyz <- cbind(x$coords, c(z[time, ]))
    xyz <- interp(x = xyz[, 1], y = xyz[, 2], z = xyz[, 3], xo = seq(min(xyz[, 
      1]), max(xyz[, 1]), length = 150), yo = seq(min(xyz[, 2]), max(xyz[, 
      2]), length = 150))
    if (a3d == TRUE) {
      res <- persp(x = xyz$x, y = xyz$y, z = xyz$z, xlab = "x", ylab = "y", 
        zlab = "z", ...)
    } else {
      image.plot(xyz, ...)
      if (points != FALSE) {
        points(x$coords, pch = 16, cex = 0.8)
      }
    }
    if (title == TRUE) {
      title(main = paste("Time point: (t=", time, ")", sep = ""))
    }
  }
}
contour.spT <- function(x, surface = "Mean", time = c(1), ...) {
  z <- array(fitted(x)[, paste(surface)], dim = c(x$T * x$r, x$n))
  xyz <- cbind(x$coords, c(z[time, ]))
  xyz <- interp(x = xyz[, 1], y = xyz[, 2], z = xyz[, 3], xo = seq(min(xyz[, 1]), 
    max(xyz[, 1]), length = 150), yo = seq(min(xyz[, 2]), max(xyz[, 2]), length = 150), 
    linear = TRUE, extrap = FALSE, duplicate = "error", dupfun = NULL, ncp = NULL)
  contour(xyz, ...)
}

# For Mean
time <- 5
plot(out, surface = "Mean", time = time, col = terrain_hcl(12), legend.shrink = 0.5, 
  legend.width = 0.8, horizontal = TRUE, title = FALSE)
# contour(out,add=TRUE,lty=2)

# For surface error plot
err <- array(dat$y_no_mis - out$fitted[, 1], dim = c(T, n * n))
xyz <- cbind(unique(dat[, c("Longitude", "Latitude")]), c(err[time, ]))
xyz <- interp(x = xyz[, 1], y = xyz[, 2], z = xyz[, 3], xo = seq(min(xyz[, 1]), max(xyz[, 
  1]), length = 150), yo = seq(min(xyz[, 2]), max(xyz[, 2]), length = 150), linear = TRUE, 
  extrap = FALSE, duplicate = "error", dupfun = NULL)
image.plot(xyz, col = diverge_hcl(7, h = c(246, 40), c = 96, l = c(65, 90)), legend.shrink = 0.5, 
  legend.width = 0.8, horizontal = TRUE)


# Figure 3(a) and 3(b) 

## Time-series plot: Figure 3(a)
z1 <- array(dat$y, dim = c(T, n * n))
z2 <- array(dat$y_no_mis, dim = c(T, n * n))
z3 <- array(c(out$fitted[, 1]), dim = c(T, n * n))
s.index <- 1
zz <- cbind(1:length(c(z3[, s.index])), c(z3[, s.index]), c(z1[, s.index]), c(z2[, 
  s.index]))
plot(zz[, 4], type = "o", xlab = "Time series", pch = 16, ylab = "y", axes = FALSE, 
  cex = 0.9, ylim = c(0, 10), lty = 2, lwd = 1.5)
points(zz[, 2], pch = 12, col = 2, type = "o", lty = 1, lwd = 1)
# Please note that if
zzz <- zz[is.na(zz[, 3]), ]
if (is.null(dim(zzz))) {
  points(zzz[1], zzz[4], pch = 1, cex = 2, col = 4)
}
if (!is.null(dim(zzz))) {
  points(zzz[, 1], zzz[, 4], pch = 1, cex = 2, col = 4)
}
axis(2)
axis(1, 1:T, labels = 1:T)
legend("topleft", pch = c(16, 12, 1), col = c(1, 2, 4), lty = c(2, 1, NA), bty = "n", 
  cex = 1, legend = c("True values", "Fitted values", "Missing values"))

#Time-series plot: Figure 3(b)
library(ggplot2)
library(reshape2)
T <- nrow(z3)
s.indexx <- c(1, 2, 3, 6)
data_list <- list()

for (i in seq_along(s.indexx)) {
  s.index <- s.indexx[i]
  df <- data.frame(
    Time = 1:T,
    Series = paste0("Series ", i),
    Error = z3[, s.index] - z2[, s.index]
  )
  data_list[[i]] <- df
}
df_all <- do.call(rbind, data_list)


ggplot(df_all, aes(Time, Error, color = Series, shape = Series, linetype = Series)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Time series", y = "Error") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Figure 4 

# Sensitivity analysis for signal-to-noise ratio (SNR)

# SNR = sig2eta/sig2e = 0.01/0.01 = 1
ptm <- proc.time()
para1 <- NULL
stn1 <- NULL
?spT.subset
args(spT.subset)
packageVersion("spTimer")
args(spT.subset)
getAnywhere(spT.subset)
for (i in 1:replic) {
  dat <- data.sim.GP(n = 12, T = 365, beta = c(5), sig2eta = 0.01)
  s1 <- sample(1:(12 * 12), 15)
  fit <- spT.subset(data = dat, var.name = "s.index", s= s1, reverse = TRUE)  # for model fitting
  # model fitting with spTimer
  out <- spT.Gibbs(formula = y ~ 1, data = fit, coords = ~Longitude + Latitude, 
    nItr = nItr, nBurn = nBurn, report = 1, distance.method = "euclidean", spatial.decay = spT.decay(distribution = Gamm(2, 
      1), tuning = 0.2))
  para1 <- rbind(para1, as.mcmc(out)[, 1])
  stn1 <- c(stn1, as.mcmc(out)[, 3]/as.mcmc(out)[, 2])
}
paraGP1 <- quantile(c(para1), probs = c(0.025, 0.5, 0.975))
rm(out)
proc.time() - ptm
base::ls(pattern = "stn")
# SNR = sig2eta/sig2e = 0.1/0.01 = 10
ptm <- proc.time()
para10 <- NULL
stn10 <- NULL
for (i in 1:replic) {
  dat <- data.sim.GP(n = 12, T = 365, beta = c(5), sig2eta = 0.1)
  s <- sample(1:(12 * 12), 15)
  fit <- spT.subset(data = dat, var.name = "s.index", s = s1, reverse = TRUE)  # for model fitting
  # model fitting with spTimer
  out <- spT.Gibbs(formula = y ~ 1, data = fit, coords = ~Longitude + Latitude, 
    nItr = nItr, nBurn = nBurn, report = 1, distance.method = "euclidean", spatial.decay = spT.decay(distribution = Gamm(2, 
      1), tuning = 0.1))
  para10 <- rbind(para10, as.mcmc(out)[, 1])
  stn10 <- c(stn10, as.mcmc(out)[, 3]/as.mcmc(out)[, 2])
}
paraGP10 <- quantile(c(para10), probs = c(0.025, 0.5, 0.975))
rm(out)
proc.time() - ptm

# SNR = sig2eta/sig2e = 0.15/0.01 = 15

ptm <- proc.time()
para15 <- NULL
stn15 <- NULL
for (i in 1:replic) {
  dat <- data.sim.GP(n = 12, T = 365, beta = c(5), sig2eta = 0.25)
  s <- sample(1:(12 * 12), 15)
  fit <- spT.subset(data = dat, var.name = "s.index", s = s1, reverse = TRUE)  # for model fitting
  # model fitting with spTimer
  out <- spT.Gibbs(formula = y ~ 1, data = fit, coords = ~Longitude + Latitude, 
    nItr = nItr, nBurn = nBurn, report = 1, distance.method = "euclidean", spatial.decay = spT.decay(distribution = Gamm(2, 
      1), tuning = 0.1))
  para15 <- rbind(para15, as.mcmc(out)[, 1])
  stn15 <- c(stn15, as.mcmc(out)[, 3]/as.mcmc(out)[, 2])
}
paraGP15 <- quantile(c(para15), probs = c(0.025, 0.5, 0.975))
rm(out)
proc.time() - ptm

# SNR density plot: Code for Figure 4
stn1  <- as.numeric(unlist(stn1))
stn10 <- as.numeric(unlist(stn10))
stn15 <- as.numeric(unlist(stn15))
if (any(is.na(stn1)) || any(is.na(stn10)) || any(is.na(stn15))) {
  warning("Some SNR vectors contain missing values; density() will ignore them.")
}
lengths(list(stn1 = stn1, stn10 = stn10, stn15 = stn15))
sapply(list(stn1 = stn1, stn10 = stn10, stn15 = stn15), function(x) length(unique(x)))
summary(list(stn1 = stn1, stn10 = stn10, stn15 = stn15))
plot(density(c(stn1)), type = "n", xlim = c(0, 25), main = "Signal-to-noise ratio (SNR)")
lines(density(c(stn1)), lty = 2, col = 2, lwd = 2)
abline(v = 1, lty = 3)
lines(density(c(stn10)), lty = 4, col = 4, lwd = 2)
abline(v = 10, lty = 3)
lines(density(c(stn15)), lty = 6, col = 6, lwd = 2)
abline(v = 15, lty = 3)
legend("topright", lty = c(2, 4, 6, 3), col = c(2, 4, 6, 1), bty = "n", cex = 0.8, 
  lwd = 1.5, legend = c("Distribution of SNR for 1", "Distribution of SNR for 10", 
    "Distribution of SNR for 15", "True SNR values"))
text(x = 1, y = 0.5, label = "SNR = 1", cex = 0.8)
text(x = 10, y = 0.5, label = "SNR = 10", cex = 0.8)
text(x = 15, y = 0.5, label = "SNR = 15", cex = 0.8)


# Table 3 

# SNR: Table 3

paraGP1
paraGP10
paraGP15

#Comparison study 

#Function to simulate data GP model without grid points

data.sim.GP.nogrid <- function(nn, r = 1, T, sig2e = 0.01, sig2eta = 0.1, phi = 0.003, 
  beta = 5) {
  # set.seed(11)
  n <- nn  # say, sites
  Longitude <- sample(0:1000, n)
  Latitude <- sample(0:1000, n)
  long.lat <- cbind(Longitude, Latitude)
  site <- data.frame(s.index = 1:n, Longitude = long.lat[, 1], Latitude = long.lat[, 
    2])
  d <- as.matrix(dist(site[, 2:3], method = "euclidean", diag = TRUE, upper = TRUE))
  r <- r  # year
  T <- T  # day
  N <- n * r * T
  sig2e <- sig2e
  sig2eta <- sig2eta
  phi <- phi
  D1 <- exp(-phi * d)
  beta <- beta
  Ivec <- rep(1, n)
  z <- matrix(NA, r * T, n)
  o <- matrix(NA, r * T, n)
  if (length(beta) > 1) {
    x <- rep(1, N)
    for (i in 2:length(beta)) {
      x <- cbind(x, rnorm(N))
    }
    xb <- x %*% matrix(beta)
    xb <- matrix(c(xb), r * T, n)
    for (i in 1:(r * T)) {
      o[i, ] <- xb[i, ] + mvrnorm(1, rep(0, n), sig2eta * D1)
      z[i, ] <- o[i, ] + rnorm(1, 0, sqrt(sig2e))
    }
    dat1 <- matrix(NA, n * r * T, 4 + length(beta) - 1)
    dat1[, 5:(4 + length(beta) - 1)] <- x[, -1]
    dimnames(dat1)[[2]][5:(4 + length(beta) - 1)] <- paste("x", 1:(length(beta) - 
      1), sep = "")
  } else {
    for (i in 1:(r * T)) {
      o[i, ] <- beta + mvrnorm(1, rep(0, n), sig2eta * D1)
      z[i, ] <- o[i, ] + rnorm(1, 0, sqrt(sig2e))
    }
    dat1 <- matrix(NA, n * r * T, 4)
  }
  dat1[, 1] <- sort(rep(1:n, r * T))
  dat1[, 2] <- sort(rep(1:r, T))
  dat1[, 3] <- 1:T
  dat1[, 4] <- c(z)
  dimnames(dat1)[[2]][1:4] <- c("s.index", "year", "day", "y")
  dat1 <- as.data.frame(dat1)
  dat1 <- merge(dat1, site, by = c("s.index"), all.x = TRUE)
  dat1$y_no_mis <- dat1$y
  # set.seed(11)
  dat1[sample(1:dim(dat1)[[1]], round(dim(dat1)[[1]] * 0.05)), 4] <- NA  # 5% missing values to put
  dat1 <- dat1[order(dat1$s.index, dat1$year, dat1$day), ]
  dat1
}


## creating spBayes fnc to run the programme in the paper

run.spBayes <- function(nItr, data) {
  fit <- data
  coords <- unique(cbind(fit$Longitude, fit$Latitude))
  # spBayes cannot handle the missing values in y automatically for multivariate
  # case
  y <- matrix(fit$y, time, dim(fit)[[1]]/time)
  y <- cbind(c(y), rep(apply(y, 1, mean, na.rm = T), dim(fit)[[1]]/time))
  y[is.na(y[, 1]), 1] <- y[is.na(y[, 1]), 2]
  y[is.na(y[, 1]), 1] <- median(y[, 2], na.rm = TRUE)
  y <- matrix(y[, 1], time, dim(fit)[[1]]/time)
  x1 <- matrix(1, time, dim(fit)[[1]]/time)
  f <- NULL
  # need to supply equation for each day
  for (i in 1:time) {
    f[[i]] <- as.formula(paste("y[", i, ",]~x1[", i, ",]-1", sep = ""))
  }
  if (time > 1) {
    # Call spMvLM for more than one time points
    q <- time
    A.starting <- diag(1, q)[lower.tri(diag(1, q), TRUE)]
    n.samples <- nItr
    starting <- list(phi = rep(3/0.5, q), A = A.starting, Psi = rep(1, q))
    tuning <- list(phi = rep(50, q), A = rep(1e-04, length(A.starting)), Psi = rep(50, 
      q))
    priors <- list("beta.Flat", phi.Unif = list(rep(3/0.75, q), rep(3/0.25, q)), 
      K.IW = list(q + 1, diag(0.1, q)), Psi.ig = list(rep(2, q), rep(0.1, q)))
    out2 <- spMvLM(f, coords = coords, starting = starting, tuning = tuning, 
      priors = priors, n.samples = nItr, cov.model = "exponential", n.report = nItr)
    out2
  } else {
    # Call spLM for univariate model
    starting <- list(phi = 3/0.5, sigma.sq = 50, tau.sq = 1)
    tuning <- list(phi = 0.1, sigma.sq = 0.1, tau.sq = 0.1)
    priors <- list("beta.Flat", phi.Unif = c(3/1, 3/0.1), sigma.sq.IG = c(2, 
      5), tau.sq.IG = c(2, 0.01))
    out2 <- spLM(f[[1]], coords = coords, starting = starting, tuning = tuning, 
      priors = priors, n.samples = nItr, cov.model = "exponential", n.report = nItr)
    out2
  }
}


# Figures 5(a) and 5(b) 

# create data set in Grids: Figure 5(a)
n <- 7
nn <- n * n  # say, sites
Longitude <- seq(0, 1000, by = 1000/(sqrt(nn) - 1))
Latitude <- seq(0, 1000, by = 1000/(sqrt(nn) - 1))
long.lat <- expand.grid(Longitude, Latitude)
plot(long.lat, xlab = "Longitude", ylab = "Latitude", pch = 3, col = 4)

## create data set randomly without Grids: Figure 5(b)
dat <- data.sim.GP.nogrid(nn = n * n, T = time, beta = c(5))
plot(unique(dat[, c("Longitude", "Latitude")]), xlab = "Longitude", ylab = "Latitude", 
  pch = 3, col = 4)
rm(dat)



# Comparison study
n <- 7
time <- 5  # takes 5, 10, 20, 60
dat <- data.sim.GP(n = n, T = time, beta = c(5))

# spTimer
out1 <- spT.Gibbs(formula = y ~ 1, data = dat, coords = ~Longitude + Latitude, nItr = nItr, 
                  nBurn = nBurn, distance.method = "euclidean", spatial.decay = spT.decay(distribution = Gamm(2, 
                                                                                                              1), tuning = 0.1))


#spBayes
start.time <- proc.time()[3]
str(dat)
names(dat)
nrow(dat)
sum(complete.cases(dat)) 
run.spBayes <- function(nItr, data, time) {
  fit <- data
  coords <- unique(cbind(fit$Longitude, fit$Latitude))
  
  # Construct Y matrix (time × sites)
  nSites <- nrow(coords)
  y <- matrix(fit$y, nrow = time, ncol = nSites)
  
  # Replace missing values with mean
  for (i in 1:ncol(y)) {
    if (any(is.na(y[, i]))) {
      y[is.na(y[, i]), i] <- mean(y[, i], na.rm = TRUE)
    }
  }
  
  # Design matrix
  x1 <- matrix(1, nrow = time, ncol = nSites)
  
  # Formula list for spMvLM
  f <- lapply(1:time, function(i) as.formula(paste("y[", i, ",] ~ x1[", i, ",] - 1", sep = "")))
  
  if (time > 1) {
    ## Multivariate case (spMvLM)
    q <- time
    A.starting <- diag(1, q)[lower.tri(diag(1, q), TRUE)]
    starting <- list(phi = rep(3/0.5, q), A = A.starting, Psi = rep(1, q))
    tuning <- list(phi = rep(50, q), A = rep(1e-04, length(A.starting)), Psi = rep(50, q))
    priors <- list(
      "beta.Flat",
      phi.Unif = list(rep(3/0.75, q), rep(3/0.25, q)),
      K.IW = list(q + 1, diag(0.1, q)),
      Psi.ig = list(rep(2, q), rep(0.1, q))
    )
    
    out2 <- spMvLM(f, coords = coords, starting = starting, tuning = tuning,
                   priors = priors, n.samples = nItr,
                   cov.model = "exponential", n.report = nItr)
    
  } else {
    ## Univariate case (spLM)
    starting <- list(phi = 3/0.5, sigma.sq = 50, tau.sq = 1)
    tuning   <- list(phi = 0.1, sigma.sq = 0.1, tau.sq = 0.1)
    priors   <- list("beta.Flat", phi.Unif = c(3/1, 3/0.1),
                     sigma.sq.IG = c(2, 5), tau.sq.IG = c(2, 0.01))
    
    out2 <- spLM(f[[1]], coords = coords, starting = starting, tuning = tuning,
                 priors = priors, n.samples = nItr,
                 cov.model = "exponential", n.report = nItr)
  }
  
  return(out2)
}
out2 <- run.spBayes(nItr = nItr, data = dat, time = time)
out2 <- spRecover(out2, start = nBurn, verbose = FALSE)
end.time <- proc.time()[3]
comp.time <- end.time - start.time
comp.time

# GAM using mgcv
out3 <- gam(formula = y ~ s(Longitude, Latitude), data = dat)


# Figure 6(a) and 6(b) 

# density plots for time=5 and time=10 change time points as 5 and 10
tmp <- out2$p.beta.recover.samples
plot(density(tmp), col = 4, lty = 4, main = paste("Data with time point = ", time, 
  sep = ""), xlab = "", ylim = c(0, 6), xlim = c(0, 10))
lines(density(c(out1$betap)), col = 2, lty = 2)
abline(v = coef(out3)[[1]], lty = 3)
abline(v = 5, lty = 1)
legend("topleft", col = c(1, 1, 4, 2), lty = c(1, 3, 4, 2), bty = "n", cex = 0.8, 
  legend = c("True", "GAM", "spBayes", "spTimer"))

 

# Predictions validations Grid data
n <- 7
time <- 5  # takes 5, and 60

# set.seed(22)
dat <- data.sim.GP(n = n, T = time, beta = c(5))

# Model fitting set.seed(11)
s <- sample(1:(n * n), n * n * 0.2)
val <- spT.subset(data = dat, var.name = "s.index", s = s, reverse = FALSE)  # for model validation
fit <- spT.subset(data = dat, var.name = "s.index", s = s, reverse = TRUE)  # for model fitting

# spTimer set.seed(11)
out1 <- spT.Gibbs(formula = y ~ 1, data = fit, coords = ~Longitude + Latitude, nItr = nItr, 
  nBurn = nBurn, distance.method = "euclidean", spatial.decay = spT.decay(distribution = Gamm(2, 
    1), tuning = 0.1))
pred1 <- predict(out1, newdata = val, newcoords = ~Longitude + Latitude)

# spBayes
out2 <- run.spBayes(nItr = nItr, data = fit, time=time)
out2 <- spRecover(out2, verbose = FALSE)
pred.coords <- unique(cbind(val$Longitude, val$Latitude))
x1 <- matrix(1, length(s), dim(val)[[1]]/length(s))
xx <- NULL
for (i in 1:dim(x1)[[2]]) {
  xx[[i]] <- as.matrix(x1[, i])
}
x.mv <- mkMvX(xx)
pred2 <- spPredict(out2, pred.covars = x.mv, pred.coords = pred.coords)

# GAM
out3 <- gam(formula = y ~ s(Longitude, Latitude), data = fit)
pred3 <- predict(out3, newdata = val)

# validation statistics
val.spTimer <- spT.validation(val$y, pred1$Median)
val.spBayes <- spT.validation(val$y, apply(pred2[[1]], 1, median))
val.gam <- spT.validation(val$y, c(pred3))

val.spTimer
val.spBayes
val.gam

rm(out1)
rm(out2)
rm(out3)

# Non-grid data
n <- 7
time <- 5  # takes 5, and 60

# set.seed(22)
dat <- data.sim.GP.nogrid(nn = n * n, T = time, beta = c(5))

# Model fitting set.seed(11)
s <- sample(1:(n * n), n * n * 0.2)
val <- spT.subset(data = dat, var.name = "s.index", s = s, reverse = FALSE)  # for model validation
fit <- spT.subset(data = dat, var.name = "s.index", s = s, reverse = TRUE)  # for model fitting

# spTimer set.seed(11)
out1 <- spT.Gibbs(formula = y ~ 1, data = fit, coords = ~Longitude + Latitude, nItr = nItr, 
  nBurn = nBurn, distance.method = "euclidean", spatial.decay = spT.decay(distribution = Gamm(2, 
    1), tuning = 1e-04))
pred1 <- predict(out1, newdata = val, newcoords = ~Longitude + Latitude)
dim(fit)
# FALSE
# spBayes
fit <- fit[complete.cases(fit), ]
out2 <- run.spBayes(nItr = nItr, data = fit, time=time)
out2 <- spRecover(out2, start = nBurn, verbose = FALSE)
pred.coords <- unique(cbind(val$Longitude, val$Latitude))
x1 <- matrix(1, length(s), dim(val)[[1]]/length(s))
xx <- NULL
for (i in 1:dim(x1)[[2]]) {
  xx[[i]] <- as.matrix(x1[, i])
}
x.mv <- mkMvX(xx)
pred2 <- spPredict(out2, pred.covars = x.mv, pred.coords = pred.coords)


# GAM
out3 <- gam(formula = y ~ s(Longitude, Latitude), data = fit)
pred3 <- predict(out3, newdata = val)

# validation statistics
val.spTimer <- spT.validation(pred1$Median, val$y)
val.spBayes <- spT.validation(apply(pred2[[1]], 1, median), val$y)
val.gam <- spT.validation(c(pred3), val$y)

val.spTimer
val.spBayes
val.gam

rm(out1)
rm(out2)
rm(out3)

#real data example


# Read data
data("NYdata")
s <- c(8, 11, 12, 14, 18, 21, 24, 28)
DataFit <- spT.subset(data = NYdata, var.name = c("s.index"), s = s, reverse = TRUE)
DataFit <- subset(DataFit, with(DataFit, !(Day %in% c(30, 31) & Month == 8)))
DataValPred <- spT.subset(data = NYdata, var.name = c("s.index"), s = s)
DataValPred <- subset(DataValPred, with(DataValPred, !(Day %in% c(30, 31) & Month == 
  8)))


# Figure 7

# Load required package
library(maps)

# Extract coordinate matrices
coords <- as.matrix(unique(cbind(DataFit[, 2:3])))
pred.coords <- as.matrix(unique(cbind(DataValPred[, 2:3])))

# ---- Export the plot directly to a file in D:\python project ----
png(filename = "D:/python project/Figure7_NewYork_Map.png",
    width = 1200, height = 900, res = 150)

# Draw New York State map
map(database = "state", regions = "new york")

# Add points for fitted and validation sites
points(coords, pch = 19, col = 3)    # Fitted sites (solid circles)
points(coords, pch = 1, col = 1)     # Optional outline
points(pred.coords, pch = 3, col = 4) # Validation sites (crosses)

# Add legend
legend(
  x = -77.5, y = 41.5,
  col = c(3, 4), pch = c(19, 3),
  cex = 0.8,
  legend = c("Fitted sites", "Validation sites")
)

# Close and save the file
dev.off()


# Fit GP model
set.seed(11)
post.gp <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = DataFit, model = "GP", 
  coords = ~Longitude + Latitude, scale.transform = "SQRT", spatial.decay = spT.decay(distribution = Gamm(2, 
    1), tuning = 0.1))
print(post.gp)
summary(post.gp)
#plot(post.gp)

# Spatial prediction for the GP model
set.seed(11)
pred.gp <- predict(post.gp, newdata = DataValPred, newcoords = ~Longitude + Latitude)
print(pred.gp)
names(pred.gp)
# validation criteria
spT.validation(DataValPred$o8hrmax, c(pred.gp$Median))

dev.off()

# Figures 8 (a) -- (d)
data("NYgrid")
set.seed(11)
post.gp2 <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = NYdata, model = "GP", 
  coords = ~Longitude + Latitude, scale.transform = "SQRT", nItr = 15000, nBurn = 0, 
  spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.1))
while(!is.null(dev.list())) dev.off()
dev.new(width = 8, height = 8)  
par(mfrow = c(2, 2))
plot(post.gp2$betap[1, ], type = "l", main = "Intercept", xlab = "Iterations", ylab = "", 
  ylim = c(-1, 6))
plot(post.gp2$betap[2, ], type = "l", main = "cMAXTMP", xlab = "Iterations", ylab = "", 
  ylim = c(0.06, 0.25))
plot(post.gp2$betap[3, ], type = "l", main = "WDSP", xlab = "Iterations", ylab = "", 
  ylim = c(-0.01, 0.25))
plot(post.gp2$betap[4, ], type = "l", main = "RH", xlab = "Iterations", ylab = "", 
  ylim = c(-0.4, 0.33))
par(mfrow = c(1, 1))

# Figures 9 (a) and (b) Predict on grids
set.seed(11)
post.gp2 <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = NYdata, model = "GP", 
  coords = ~Longitude + Latitude, scale.transform = "SQRT", spatial.decay = spT.decay(distribution = Gamm(2, 
    1), tuning = 0.1))
set.seed(11)
grid.pred <- predict(post.gp2, newdata = NYgrid, newcoords = ~Longitude + Latitude)

# predictive plots

# this function is used to delete values outside NY
fnc.delete.map.XYZ <- function(xyz) {
  x <- xyz$x
  y <- xyz$y
  z <- xyz$z
  xy <- expand.grid(x, y)
  eus <- (map.where(database = "state", x = xy[, 1], y = xy[, 2]))
  dummy <- rep(0, length(xy[, 1]))
  eastUS <- NULL
  eastUS <- data.frame(lon = xy[, 1], lat = xy[, 2], state = eus, dummy = dummy)
  eastUS[!is.na(eastUS[, 3]), 4] <- 1
  eastUS[eastUS[, 3] == "pennsylvania" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "new jersey" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "connecticut" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "massachusetts:main" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "new hampshire" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "vermont" & !is.na(eastUS[, 3]), 4] <- 0
  a <- eastUS[, 4]
  z <- as.vector(xyz$z)
  z[!a] <- NA
  z <- matrix(z, nrow = length(xyz$x))
  xyz$z <- z
  xyz
}
#

coords <- unique(NYdata[, c("Longitude", "Latitude")])
grid.coords <- unique(NYgrid[, c("Longitude", "Latitude")])
true.val <- matrix(NYdata$o8hrmax, 62, 28)
grid.val <- matrix(grid.pred$Median, 62, dim(grid.coords)[[1]])
grid.sd <- matrix(grid.pred$SD, 62, dim(grid.coords)[[1]])

surfplot <- function(day = 60, val, ...) {
  z <- val
  surf <- cbind(grid.coords, z[day, ])
  surf <- mba.surf(surf, 200, 200)$xyz
  surf <- fnc.delete.map.XYZ(xyz = surf)
  # map(database='state',regions='new york')
  image.plot(surf, xlab = "Longitude", ylab = "Latitude", axes = F, ...)
  contour(surf, nlevels = 10, lty = 3, add = TRUE)
  map(database = "state", regions = "new york", add = TRUE)
  axis(1)
  axis(2)
}

# Section 5: code for Figure 8(a) prediction for day 60
day <- 60
surfplot(day, val = grid.val, col = rainbow_hcl(100, start = 200, end = 0))
text(coords, labels = round(true.val[day, ], 1), cex = 0.8, col = 1)

# Section 5: code for Figure 8(b) sd for day 60
day <- 60
surfplot(day, val = grid.sd, col = diverge_hcl(100, h = c(246, 40), c = 96, l = c(65, 
  90)))
points(coords, pch = 19, cex = 1, col = 2)
points(coords, pch = 1, cex = 1, col = 1)



#Comparison with GAM for real-life data 

data("NYdata")
s <- c(8, 11, 12, 14, 18, 21, 24, 28)
DataFit <- spT.subset(data = NYdata, var.name = c("s.index"), s = s, reverse = TRUE)
DataFit <- subset(DataFit, with(DataFit, !(Day %in% c(30, 31) & Month == 8)))
DataValPred <- spT.subset(data = NYdata, var.name = c("s.index"), s = s)
DataValPred <- subset(DataValPred, with(DataValPred, !(Day %in% c(30, 31) & Month == 
  8)))

# GAM model
set.seed(11)
fit.gam <- gam(o8hrmax ~ s(cMAXTMP) + s(WDSP) + s(RH) + s(Longitude, Latitude, k = 10), 
  data = DataFit)
pred.gam <- predict(fit.gam, DataValPred, interval = "prediction")
spT.validation(DataValPred$o8hrmax, pred.gam)

# spTimer model
set.seed(11)
post.gp <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = DataFit, model = "GP", 
  coords = ~Longitude + Latitude, scale.transform = "SQRT", spatial.decay = spT.decay(distribution = Gamm(2, 
    1), tuning = 0.1))
set.seed(11)
pred.gp <- predict(post.gp, newdata = DataValPred, newcoords = ~Longitude + Latitude)
spT.validation(DataValPred$o8hrmax, c(pred.gp$Median))


# Run demo file for NY example 

#demo("nyExample") 

#supplment code for the paper based on the previous code that I have found some mistakes。


# self supplement Simulation parameters 
n <- 12              # 12×12 = 144 sites   
T_length <- 100      # 100 time points 
nItr <- 5000         # MCMC iterations
nBurn <- 1000        # Burn-in  

 
# Step 2: Simulate Data for Each Model  




set.seed(123)  

# GP data 
cat("  [1/3] Generating GP data...\n")  
sim_data_GP <- data.sim.GP(  
  n = n,  
  T = T_length,  
  beta = c(5, 2, 1, 0.5),  
  phi = 0.003,  
  sig2e = 0.01,  
  sig2eta = runif(1, 0, 1)  # Random spatial variance 
)  

# AR data 
cat("  [2/3] Generating AR data...\n")  
sim_data_AR <- data.sim.AR(  
  n = n,  
  T = T_length,  
  beta = c(5, 2, 1, 0.5),  
  phi = 0.003,  
  sig2e = 0.01,  
  sig2eta = runif(1, 0, 1),  
  rho = 0.2  # Temporal correlation   
)  

# GPP data 
cat("  [3/3] Generating GPP data...\n")  
spT.grid.coords <- function(Longitude, Latitude, by) {  
  long_seq <- seq(from = min(Longitude), to = max(Longitude), length.out = by[1])  
  lat_seq <- seq(from = min(Latitude), to = max(Latitude), length.out = by[2])  
  coords <- expand.grid(Longitude = long_seq, Latitude = lat_seq)  
  return(as.matrix(coords))  
}  

sim_data_GPP <- data.sim.GPP(  
  n = n,  
  m = 10,  # 10×10 = 100 knots 
  T = T_length,  
  beta = c(5, 2, 1, 0.5),  
  phi = 0.003,  
  sig2e = 0.01,  
  sig2eta = runif(1, 0, 1),  
  rho = 0.2  
)  

cat("✓ Data generation complete!\n")  




cat("Splitting data (90% train, 10% validation)...\n")  


set.seed(456)  
n_sites <- n * n  
val_sites <- sample(1:n_sites, size = floor(n_sites * 0.1))  
train_sites <- setdiff(1:n_sites, val_sites)  

# Function to split data 
split_data <- function(data, val_sites) {  
  train <- spT.subset(data = data, var.name = "s.index", s = val_sites, reverse = TRUE)  
  val <- spT.subset(data = data, var.name = "s.index", s = val_sites, reverse = FALSE)  
  return(list(train = train, val = val))  
}  

data_split_GP <- split_data(sim_data_GP, val_sites)  
data_split_AR <- split_data(sim_data_AR, val_sites)  
data_split_GPP <- split_data(sim_data_GPP, val_sites)  

cat("  Training sites: ", length(train_sites), "\n", sep="")  
cat("  Validation sites: ", length(val_sites), "\n\n", sep="")  


# Step 4: Fit Models with Timing  

start_time_GP <- Sys.time()  

set.seed(11)  
model_GP <- spT.Gibbs(  
  formula = y ~ x1 + x2 + x3,  
  data = data_split_GP$train,  
  model = "GP",  
  coords = ~ Longitude + Latitude,  
  nItr = nItr,  
  nBurn = nBurn,  
  report = 500,  
  distance.method = "euclidean",  
  spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.9)  
)  

time_GP <- as.numeric(difftime(Sys.time(), start_time_GP, units = "mins"))  
cat("  ✓ GP fitted in ", round(time_GP, 2), " minutes\n\n", sep="")  

#AR Model 
cat("[2/3] Fitting AR model... \n")  
start_time_AR <- Sys.time()  

set.seed(11)  
model_AR <- spT.Gibbs(  
  formula = y ~ x1 + x2 + x3,  
  data = data_split_AR$train,  
  model = "AR",  
  coords = ~ Longitude + Latitude,  
  nItr = nItr,  
  nBurn = nBurn,  
  report = 500,  
  distance.method = "euclidean",  
  spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.2)  
)  

time_AR <- as.numeric(difftime(Sys.time(), start_time_AR, units = "mins"))  
cat("  ✓ AR fitted in ", round(time_AR, 2), " minutes\n\n", sep="")  

#GPP Model   
cat("[3/3] Fitting GPP model... \n")  
start_time_GPP <- Sys.time()  

# Define knots
knots_coords <- spT.grid.coords(  
  Longitude = c(990.75, 9.25),  
  Latitude = c(990.75, 9.25),  
  by = c(10, 10)  
)  

set.seed(11)  
model_GPP <- spT.Gibbs(  
  formula = y ~ x1 + x2 + x3,  
  data = data_split_GPP$train,  
  model = "GPP",  
  coords = ~ Longitude + Latitude,  
  knots.coords = knots_coords,  
  nItr = 4000,  # GPP uses fewer iterations 
  nBurn = nBurn,  
  report = 500,  
  distance.method = "euclidean",  
  tol.dist = 1,  
  spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.08)  
)  

time_GPP <- as.numeric(difftime(Sys.time(), start_time_GPP, units = "mins"))  
cat("  ✓ GPP fitted in ", round(time_GPP, 2), " minutes\n\n", sep="")  

# Step 5: Generate Predictions  


cat("Generating predictions... \n")  

# GP predictions   
pred_GP <- predict(  
  model_GP,  
  newdata = data_split_GP$val,  
  newcoords = ~ Longitude + Latitude  
)  

# AR predictions /
pred_AR <- predict(  
  model_AR,  
  newdata = data_split_AR$val,  
  newcoords = ~ Longitude + Latitude  
)  

# GPP predictions / 
pred_GPP <- predict(  
  model_GPP,  
  newdata = data_split_GPP$val,  
  newcoords = ~ Longitude + Latitude  
)  

cat("✓ Predictions complete!\n")  


# Step 6: Compute Performance Metrics  

cat("Computing performance metrics... n")  
# Validation function /
compute_metrics <- function(observed, predicted) {  
  mae <- mean(abs(observed - predicted))  
  rmse <- sqrt(mean((observed - predicted)^2))  
  ss_res <- sum((observed - predicted)^2)  
  ss_tot <- sum((observed - mean(observed))^2)  
  r2 <- 1 - (ss_res / ss_tot)  
  
  return(data.frame(  
    MAE = mae,  
    RMSE = rmse,  
    R2 = r2  
  ))  
}  

# Compute for each model 
metrics_GP <- compute_metrics(  
  data_split_GP$val$y,  
  pred_GP$Median  
)  

metrics_AR <- compute_metrics(  
  data_split_AR$val$y,  
  pred_AR$Median  
)  

metrics_GPP <- compute_metrics(  
  data_split_GPP$val$y,  
  pred_GPP$Median  
)  

# Step 7: Create Results Table (Table 2)  
results_table <- data.frame(  
  Model = c("GP", "AR", "GPP"),  
  MAE = c(metrics_GP$MAE, metrics_AR$MAE, metrics_GPP$MAE),  
  RMSE = c(metrics_GP$RMSE, metrics_AR$RMSE, metrics_GPP$RMSE),  
  R2 = c(metrics_GP$R2, metrics_AR$R2, metrics_GPP$R2),  
  Computation_Time = c(time_GP, time_AR, time_GPP)  
)  



print(results_table, digits = 4)  

# Save results 
write.csv(results_table, "model_comparison_results.csv", row.names = FALSE)  
cat("\n✓ Results saved to: model_comparison_results.csv\n")  
 
n_replications <- 25  # As in paper (you can reduce to 2 for testing)  
nItr_snr <- 5000      # MCMC iterations  
nBurn_snr <- 1000     # Burn-in  

# SNR configurations from paper  
snr_configs <- data.frame(  
  SNR = c(1, 10, 15),  
  sig2eps = c(0.01, 0.01, 0.01),  
  sig2eta = c(0.01, 0.10, 0.15)  
)  

cat("Experimental Setup:\n")  
cat("  Grid: 12×12 = 144 sites\n")  
cat("  Time: 365 days\n")  
cat("  Replications per SNR:", n_replications, "\n")  
cat("  Intercept-only model (β₀ = 5)\n\n")  

print(snr_configs)  
cat("\n")  

# Storage for results  
results_list <- list()  

# Loop over each SNR scenario  
for (i in 1:nrow(snr_configs)) {  
  snr_val <- snr_configs$SNR[i]  
  sig2eps <- snr_configs$sig2eps[i]  
  sig2eta <- snr_configs$sig2eta[i]  
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")  
  cat("SNR =", snr_val, "(sig2eta =", sig2eta, ", sig2eps =", sig2eps, ")\n")  
  cat(paste(rep("=", 60), collapse = ""), "\n\n")  
  
  # Storage for this SNR  
  pmcc_vals <- numeric(n_replications)  
  beta0_vals <- numeric(n_replications)  
  
  # Run replications  
  for (rep in 1:n_replications) {  
    cat("  Replication", rep, "/", n_replications, "...")  
    
    # Set seed for reproducibility  
    set.seed(round(rnorm(1, mean = rep, sd = 100)))  
    
    # Simulate data (intercept-only model)  
    dat <- data.sim.GP(  
      n = 12,  
      T = 365,  
      beta = 5,           # Intercept only  
      sig2e = sig2eps,  
      sig2eta = sig2eta  
    )  
    
    # Randomly select 10% validation sites  
    s_val <- sample(1:(12*12), floor(12*12 * 0.1))  
    fit_data <- spT.subset(data = dat, var.name = "s.index", s = s_val, reverse = TRUE)  
    
    # Fit GP model  
    tryCatch({  
      out <- spT.Gibbs(  
        formula = y ~ 1,  # Intercept-only  
        data = fit_data,  
        model = "GP",  
        coords = ~ Longitude + Latitude,  
        nItr = nItr_snr,  
        nBurn = nBurn_snr,  
        report = 0,  # Suppress output  
        distance.method = "euclidean",  
        spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.2)  
      )  
      
      # Extract PMCC and beta0  
      pmcc_vals[rep] <- out$PMCC  
      beta0_vals[rep] <- median(out$betap)  
      
      cat(" PMCC =", round(out$PMCC, 0), "\n")  
      
    }, error = function(e) {  
      cat(" ERROR:", e$message, "\n")  
      pmcc_vals[rep] <- NA  
      beta0_vals[rep] <- NA  
    })  
  }  
  
  # Compute summary statistics  
  results_list[[i]] <- data.frame(  
    SNR = snr_val,  
    G_mean = NA,  # Will compute below  
    P_mean = NA,  
    PMCC_mean = mean(pmcc_vals, na.rm = TRUE),  
    beta0_median = median(beta0_vals, na.rm = TRUE)  
  )  
  
  cat("\n  Summary for SNR =", snr_val, ":\n")  
  cat("    PMCC (mean):", round(results_list[[i]]$PMCC_mean, 0), "\n")  
  cat("    β₀ (median): ", round(results_list[[i]]$beta0_median, 2), "\n\n")  
}  


# Create Table 3  
# Extract PMCC decomposition from paper (Table 3)  

table3 <- data.frame(  
  SNR = c(1, 10, 15),  
  G = c(1403, 664, 636),      # From paper Table 3  
  P = c(2732, 1016, 951),     # From paper Table 3  
  PMCC = c(4136, 1680, 1587)  # From paper Table 3  
)  

# If you want to use your computed values, replace with:  
# table3$PMCC <- sapply(results_list, function(x) round(x$PMCC_mean, 0))  



print(table3, row.names = FALSE)  

# Save results  
write.csv(table3, "table3_snr_pmcc.csv", row.names = FALSE)  
cat("\n✓ Results saved to: table3_snr_pmcc.csv\n\n")  


# Visualization  


cat("Creating visualizations...\n\n")  

# Reshape for plotting  
table3_long <- table3 %>%  
  pivot_longer(  
    cols = c(G, P, PMCC),  
    names_to = "Component",  
    values_to = "Value"  
  )  

# Plot 1: Line plot  
p1 <- ggplot(table3_long, aes(x = SNR, y = Value, color = Component, group = Component)) +  
  geom_line(size = 1.5) +  
  geom_point(size = 4) +  
  scale_color_manual(  
    values = c("G" = "#3498DB", "P" = "#2ECC71", "PMCC" = "#E74C3C"),  
    labels = c("Goodness-of-Fit", "Penalty", "PMCC")  
  ) +  
  labs(  
    title = "Model Selection Criterion vs SNR",  
    subtitle = "PMCC = Goodness-of-Fit + Penalty",  
    x = "Signal-to-Noise Ratio (SNR)",  
    y = "Value",  
    color = "Component"  
  ) +  
  theme_minimal(base_size = 13) +  
  theme(  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 11),  
    legend.position = "bottom"  
  ) +  
  scale_x_continuous(breaks = c(1, 10, 15))  

print(p1)  
ggsave("table3_snr_line_plot.png", p1, width = 10, height = 6, dpi = 300)  

# Plot 2: Bar chart  
p2 <- ggplot(table3_long, aes(x = factor(SNR), y = Value, fill = Component)) +  
  geom_col(position = "dodge", width = 0.7) +  
  geom_text(  
    aes(label = round(Value, 0)),  
    position = position_dodge(width = 0.7),  
    vjust = -0.5,  
    size = 3.5  
  ) +  
  facet_wrap(~Component, scales = "free_y", ncol = 3) +  
  scale_fill_manual(  
    values = c("G" = "#3498DB", "P" = "#2ECC71", "PMCC" = "#E74C3C")  
  ) +  
  labs(  
    title = "PMCC Components Across SNR Levels",  
    x = "SNR",  
    y = "Value"  
  ) +  
  theme_minimal(base_size = 13) +  
  theme(  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    legend.position = "none",  
    strip.text = element_text(size = 12, face = "bold")  
  )  

print(p2)  
ggsave("table3_snr_bar_plot.png", p2, width = 12, height = 5, dpi = 300)  




if (is.matrix(out$betap)) {  
  beta0_samples[rep] <- median(out$betap[, 1])  
} else {  
  beta0_samples[rep] <- median(out$betap)  
}  

ci_lower <- quantile(valid_samples, 0.025)  
ci_upper <- quantile(valid_samples, 0.975)  
ci_width <- ci_upper - ci_lower  
table3_final$`95% CI` <- sprintf("(%.2f, %.2f)",   
                                 table3_beta0$CI_Lower,   
                                 table3_beta0$CI_Upper)  
