
# Robustness analysis: simulate data and fit GP models under different SNR (signal-to-noise ratio) values
setwd('D:/python project')
source("D:/python project/data_simulation.R")
source("D:/python project/model_fitting.R")
library('spTimer')
# Set SNR test values
snr_values <- c(1, 10, 15)
results <- list()

for (snr in snr_values) {
  cat("SNR =", snr, "\n")
  
  # Set signal and noise variances
  sig2eta <- 1
  sig2eps <- sig2eta / snr
  
  # Simulate data
  sim_data <- data.sim.GP(sig2eps = sig2eps, sig2eta = sig2eta)
  
  # Fit the GP model
  model_result <- fit_model(model_type = "GP", data.sim = sim_data, nItr = 200, nBurn = 100)
  
  # Save results
  results[[paste0("SNR_", snr)]] <- model_result
}

# Print all model results
print(results)
# Create summary of results
snr_values <- c(1, 10, 15)
gof <- c(1403.39, 664.16, 636.12)
penalty <- c(2732.49, 1016.29, 950.62)
pmcc <- c(4135.88, 1680.45, 1586.74)

# Set up plot layout
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Plot 1: PMCC vs SNR
plot(snr_values, pmcc, type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "SNR", ylab = "PMCC", main = "Model Selection Criterion vs SNR",
     cex = 1.5, cex.lab = 1.2)
grid()
text(snr_values, pmcc, labels = round(pmcc, 1), pos = 3, cex = 0.9)

# Plot 2: Goodness of Fit vs SNR
plot(snr_values, gof, type = "b", pch = 19, col = "blue", lwd = 2,
     xlab = "SNR", ylab = "Goodness of Fit", main = "Goodness of Fit vs SNR",
     cex = 1.5, cex.lab = 1.2)
grid()
text(snr_values, gof, labels = round(gof, 1), pos = 3, cex = 0.9)

# Plot 3: Penalty vs SNR
plot(snr_values, penalty, type = "b", pch = 19, col = "darkgreen", lwd = 2,
     xlab = "SNR", ylab = "Penalty", main = "Penalty Term vs SNR",
     cex = 1.5, cex.lab = 1.2)
grid()
text(snr_values, penalty, labels = round(penalty, 1), pos = 3, cex = 0.9)

# Plot 4: Stacked bar chart of components
barplot(rbind(gof, penalty), beside = TRUE, col = c("blue", "darkgreen"),
        names.arg = snr_values, xlab = "SNR", ylab = "Value",
        main = "PMCC Components", legend = c("Goodness of Fit", "Penalty"),
        args.legend = list(x = "topright"))

