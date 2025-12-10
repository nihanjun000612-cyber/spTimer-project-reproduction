
library(dplyr)
library(ggplot2)
library(spTimer)

# Function to compute validation metrics
compute_metrics <- function(observed, predicted) {
  mae <- mean(abs(observed - predicted))
  rmse <- sqrt(mean((observed - predicted)^2))
  r2 <- 1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
  return(data.frame(MAE = mae, RMSE = rmse, R2 = r2))
}

# Plot model comparison
plot_model_metrics <- function(results_df) {
  results_long <- results_df %>%
    tidyr::pivot_longer(cols = c(MAE, RMSE, R2), names_to = "Metric", values_to = "Value")
  
  ggplot(results_long, aes(x = Model, y = Value, fill = Model)) +
    geom_bar(stat = "identity", width = 0.6) +
    facet_wrap(~Metric, scales = "free_y") +
    geom_text(aes(label = round(Value, 4)), vjust = -0.3, size = 3.5) +
    labs(title = "Model Performance Comparison (GP vs AR vs GPP)", 
         y = "Metric Value",
         subtitle = "Lower MAE/RMSE and Higher R² indicate better performance") +
    theme_minimal() +
    theme(legend.position = "bottom",
          text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
}


# Set working directory and load functions
setwd('D:/python project')
source("D:/python project/data_simulation.R")
source("D:/python project/model_fitting.R")

# Simulate the same data once
set.seed(123)
sim_data <- data.sim.GP(sig2eps = 1, sig2eta = 1)  # SNR = 1


# Visualize site locations
plot(sim_data$coords, pch = 19, col = "blue", cex = 1.5,
     main = "Simulated Site Locations",
     xlab = "Longitude", ylab = "Latitude")
text(sim_data$coords, labels = 1:nrow(sim_data$coords), pos = 3, cex = 0.8)
grid()

# Fit GP Model
cat("\n========================================\n")
cat("Fitting GP Model\n")
cat("========================================\n")
model_GP <- fit_model("GP", data.sim = sim_data, nItr = 200, nBurn = 100)
cat("GP model fitting completed.\n")

# Fit AR Model
cat("\n========================================\n")
cat("Fitting AR Model\n")
cat("========================================\n")
model_AR <- fit_model("AR", data.sim = sim_data, nItr = 200, nBurn = 100)
cat("AR model fitting completed.\n")

# Fit GPP Model
model_GPP <- fit_model("GPP", data.sim = sim_data, nItr = 200, nBurn = 100, n_knots = 4)
cat("GPP model fitting completed.\n")

# Extract observed values
obs_Y <- as.vector(t(sim_data$Y))

# Debug: Check model structure
cat("\n=== Debugging Model Output Structure ===\n")
cat("GP model structure:\n")
cat("  Names:", names(model_GP), "\n")
cat("  Class:", class(model_GP), "\n")

# Extract fitted values - try multiple approaches
cat("\n=== Extracting Predictions ===\n")

# Method 1: Try fitted.values
if (!is.null(model_GP$fitted.values) && length(model_GP$fitted.values) > 0) {
  pred_GP <- model_GP$fitted.values
  pred_AR <- model_AR$fitted.values
  pred_GPP <- model_GPP$fitted.values
  cat("Using fitted.values\n")
  
  # Method 2: Try posterior means
} else if (!is.null(model_GP$fitted) && length(model_GP$fitted) > 0) {
  pred_GP <- model_GP$fitted
  pred_AR <- model_AR$fitted
  pred_GPP <- model_GPP$fitted
  cat(" Using fitted\n")
  
  # Method 3: Compute from posterior samples
} else {
  cat("Computing predictions from posterior samples...\n")
  
  # For GP model
  if (!is.null(model_GP$mcmc.beta)) {
    beta_GP <- colMeans(model_GP$mcmc.beta)
    X_matrix <- sim_data$X
    pred_GP <- as.vector(X_matrix %*% beta_GP)
    pred_GP <- rep(pred_GP, each = ncol(sim_data$Y))
  } else {
    stop("Cannot extract predictions from GP model")
  }
  
  # For AR model
  if (!is.null(model_AR$mcmc.beta)) {
    beta_AR <- colMeans(model_AR$mcmc.beta)
    pred_AR <- as.vector(X_matrix %*% beta_AR)
    pred_AR <- rep(pred_AR, each = ncol(sim_data$Y))
  } else {
    stop("Cannot extract predictions from AR model")
  }
  
  # For GPP model
  if (!is.null(model_GPP$mcmc.beta)) {
    beta_GPP <- colMeans(model_GPP$mcmc.beta)
    pred_GPP <- as.vector(X_matrix %*% beta_GPP)
    pred_GPP <- rep(pred_GPP, each = ncol(sim_data$Y))
  } else {
    stop("Cannot extract predictions from GPP model")
  }
  
  cat(" Using posterior mean coefficients\n")
}

# Verify lengths
cat("\nData lengths:\n")
cat("  Observed:", length(obs_Y), "\n")
cat("  GP predictions:", length(pred_GP), "\n")
cat("  AR predictions:", length(pred_AR), "\n")
cat("  GPP predictions:", length(pred_GPP), "\n")

# Make sure lengths match
if (length(pred_GP) != length(obs_Y)) {
  cat("Warning: Length mismatch, attempting to fix...\n")
  
  # Try to get the right structure
  n_sites <- nrow(sim_data$coords)
  n_times <- ncol(sim_data$Y)
  
  # Repeat predictions for each time point if needed
  if (length(pred_GP) == n_sites) {
    pred_GP <- rep(pred_GP, each = n_times)
    pred_AR <- rep(pred_AR, each = n_times)
    pred_GPP <- rep(pred_GPP, each = n_times)
  }
  
  # Truncate if too long
  if (length(pred_GP) > length(obs_Y)) {
    pred_GP <- pred_GP[1:length(obs_Y)]
    pred_AR <- pred_AR[1:length(obs_Y)]
    pred_GPP <- pred_GPP[1:length(obs_Y)]
  }
}


gp_metrics <- compute_metrics(obs_Y, pred_GP)
gp_metrics$Model <- "GP"

ar_metrics <- compute_metrics(obs_Y, pred_AR)
ar_metrics$Model <- "AR"

gpp_metrics <- compute_metrics(obs_Y, pred_GPP)
gpp_metrics$Model <- "GPP"

results_df <- rbind(gp_metrics, ar_metrics, gpp_metrics)

# Reorder columns
results_df <- results_df[, c("Model", "MAE", "RMSE", "R2")]



# Identify best model for each metric
best_mae <- results_df[which.min(results_df$MAE), "Model"]
best_rmse <- results_df[which.min(results_df$RMSE), "Model"]
best_r2 <- results_df[which.max(results_df$R2), "Model"]


# Save results
write.csv(results_df, "model_comparison_results.csv", row.names = FALSE)
cat(" Results saved to: model_comparison_results.csv\n\n")

p <- plot_model_metrics(results_df)
print(p)

ggsave("model_comparison_plot.png", plot = p, width = 12, height = 5, dpi = 300)
cat(" Plot saved to: model_comparison_plot.png\n\n")

# Additional visualization: Observed vs Predicted (only if predictions available)
if (length(pred_GP) == length(obs_Y)) {
  comparison_df <- data.frame(
    Observed = rep(obs_Y, 3),
    Predicted = c(pred_GP, pred_AR, pred_GPP),
    Model = factor(rep(c("GP", "AR", "GPP"), each = length(obs_Y)),
                   levels = c("GP", "AR", "GPP"))
  )
  
  p_scatter <- ggplot(comparison_df, aes(x = Observed, y = Predicted, color = Model)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.8) +
    facet_wrap(~Model, ncol = 3) +
    labs(title = "Observed vs Predicted Values",
         subtitle = "Points closer to diagonal line indicate better predictions",
         x = "Observed Values", y = "Predicted Values") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  print(p_scatter)
  ggsave("observed_vs_predicted.png", plot = p_scatter, width = 12, height = 4, dpi = 300)
  cat(" Scatter plot saved to: observed_vs_predicted.png\n\n")
} else {
  cat("Skipping scatter plot due to dimension mismatch\n\n")
}
cat("  1. model_comparison_results.csv\n")
cat("  2. model_comparison_plot.png\n")
if (length(pred_GP) == length(obs_Y)) {
  cat("  3. observed_vs_predicted.png\n")
}
cat("\n")