

library(MASS)
library(spTimer)
library(dplyr)
library(ggplot2)

setwd('D:/python project')
source("D:/python project/data_simulation.R")
source("D:/python project/nonseparable_cov.R")
source("D:/python project/model_fitting.R")

compute_metrics <- function(observed, predicted) {
  min_len <- min(length(observed), length(predicted))
  observed <- observed[1:min_len]
  predicted <- predicted[1:min_len]
  
  mae <- mean(abs(observed - predicted))
  rmse <- sqrt(mean((observed - predicted)^2))
  r2 <- 1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
  
  return(data.frame(MAE = mae, RMSE = rmse, R2 = r2))
}

plot_metrics_balanced <- function(results_df) {
  results_long <- results_df %>%
    tidyr::pivot_longer(cols = c(MAE, RMSE, R2), 
                        names_to = "Metric", 
                        values_to = "Value")
  
  results_long$Model_Label <- ifelse(
    grepl("Nonseparable", results_long$Model),
    "Nonseparable\n(Misspecified)",
    "Separable\n(Correct)"
  )
  
  ggplot(results_long, aes(x = Model_Label, y = Value, fill = Model)) +
    geom_col(width = 0.7, color = NA) +
    geom_text(aes(label = sprintf("%.3f", Value)), 
              vjust = -0.5, size = 5, fontface = "bold", color = "black") +
    facet_wrap(~Metric, scales = "free_y", ncol = 3) +
    scale_fill_manual(
      values = c("Nonseparable (Misspecified)" = "#FF7675", 
                 "Separable (Correct)" = "#00CEC9"),
      name = NULL
    ) +
    labs(title = "Model Performance Comparison",
         subtitle = "Impact of Covariance Misspecification",
         x = NULL, y = "Metric Value") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 13),
      strip.text = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
}

# IMPROVED: Extract predictions with full diagnostics
extract_predictions_diagnostic <- function(model, sim_data, model_name) {
  n_sites <- nrow(sim_data$coords)
  n_times <- ncol(sim_data$Y)
  
  cat("\n  ===", model_name, "===\n")
  cat("  Checking model structure...\n")
  
  # Check what's in the model
  cat("  Model components:", paste(names(model), collapse = ", "), "\n")
  
  # Check fitted values
  if (!is.null(model$fitted)) {
    fitted_data <- model$fitted
    cat("  model$fitted class:", class(fitted_data), "\n")
    
    if (is.matrix(fitted_data)) {
      cat("  Matrix dimensions:", nrow(fitted_data), "×", ncol(fitted_data), "\n")
      
      if (ncol(fitted_data) == 2) {
        col1 <- fitted_data[, 1]
        col2 <- fitted_data[, 2]
        
        cat("  Column 1: range = [", paste(round(range(col1), 3), collapse = ", "), 
            "], sd =", round(sd(col1), 4), "\n")
        cat("  Column 2: range = [", paste(round(range(col2), 3), collapse = ", "), 
            "], sd =", round(sd(col2), 4), "\n")
        
        # Choose column with more variation
        if (sd(col2) > sd(col1) * 1.1) {
          pred <- col2
          cat("  → Using Column 2 (more variation)\n")
        } else {
          pred <- col1
          cat("  → Using Column 1\n")
        }
        
        # Check if predictions are reasonable
        obs_range <- range(as.vector(sim_data$Y))
        pred_range <- range(pred)
        
        cat("  Observed data range: [", paste(round(obs_range, 2), collapse = ", "), "]\n")
        cat("  Predicted range:     [", paste(round(pred_range, 2), collapse = ", "), "]\n")
        
        # Warning if predictions are too flat
        if (sd(pred) < sd(as.vector(sim_data$Y)) * 0.1) {
          cat("  ⚠ WARNING: Predictions have very low variance!\n")
          cat("  This suggests the model is only capturing fixed effects.\n")
        }
        
        return(pred)
      }
    }
  }
  
  stop("Cannot extract predictions!")
}

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║     Nonseparable Covariance Test (FULL VERSION)      ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

# Step 1: Simulate data
cat("Step 1: Simulating separable GP data...\n")
set.seed(123)
sim_data_sep <- data.sim.GP(sig2eps = 1, sig2eta = 1)

coords <- sim_data_sep$coords
n_sites <- nrow(coords)
n_times <- ncol(sim_data_sep$Y)

cat("  Sites:", n_sites, ", Times:", n_times, ", Total:", n_sites * n_times, "\n")
cat("  Observed Y: mean =", round(mean(sim_data_sep$Y), 2), 
    ", sd =", round(sd(as.vector(sim_data_sep$Y)), 2), "\n\n")

# Step 2: Create nonseparable data
cat("Step 2: Creating nonseparable covariance data...\n")

time_points <- 1:n_times
cov_mat_nonsep <- nonsep_cov_fast(coords, time_points, phi_s = 0.8, phi_t = 0.08, sig2 = 2)

set.seed(456)
eta_nonsep <- mvrnorm(1, mu = rep(0, nrow(cov_mat_nonsep)), Sigma = cov_mat_nonsep)
eta_matrix <- matrix(eta_nonsep, nrow = n_sites, ncol = n_times, byrow = FALSE)

X_matrix <- sim_data_sep$X
beta_true <- c(5, 2, 1, 0.5)
mu_vector <- as.vector(X_matrix %*% beta_true)
mu_matrix <- matrix(rep(mu_vector, n_times), nrow = n_sites, ncol = n_times)

set.seed(789)
epsilon_matrix <- matrix(rnorm(n_sites * n_times, 0, sqrt(1)), nrow = n_sites, ncol = n_times)

Y_nonsep <- mu_matrix + eta_matrix + epsilon_matrix

sim_data_nonsep <- sim_data_sep
sim_data_nonsep$Y <- Y_nonsep

cat("  Nonseparable Y: mean =", round(mean(Y_nonsep), 2), 
    ", sd =", round(sd(as.vector(Y_nonsep)), 2), "\n\n")

# Step 3: Fit models
cat("Step 3: Fitting GP models...\n")

cat("  [1/2] Fitting on separable data...\n")
model_sep <- fit_model("GP", data.sim = sim_data_sep, nItr = 500, nBurn = 200)

cat("\n  [2/2] Fitting on nonseparable data...\n")
model_nonsep <- fit_model("GP", data.sim = sim_data_nonsep, nItr = 500, nBurn = 200)

cat("\n✓ Both models fitted\n")

# Step 4: Extract predictions WITH DIAGNOSTICS
cat("\n")
cat(paste(rep("=", 60), collapse = ""))
cat("\nStep 4: DIAGNOSTIC - Extracting Predictions\n")
cat(paste(rep("=", 60), collapse = ""))
cat("\n")

obs_sep <- as.vector(t(sim_data_sep$Y))
obs_nonsep <- as.vector(t(sim_data_nonsep$Y))

pred_sep <- extract_predictions_diagnostic(model_sep, sim_data_sep, "SEPARABLE MODEL")
pred_nonsep <- extract_predictions_diagnostic(model_nonsep, sim_data_nonsep, "NONSEPARABLE MODEL")

cat("\n")
cat(paste(rep("=", 60), collapse = ""))
cat("\n")

# Step 5: Compute metrics
cat("\nStep 5: Computing performance metrics...\n")

metrics_nonsep <- compute_metrics(obs_nonsep, pred_nonsep)
metrics_nonsep$Model <- "Nonseparable (Misspecified)"

metrics_sep <- compute_metrics(obs_sep, pred_sep)
metrics_sep$Model <- "Separable (Correct)"

results_df <- rbind(metrics_nonsep, metrics_sep)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║                    Results                            ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

print(results_df, row.names = FALSE, digits = 4)

mae_change <- (metrics_nonsep$MAE - metrics_sep$MAE) / metrics_sep$MAE * 100
cat("\nPerformance Impact: MAE changed by", sprintf("%+.2f%%", mae_change), "\n\n")

# Step 6: Create ALL visualizations
cat("Step 6: Creating visualizations...\n\n")

# Plot 1: Metrics
cat("  [1/5] Metric comparison...\n")
p1 <- plot_metrics_balanced(results_df)
print(p1)
ggsave("nonsep_metrics_comparison.png", p1, width = 14, height = 5.5, dpi = 300, bg = "white")
cat("      ✓ Saved\n\n")

# Plot 2: First scatter (the "flat line" problem)
cat("  [2/5] Initial scatter plot (diagnostic view)...\n")

diagnostic_df <- data.frame(
  Observed = c(obs_nonsep, obs_sep),
  Predicted = c(pred_nonsep, pred_sep),
  Model = rep(c("Nonseparable Data\n(Misspecified)", "Separable Data\n(Correct)"),
              c(length(obs_nonsep), length(obs_sep)))
)

p2_diagnostic <- ggplot(diagnostic_df, aes(x = Observed, y = Predicted, color = Model)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 1) +
  facet_wrap(~Model, ncol = 2) +
  scale_color_manual(values = c("#FF7675", "#00CEC9")) +
  labs(
    title = "Observed vs Predicted: Impact of Covariance Misspecification",
    subtitle = "Note: Flat predictions indicate model is only capturing fixed effects",
    x = "Observed Values", y = "Predicted Values"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "red"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

print(p2_diagnostic)
ggsave("nonsep_obs_vs_pred_diagnostic.png", p2_diagnostic, 
       width = 12, height = 6, dpi = 300, bg = "white")
cat("      ✓ Saved (diagnostic version)\n\n")

# Plot 3: Better scatter (zoomed in)
cat("  [3/5] Detailed scatter plot (zoomed)...\n")

cor_nonsep <- cor(obs_nonsep, pred_nonsep)
cor_sep <- cor(obs_sep, pred_sep)

p3 <- ggplot(diagnostic_df, aes(x = Observed, y = Predicted, color = Model)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 1.2) +
  facet_wrap(~Model, ncol = 2, scales = "free") +
  scale_color_manual(values = c("#FF7675", "#00CEC9")) +
  labs(
    title = "Observed vs Predicted (Zoomed View)",
    subtitle = sprintf("Correlation: Nonseparable = %.3f | Separable = %.3f", cor_nonsep, cor_sep),
    x = "Observed Values", y = "Predicted Values"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

print(p3)
ggsave("nonsep_obs_vs_pred_zoomed.png", p3, width = 12, height = 6, dpi = 300, bg = "white")
cat("      ✓ Saved (zoomed version)\n\n")

# Plot 4: Residuals
cat("  [4/5] Residual analysis...\n")

residual_df <- data.frame(
  Index = c(1:length(obs_nonsep), 1:length(obs_sep)),
  Residual = c(obs_nonsep - pred_nonsep, obs_sep - pred_sep),
  Model = rep(c("Nonseparable (Misspecified)", "Separable (Correct)"), 
              c(length(obs_nonsep), length(obs_sep)))
)

p4 <- ggplot(residual_df, aes(x = Index, y = Residual, color = Model)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_smooth(se = TRUE, alpha = 0.2, size = 1.2, method = "loess") +
  facet_wrap(~Model, ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("#FF7675", "#00CEC9")) +
  labs(
    title = "Residual Analysis",
    subtitle = "Systematic patterns indicate model misspecification",
    x = "Observation Index", y = "Residual"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

print(p4)
ggsave("nonsep_residuals.png", p4, width = 12, height = 7.5, dpi = 300, bg = "white")
cat("      ✓ Saved\n\n")

# Plot 5: Covariance comparison
cat("  [5/5] Covariance structures...\n")

png("nonsep_covariance_comparison.png", width = 2200, height = 950, res = 180)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 2), cex.main = 1.6, cex.lab = 1.3, cex.axis = 1.1)

n_show <- min(150, nrow(cov_mat_nonsep))

# Nonseparable
cor_nonsep_mat <- cov2cor(cov_mat_nonsep[1:n_show, 1:n_show])
image(1:n_show, 1:n_show, cor_nonsep_mat,
      main = "Nonseparable Correlation Structure",
      xlab = "Space-Time Index", ylab = "Space-Time Index",
      col = heat.colors(100))
box(lwd = 2.5)

# Separable
dist_s <- as.matrix(dist(coords))
dist_t <- as.matrix(dist(time_points))
cov_s <- exp(-1 * dist_s)
cov_t <- exp(-0.1 * dist_t)
cov_sep_mat <- kronecker(cov_t, cov_s)
cor_sep_mat <- cov2cor(cov_sep_mat[1:n_show, 1:n_show])

image(1:n_show, 1:n_show, cor_sep_mat,
      main = "Separable Correlation Structure",
      xlab = "Space-Time Index", ylab = "Space-Time Index",
      col = heat.colors(100))
box(lwd = 2.5)

dev.off()
cat("      ✓ Saved\n\n")

write.csv(results_df, "nonsep_test_results.csv", row.names = FALSE)

cat("Generated 5 plots:\n")
cat("  1. nonsep_metrics_comparison.png         - Performance metrics\n")
cat("  2. nonsep_obs_vs_pred_diagnostic.png     - Full range view (shows flat predictions)\n")
cat("  3. nonsep_obs_vs_pred_zoomed.png         - Zoomed view (better detail)\n")
cat("  4. nonsep_residuals.png                  - Residual patterns\n")
cat("  5. nonsep_covariance_comparison.png      - Covariance structures\n\n")

cat("IMPORTANT NOTE:\n")
cat("  If predictions are flat (no variation), this means:\n")
cat("  → The GP model is only capturing fixed effects (β)\n")
cat("  → Spatial-temporal random effects (η) are not being predicted\n")
cat("  → This could be due to model convergence issues or data structure\n\n")

cat("Performance Impact: ", sprintf("%+.2f%%", mae_change), " MAE change\n\n")