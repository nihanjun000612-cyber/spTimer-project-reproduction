
  
# ==== Grid vs Random Sampling with Advanced Diagnostics (FULLY FIXED) ====

library(MASS)
library(spTimer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

setwd('D:/python project')
source("D:/python project/data_simulation.R")
source("D:/python project/model_fitting.R")

# ============================================================
# Helper Functions
# ============================================================

compute_metrics <- function(observed, predicted) {
  min_len <- min(length(observed), length(predicted))
  observed <- observed[1:min_len]
  predicted <- predicted[1:min_len]
  
  mae <- mean(abs(observed - predicted))
  rmse <- sqrt(mean((observed - predicted)^2))
  r2 <- 1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
  
  return(data.frame(MAE = mae, RMSE = rmse, R2 = r2))
}

generate_grid_coords <- function(n_sites) {
  n_side <- ceiling(sqrt(n_sites))
  grid <- expand.grid(
    lon = seq(0, 1, length.out = n_side),
    lat = seq(0, 1, length.out = n_side)
  )
  return(grid[1:n_sites, ])
}

generate_random_coords <- function(n_sites, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data.frame(
    lon = runif(n_sites, 0, 1),
    lat = runif(n_sites, 0, 1)
  )
}

extract_predictions_safe <- function(model, sim_data, model_name) {
  cat("  Extracting predictions for", model_name, "...\n")
  
  if (!is.null(model$fitted) && is.matrix(model$fitted)) {
    if (ncol(model$fitted) == 2) {
      col1 <- model$fitted[, 1]
      col2 <- model$fitted[, 2]
      
      if (sd(col2) > sd(col1) * 1.1) {
        pred <- col2
        cat("    Using Column 2 (sd =", round(sd(col2), 4), ")\n")
      } else {
        pred <- col1
        cat("    Using Column 1 (sd =", round(sd(col1), 4), ")\n")
      }
      
      return(pred)
    }
  }
  
  stop("Cannot extract predictions from model!")
}

# ============================================================
# Visualization Functions
# ============================================================

plot_metrics_comparison <- function(metrics_df) {
  metrics_long <- metrics_df %>%
    pivot_longer(cols = c(MAE, RMSE, R2), 
                 names_to = "Metric", 
                 values_to = "Value")
  
  ggplot(metrics_long, aes(x = Model, y = Value, fill = Model)) +
    geom_col(width = 0.6, color = "black", size = 0.3) +
    geom_text(aes(label = sprintf("%.3f", Value)), 
              vjust = -0.5, size = 5, fontface = "bold") +
    facet_wrap(~Metric, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("Grid Sites" = "#6C5CE7", 
                                 "Random Sites" = "#00B894")) +
    labs(title = "Grid vs Random Sampling: Model Performance",
         subtitle = "Impact of spatial sampling design on GP predictions",
         x = NULL, y = "Metric Value") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      strip.text = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
}

plot_site_locations <- function(grid_coords, random_coords) {
  coords_df <- rbind(
    data.frame(grid_coords, Type = "Grid"),
    data.frame(random_coords, Type = "Random")
  )
  
  ggplot(coords_df, aes(x = lon, y = lat, color = Type, shape = Type)) +
    geom_point(size = 4, alpha = 0.8) +
    facet_wrap(~Type, ncol = 2) +
    scale_color_manual(values = c("Grid" = "#6C5CE7", "Random" = "#00B894")) +
    scale_shape_manual(values = c("Grid" = 15, "Random" = 16)) +
    labs(title = "Spatial Sampling Designs",
         subtitle = sprintf("Grid (regular) vs Random (uniform) - %d sites each", 
                            nrow(grid_coords)),
         x = "Longitude", y = "Latitude") +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}

plot_obs_vs_pred <- function(obs_grid, pred_grid, obs_random, pred_random) {
  comparison_df <- data.frame(
    Observed = c(obs_grid, obs_random),
    Predicted = c(pred_grid, pred_random),
    Model = rep(c("Grid Sites", "Random Sites"),
                c(length(obs_grid), length(obs_random)))
  )
  
  cor_grid <- cor(obs_grid, pred_grid)
  cor_random <- cor(obs_random, pred_random)
  
  ggplot(comparison_df, aes(x = Observed, y = Predicted, color = Model)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "black", size = 1) +
    facet_wrap(~Model, ncol = 2) +
    scale_color_manual(values = c("Grid Sites" = "#6C5CE7", 
                                  "Random Sites" = "#00B894")) +
    labs(
      title = "Observed vs Predicted Values",
      subtitle = sprintf("Correlation: Grid = %.3f | Random = %.3f", 
                         cor_grid, cor_random),
      x = "Observed Values", y = "Predicted Values"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
}

plot_distance_comparison <- function(grid_coords, random_coords) {
  dist_grid <- as.matrix(dist(grid_coords))
  dist_random <- as.matrix(dist(random_coords))
  
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 2))
  
  image(1:nrow(dist_grid), 1:ncol(dist_grid), dist_grid,
        main = "Grid: Pairwise Distances",
        xlab = "Site Index", ylab = "Site Index",
        col = hcl.colors(100, "YlOrRd", rev = TRUE),
        cex.main = 1.4, cex.lab = 1.2)
  box(lwd = 2)
  
  image(1:nrow(dist_random), 1:ncol(dist_random), dist_random,
        main = "Random: Pairwise Distances",
        xlab = "Site Index", ylab = "Site Index",
        col = hcl.colors(100, "YlOrRd", rev = TRUE),
        cex.main = 1.4, cex.lab = 1.2)
  box(lwd = 2)
}

plot_residual_surface <- function(observed, predicted, coords, n_sites, n_times, 
                                  title = "Residual Surface (Time-Averaged)") {
  obs_matrix <- matrix(observed, nrow = n_sites, ncol = n_times, byrow = FALSE)
  pred_matrix <- matrix(predicted, nrow = n_sites, ncol = n_times, byrow = FALSE)
  
  residuals_avg <- rowMeans(obs_matrix - pred_matrix)
  
  df <- data.frame(
    lon = coords[, 1],
    lat = coords[, 2],
    residual = residuals_avg
  )
  
  ggplot(df, aes(x = lon, y = lat, fill = residual)) +
    geom_tile(width = 0.05, height = 0.05) +
    geom_point(shape = 21, size = 4, color = "black", stroke = 1) +
    scale_fill_gradient2(
      low = "#2E86AB", high = "#A23B72", mid = "white", 
      midpoint = 0,
      name = "Avg\nResidual"
    ) +
    labs(title = title, 
         x = "Longitude", y = "Latitude") +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right"
    )
}

plot_time_series_sites <- function(observed, predicted, n_sites, n_times, 
                                   site_indices = c(1, 5, 10),
                                   title = "Time Series at Selected Sites") {
  obs_matrix <- matrix(observed, nrow = n_sites, ncol = n_times, byrow = FALSE)
  pred_matrix <- matrix(predicted, nrow = n_sites, ncol = n_times, byrow = FALSE)
  
  df_list <- lapply(site_indices, function(i) {
    data.frame(
      Time = 1:n_times,
      Observation = obs_matrix[i, ],
      Prediction = pred_matrix[i, ],
      Site = paste("Site", i)
    )
  })
  
  df <- do.call(rbind, df_list)
  df_melt <- melt(df, id.vars = c("Time", "Site"))
  
  ggplot(df_melt, aes(x = Time, y = value, color = variable, linetype = variable)) +
    geom_line(size = 1.2) +
    geom_point(size = 2, alpha = 0.6) +
    facet_wrap(~Site, ncol = 1, scales = "free_y") +
    scale_color_manual(
      values = c("Observation" = "#E63946", "Prediction" = "#457B9D"),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c("Observation" = "solid", "Prediction" = "dashed"),
      name = NULL
    ) +
    labs(title = title, x = "Time", y = "Value") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    )
}

# ============================================================
# NEW: Safe parameter extraction and plotting
# ============================================================

plot_param_density_beta <- function(model_grid, model_random) {
  cat("    Checking beta availability...\n")
  
  # Check if beta exists
  has_beta_grid <- !is.null(model_grid$output) && 
    !is.null(model_grid$output$beta) && 
    length(model_grid$output$beta) > 0
  
  has_beta_random <- !is.null(model_random$output) && 
    !is.null(model_random$output$beta) && 
    length(model_random$output$beta) > 0
  
  if (!has_beta_grid || !has_beta_random) {
    cat("    ⚠ Beta coefficients not available in model output\n")
    cat("    Skipping beta plot\n")
    return(NULL)
  }
  
  # Extract beta
  beta_grid <- as.data.frame(model_grid$output$beta)
  beta_random <- as.data.frame(model_random$output$beta)
  
  if (nrow(beta_grid) == 0 || nrow(beta_random) == 0) {
    cat("    ⚠ Beta data frame is empty\n")
    return(NULL)
  }
  
  cat("    ✓ Beta found: Grid (", nrow(beta_grid), " samples), ",
      "Random (", nrow(beta_random), " samples)\n")
  
  beta_grid$Model <- "Grid"
  beta_random$Model <- "Random"
  
  beta_combined <- rbind(beta_grid, beta_random)
  beta_melt <- melt(beta_combined, id.vars = "Model")
  
  ggplot(beta_melt, aes(x = value, fill = Model)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~variable, scales = "free", ncol = 2) +
    scale_fill_manual(values = c("Grid" = "#6C5CE7", "Random" = "#00B894")) +
    labs(
      title = "Posterior Density of Beta Coefficients",
      subtitle = "Comparison between Grid and Random sampling",
      x = "Value", y = "Density"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

plot_param_density_phi <- function(model_grid, model_random) {
  cat("    Checking phi availability...\n")
  
  # Check if phi exists
  has_phi_grid <- !is.null(model_grid$output) && 
    !is.null(model_grid$output$phi) && 
    length(model_grid$output$phi) > 0
  
  has_phi_random <- !is.null(model_random$output) && 
    !is.null(model_random$output$phi) && 
    length(model_random$output$phi) > 0
  
  if (!has_phi_grid || !has_phi_random) {
    cat("    ⚠ Phi (spatial decay) not available in model output\n")
    cat("    Skipping phi plot\n")
    return(NULL)
  }
  
  phi_grid <- model_grid$output$phi
  phi_random <- model_random$output$phi
  
  cat("    ✓ Phi found: Grid (", length(phi_grid), " samples), ",
      "Random (", length(phi_random), " samples)\n")
  
  df <- data.frame(
    Phi = c(phi_grid, phi_random),
    Model = rep(c("Grid", "Random"), 
                c(length(phi_grid), length(phi_random)))
  )
  
  ggplot(df, aes(x = Phi, fill = Model)) +
    geom_density(alpha = 0.6) +
    geom_vline(aes(xintercept = mean(phi_grid), color = "Grid"), 
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = mean(phi_random), color = "Random"), 
               linetype = "dashed", size = 1) +
    scale_fill_manual(values = c("Grid" = "#6C5CE7", "Random" = "#00B894")) +
    scale_color_manual(values = c("Grid" = "#6C5CE7", "Random" = "#00B894"),
                       name = "Mean") +
    labs(
      title = "Posterior Density of Spatial Decay (φ)",
      subtitle = sprintf("Mean: Grid = %.3f | Random = %.3f", 
                         mean(phi_grid), mean(phi_random)),
      x = "Phi (Spatial Decay)", y = "Density"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom"
    )
}

# ============================================================
# Main Analysis
# ============================================================

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║   Grid vs Random Sampling (Full Diagnostics)         ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("Step 1: Generating base data to determine dimensions...\n")
set.seed(123)
temp_base <- data.sim.GP(sig2eps = 1, sig2eta = 1)

n_sites <- nrow(temp_base$coords)
n_times <- ncol(temp_base$Y)

cat("  Detected dimensions:", n_sites, "sites ×", n_times, "times\n\n")

cat("Step 2: Generating spatial sampling designs...\n")
grid_coords <- generate_grid_coords(n_sites)
random_coords <- generate_random_coords(n_sites, seed = 456)

cat("  Grid:   ", nrow(grid_coords), "sites\n")
cat("  Random: ", nrow(random_coords), "sites\n\n")

cat("Step 3: Creating datasets with different sampling designs...\n")

set.seed(123)
temp_grid <- data.sim.GP(sig2eps = 1, sig2eta = 1)

sim_data_grid <- list()
sim_data_grid$coords <- as.matrix(grid_coords)
sim_data_grid$X <- temp_grid$X
sim_data_grid$Y <- temp_grid$Y

cat("  Grid data: ", nrow(sim_data_grid$Y), "×", ncol(sim_data_grid$Y), "\n")

set.seed(789)
temp_random <- data.sim.GP(sig2eps = 1, sig2eta = 1)

sim_data_random <- list()
sim_data_random$coords <- as.matrix(random_coords)
sim_data_random$X <- temp_random$X
sim_data_random$Y <- temp_random$Y

cat("  Random data:", nrow(sim_data_random$Y), "×", ncol(sim_data_random$Y), "\n\n")

cat("✓ Dimension verification passed\n\n")

cat("Step 4: Fitting GP models...\n")

cat("  [1/2] Grid sampling...\n")
model_grid <- fit_model("GP", data.sim = sim_data_grid, nItr = 500, nBurn = 200)

cat("\n  [2/2] Random sampling...\n")
model_random <- fit_model("GP", data.sim = sim_data_random, nItr = 500, nBurn = 200)

cat("\n✓ Both models fitted\n\n")

cat("Step 5: Extracting predictions...\n")

obs_grid <- as.vector(t(sim_data_grid$Y))
obs_random <- as.vector(t(sim_data_random$Y))

pred_grid <- extract_predictions_safe(model_grid, sim_data_grid, "Grid")
pred_random <- extract_predictions_safe(model_random, sim_data_random, "Random")

cat("\n")

cat("Step 6: Computing performance metrics...\n")

metrics_grid <- compute_metrics(obs_grid, pred_grid)
metrics_grid$Model <- "Grid Sites"

metrics_random <- compute_metrics(obs_random, pred_random)
metrics_random$Model <- "Random Sites"

metrics_df <- rbind(metrics_grid, metrics_random)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║                    Results                            ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

print(metrics_df, row.names = FALSE, digits = 4)

mae_diff <- (metrics_grid$MAE - metrics_random$MAE) / metrics_random$MAE * 100
rmse_diff <- (metrics_grid$RMSE - metrics_random$RMSE) / metrics_random$RMSE * 100

cat("\nPerformance Difference (Grid vs Random):\n")
cat("  MAE:  ", sprintf("%+.2f%%", mae_diff), "\n")
cat("  RMSE: ", sprintf("%+.2f%%", rmse_diff), "\n\n")

cat("Step 7: Creating visualizations...\n\n")

cat("  [1/10] Performance metrics...\n")
p1 <- plot_metrics_comparison(metrics_df)
print(p1)
ggsave("grid_vs_random_metrics.png", p1, width = 14, height = 5.5, dpi = 300, bg = "white")
cat("      ✓ Saved\n\n")

cat("  [2/10] Spatial sampling designs...\n")
p2 <- plot_site_locations(grid_coords, random_coords)
print(p2)
ggsave("grid_vs_random_locations.png", p2, width = 12, height = 6, dpi = 300, bg = "white")
cat("      ✓ Saved\n\n")

cat("  [3/10] Observed vs Predicted...\n")
p3 <- plot_obs_vs_pred(obs_grid, pred_grid, obs_random, pred_random)
print(p3)
ggsave("grid_vs_random_predictions.png", p3, width = 12, height = 6, dpi = 300, bg = "white")
cat("      ✓ Saved\n\n")

cat("  [4/10] Distance structure comparison...\n")
png("grid_vs_random_distances.png", width = 2200, height = 950, res = 180)
plot_distance_comparison(grid_coords, random_coords)
dev.off()
cat("      ✓ Saved\n\n")

cat("  [5/10] Residual surface (Grid)...\n")
p5 <- plot_residual_surface(obs_grid, pred_grid, grid_coords, n_sites, n_times,
                            title = "Grid Sampling: Residual Surface")
print(p5)
ggsave("grid_residual_surface.png", p5, width = 10, height = 8, dpi = 300, bg = "white")
cat("      ✓ Saved\n\n")

cat("  [6/10] Residual surface (Random)...\n")
p6 <- plot_residual_surface(obs_random, pred_random, random_coords, n_sites, n_times,
                            title = "Random Sampling: Residual Surface")
print(p6)
ggsave("random_residual_surface.png", p6, width = 10, height = 8, dpi = 300, bg = "white")
cat("      ✓ Saved\n\n")

cat("  [7/10] Time series (Grid)...\n")
p7 <- plot_time_series_sites(obs_grid, pred_grid, n_sites, n_times,
                             site_indices = c(1, floor(n_sites/2), n_sites),
                             title = "Grid Sampling: Time Series")
print(p7)
ggsave("grid_time_series.png", p7, width = 12, height = 8, dpi = 300, bg = "white")
cat("      ✓ Saved\n\n")

cat("  [8/10] Time series (Random)...\n")
p8 <- plot_time_series_sites(obs_random, pred_random, n_sites, n_times,
                             site_indices = c(1, floor(n_sites/2), n_sites),
                             title = "Random Sampling: Time Series")
print(p8)
ggsave("random_time_series.png", p8, width = 12, height = 8, dpi = 300, bg = "white")
cat("      ✓ Saved\n\n")

# Beta and Phi plots with safe extraction
cat("  [9/10] Beta coefficient posteriors...\n")
p9 <- plot_param_density_beta(model_grid, model_random)
if (!is.null(p9)) {
  print(p9)
  ggsave("posterior_beta.png", p9, width = 12, height = 8, dpi = 300, bg = "white")
  cat("      ✓ Saved\n\n")
} else {
  cat("      ⊘ Skipped (data unavailable)\n\n")
}

cat("  [10/10] Spatial decay parameter (φ)...\n")
p10 <- plot_param_density_phi(model_grid, model_random)
if (!is.null(p10)) {
  print(p10)
  ggsave("posterior_phi.png", p10, width = 10, height = 6, dpi = 300, bg = "white")
  cat("      ✓ Saved\n\n")
} else {
  cat("      ⊘ Skipped (data unavailable)\n\n")
}

write.csv(metrics_df, "grid_vs_random_results.csv", row.names = FALSE)

cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║              Analysis Complete!                       ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

n_plots <- 8 + (!is.null(p9)) + (!is.null(p10))

cat("Generated", n_plots, "visualizations:\n")
cat("  1. grid_vs_random_metrics.png       - Performance comparison\n")
cat("  2. grid_vs_random_locations.png     - Spatial designs\n")
cat("  3. grid_vs_random_predictions.png   - Prediction accuracy\n")
cat("  4. grid_vs_random_distances.png     - Distance matrices\n")
cat("  5. grid_residual_surface.png        - Grid residual map\n")
cat("  6. random_residual_surface.png      - Random residual map\n")
cat("  7. grid_time_series.png             - Grid temporal dynamics\n")
cat("  8. random_time_series.png           - Random temporal dynamics\n")
if (!is.null(p9)) cat("  9. posterior_beta.png               - Beta posteriors\n")
if (!is.null(p10)) cat(" 10. posterior_phi.png                - Spatial decay posteriors\n")

cat("\nKey Findings:\n")
if (abs(mae_diff) < 5) {
  cat("  → Minimal difference between sampling designs (<5%)\n")
} else if (mae_diff < 0) {
  cat("  → Grid sampling shows", sprintf("%.1f%%", abs(mae_diff)), "better MAE\n")
} else {
  cat("  → Random sampling shows", sprintf("%.1f%%", mae_diff), "better MAE\n")
}

if (!is.null(p10)) {
  cat("\nParameter Estimates:\n")
  cat("  Grid φ:   ", sprintf("%.3f (SD: %.3f)", 
                              mean(model_grid$output$phi), 
                              sd(model_grid$output$phi)), "\n")
  cat("  Random φ: ", sprintf("%.3f (SD: %.3f)", 
                              mean(model_random$output$phi), 
                              sd(model_random$output$phi)), "\n")
}

cat("\n")