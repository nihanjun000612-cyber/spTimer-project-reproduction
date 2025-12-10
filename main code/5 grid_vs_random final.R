#  Grid vs Random Sampling Comparison (FIXED)
# You can change the address to yours
library(MASS)
library(spTimer)
library(ggplot2)
library(dplyr)

setwd('D:/python project')
source("D:/python project/data_simulation.R")
source("D:/python project/model_fitting.R")

# Helper Functions
compute_metrics <- function(observed, predicted) {
  min_len <- min(length(observed), length(predicted))
  observed <- observed[1:min_len]
  predicted <- predicted[1:min_len]
  
  mae <- mean(abs(observed - predicted))
  rmse <- sqrt(mean((observed - predicted)^2))
  r2 <- 1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
  
  return(data.frame(MAE = mae, RMSE = rmse, R2 = r2))
}

# Generate grid coordinates
generate_grid_coords <- function(n_sites) {
  # Find closest square grid
  n_side <- ceiling(sqrt(n_sites))
  grid <- expand.grid(
    lon = seq(0, 1, length.out = n_side),
    lat = seq(0, 1, length.out = n_side)
  )
  # Take only first n_sites
  return(grid[1:n_sites, ])
}

# Generate random coordinates
generate_random_coords <- function(n_sites, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data.frame(
    lon = runif(n_sites, 0, 1),
    lat = runif(n_sites, 0, 1)
  )
}

# Extract predictions safely
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

# Visualization Functions

plot_metrics_comparison <- function(metrics_df) {
  metrics_long <- metrics_df %>%
    tidyr::pivot_longer(cols = c(MAE, RMSE, R2), 
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

# Main Analysis
# Step 1: Simulate base data first (to get dimensions)
set.seed(123)
sim_data_base <- data.sim.GP(sig2eps = 1, sig2eta = 1)

n_sites <- nrow(sim_data_base$coords)
n_times <- ncol(sim_data_base$Y)

cat("  Data dimensions: ", n_sites, "sites ×", n_times, "times\n")
cat("  Total observations:", n_sites * n_times, "\n\n")

# Step 2: Generate matching coordinates

# Grid coordinates (matching n_sites)
grid_coords <- generate_grid_coords(n_sites)

# Random coordinates (matching n_sites)
random_coords <- generate_random_coords(n_sites, seed = 456)

cat("  Grid:   ", nrow(grid_coords), "sites in regular pattern\n")
cat("  Random: ", nrow(random_coords), "sites uniformly distributed\n\n")

# Step 3: Create two datasets with different coordinates

# Grid version
sim_data_grid <- sim_data_base
sim_data_grid$coords <- as.matrix(grid_coords)

# Random version
sim_data_random <- sim_data_base
sim_data_random$coords <- as.matrix(random_coords)

cat("  ✓ Grid data prepared\n")
cat("  ✓ Random data prepared\n\n")

# Verify dimensions
cat("Verification:\n")
cat("  Grid coords:   ", nrow(sim_data_grid$coords), "rows\n")
cat("  Grid Y matrix: ", nrow(sim_data_grid$Y), "×", ncol(sim_data_grid$Y), "\n")
cat("  Random coords: ", nrow(sim_data_random$coords), "rows\n")
cat("  Random Y matrix:", nrow(sim_data_random$Y), "×", ncol(sim_data_random$Y), "\n\n")

# Step 4: Fit GP models
model_grid <- fit_model("GP", data.sim = sim_data_grid, nItr = 500, nBurn = 200)
model_random <- fit_model("GP", data.sim = sim_data_random, nItr = 500, nBurn = 200)

# Step 5: Extract predictions

obs_grid <- as.vector(t(sim_data_grid$Y))
obs_random <- as.vector(t(sim_data_random$Y))

pred_grid <- extract_predictions_safe(model_grid, sim_data_grid, "Grid")
pred_random <- extract_predictions_safe(model_random, sim_data_random, "Random")

# Step 6: Compute metrics


metrics_grid <- compute_metrics(obs_grid, pred_grid)
metrics_grid$Model <- "Grid Sites"

metrics_random <- compute_metrics(obs_random, pred_random)
metrics_random$Model <- "Random Sites"

metrics_df <- rbind(metrics_grid, metrics_random)

print(metrics_df, row.names = FALSE, digits = 4)

# Performance comparison
mae_diff <- (metrics_grid$MAE - metrics_random$MAE) / metrics_random$MAE * 100
rmse_diff <- (metrics_grid$RMSE - metrics_random$RMSE) / metrics_random$RMSE * 100

cat("\nPerformance Difference (Grid vs Random):\n")
cat("  MAE:  ", sprintf("%+.2f%%", mae_diff), "\n")
cat("  RMSE: ", sprintf("%+.2f%%", rmse_diff), "\n\n")

# Step 7: Visualizations

# Plot 1: Metrics
p1 <- plot_metrics_comparison(metrics_df)
print(p1)
ggsave("grid_vs_random_metrics.png", p1, width = 14, height = 5.5, dpi = 300, bg = "white")

# Plot 2: Site locations
p2 <- plot_site_locations(grid_coords, random_coords)
print(p2)
ggsave("grid_vs_random_locations.png", p2, width = 12, height = 6, dpi = 300, bg = "white")

# Plot 3: Obs vs Pred
p3 <- plot_obs_vs_pred(obs_grid, pred_grid, obs_random, pred_random)
print(p3)
ggsave("grid_vs_random_predictions.png", p3, width = 12, height = 6, dpi = 300, bg = "white")

# Plot 4: Distance matrices
png("grid_vs_random_distances.png", width = 2200, height = 950, res = 180)
plot_distance_comparison(grid_coords, random_coords)
dev.off()

# Save results
write.csv(metrics_df, "grid_vs_random_results.csv", row.names = FALSE)


# Spatial coverage analysis
grid_dists <- as.vector(dist(grid_coords))
random_dists <- as.vector(dist(random_coords))

cat("Distance Statistics:\n")
cat("  Grid -   Mean:", round(mean(grid_dists), 3), 
    "| SD:", round(sd(grid_dists), 3), "\n")
cat("  Random - Mean:", round(mean(random_dists), 3), 
    "| SD:", round(sd(random_dists), 3), "\n\n")

if (sd(grid_dists) < sd(random_dists)) {
  cat("  → Grid has more uniform spacing (lower SD)\n")
} else {
  cat("  → Random has more variable spacing\n")
}

cat("\n")