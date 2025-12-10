# Nonseparable Test with Robust Prediction Extraction 
#YOU can change the address of source to yours
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

plot_model_metrics <- function(results_df) {
  results_long <- results_df %>%
    tidyr::pivot_longer(cols = c(MAE, RMSE, R2), names_to = "Metric", values_to = "Value")
  
  ggplot(results_long, aes(x = Model, y = Value, fill = Model)) +
    geom_bar(stat = "identity", width = 0.6) +
    facet_wrap(~Metric, scales = "free_y") +
    geom_text(aes(label = round(Value, 4)), vjust = -0.3, size = 3.5) +
    labs(title = "Model Performance: Separable vs Nonseparable Covariance", 
         y = "Metric Value") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Robust prediction extraction function
extract_predictions_robust <- function(model, sim_data) {
  n_sites <- nrow(sim_data$coords)
  n_times <- ncol(sim_data$Y)
  expected_length <- n_sites * n_times
  
  cat("\n  Extracting predictions:\n")
  cat("    Expected length:", expected_length, "(", n_sites, "sites x", n_times, "times)\n")
  
  # Check model$fitted
  if (!is.null(model$fitted)) {
    fitted_data <- model$fitted
    
    cat("    model$fitted class:", class(fitted_data), "\n")
    
    # Case 1: It's a matrix
    if (is.matrix(fitted_data)) {
      cat("    model$fitted is matrix:", nrow(fitted_data), "x", ncol(fitted_data), "\n")
      
      # Check if it matches our data structure
      if (nrow(fitted_data) == n_sites && ncol(fitted_data) == n_times) {
        # Perfect match: sites x times
        pred <- as.vector(t(fitted_data))  # Transpose and vectorize
        cat("    Converted matrix (sites x times) to vector\n")
        return(pred)
        
      } else if (nrow(fitted_data) == n_times && ncol(fitted_data) == n_sites) {
        # Transposed: times x sites
        pred <- as.vector(fitted_data)
        cat("    Converted matrix (times x sites) to vector\n")
        return(pred)
        
      } else if (ncol(fitted_data) == 1) {
        # Column vector
        pred <- as.vector(fitted_data)
        cat("    Converted column vector to vector, length:", length(pred), "\n")
        
        # Check if length matches
        if (length(pred) == expected_length) {
          return(pred)
        } else if (length(pred) == n_sites) {
          # Per-site predictions, replicate for each time
          pred <- rep(pred, times = n_times)
          cat("    Replicated site predictions for", n_times, "times\n")
          return(pred)
        }
      } else {
        # Unknown matrix structure
        cat("    Unknown matrix structure\n")
        pred <- as.vector(fitted_data)
        cat("    Vectorized to length:", length(pred), "\n")
        
        if (length(pred) == expected_length) {
          return(pred)
        }
      }
    } 
    # Case 2: It's already a vector
    else if (is.vector(fitted_data) || is.numeric(fitted_data)) {
      cat("    model$fitted is vector, length:", length(fitted_data), "\n")
      
      if (length(fitted_data) == expected_length) {
        cat("   Length matches!\n")
        return(fitted_data)
        
      } else if (length(fitted_data) == n_sites) {
        # Per-site predictions
        pred <- rep(fitted_data, times = n_times)
        cat("    Replicated", n_sites, "site predictions for", n_times, "times\n")
        return(pred)
        
      } else if (length(fitted_data) == n_times) {
        # Per-time predictions (unusual)
        pred <- rep(fitted_data, each = n_sites)
        cat("    Replicated", n_times, "time predictions for", n_sites, "sites\n")
        return(pred)
      }
    }
    
    # Case 3: Array or other structure
    else {
      cat("    model$fitted is:", class(fitted_data), "\n")
      pred <- as.vector(fitted_data)
      cat("    Forced to vector, length:", length(pred), "\n")
      
      if (length(pred) == expected_length) {
        return(pred)
      }
    }
  }
  
  # Fallback: use posterior mean of parameters
  cat("    Using fallback: posterior mean of beta only\n")
  
  if (!is.null(model$betap)) {
    # betap is posterior samples of beta
    beta_mean <- colMeans(model$betap)
    X_matrix <- sim_data$X
    mu_pred <- as.vector(X_matrix %*% beta_mean)
    pred <- rep(mu_pred, each = n_times)
    cat("    Generated predictions from beta, length:", length(pred), "\n")
    return(pred)
  }
  
  stop("ERROR: Cannot extract predictions from model!")
}


# Step 1: Simulate separable data
set.seed(123)
sim_data_sep <- data.sim.GP(sig2eps = 1, sig2eta = 1)

coords <- sim_data_sep$coords
n_sites <- nrow(coords)
n_times <- ncol(sim_data_sep$Y)

cat("  Sites:", n_sites, ", Times:", n_times, "\n\n")

# Step 2: Create nonseparable data

time_points <- 1:n_times

cov_mat_nonsep <- nonsep_cov_fast(
  coords = coords,
  time = time_points,
  phi_s = 0.8,
  phi_t = 0.08,
  sig2 = 2
)

min_eig <- min(eigen(cov_mat_nonsep, symmetric = TRUE, only.values = TRUE)$values)

if (min_eig <= 0) {
  stop("ERROR: Covariance matrix is not positive definite!")
}

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

cat("  Separable Y:    mean =", round(mean(sim_data_sep$Y), 2), 
    ", sd =", round(sd(as.vector(sim_data_sep$Y)), 2), "\n")
cat("  Nonseparable Y: mean =", round(mean(Y_nonsep), 2), 
    ", sd =", round(sd(as.vector(Y_nonsep)), 2), "\n\n")

# Step 3: Fit models
model_sep <- fit_model("GP", data.sim = sim_data_sep, nItr = 500, nBurn = 200)
model_nonsep <- fit_model("GP", data.sim = sim_data_nonsep, nItr = 500, nBurn = 200)

# Step 4: Extract predictions (ROBUST)
cat("Step 4: Extracting predictions with robust method...\n")

obs_sep <- as.vector(t(sim_data_sep$Y))
obs_nonsep <- as.vector(t(sim_data_nonsep$Y))

cat("\nFor SEPARABLE model:")
pred_sep <- extract_predictions_robust(model_sep, sim_data_sep)

cat("\nFor NONSEPARABLE model:")
pred_nonsep <- extract_predictions_robust(model_nonsep, sim_data_nonsep)

# Verify lengths
cat("\n Prediction extraction complete\n")
cat("  Separable:    obs =", length(obs_sep), ", pred =", length(pred_sep), "\n")
cat("  Nonseparable: obs =", length(obs_nonsep), ", pred =", length(pred_nonsep), "\n\n")

# Step 5: Compute metrics
cat("Step 5: Computing metrics...\n")

metrics_sep <- compute_metrics(obs_sep, pred_sep)
metrics_sep$Model <- "GP on Separable Data"

metrics_nonsep <- compute_metrics(obs_nonsep, pred_nonsep)
metrics_nonsep$Model <- "GP on Nonseparable Data"

results_df <- rbind(metrics_sep, metrics_nonsep)

print(results_df, row.names = FALSE, digits = 4)

mae_change <- (metrics_nonsep$MAE - metrics_sep$MAE) / metrics_sep$MAE * 100
cat("\nPerformance change:", sprintf("%+.2f%%", mae_change), "MAE\n\n")

# Step 6: Visualizations

p1 <- plot_model_metrics(results_df)
ggsave("nonsep_metrics_comparison.png", plot = p1, width = 12, height = 5, dpi = 300)

comparison_df <- data.frame(
  Observed = c(obs_sep, obs_nonsep),
  Predicted = c(pred_sep, pred_nonsep),
  Model = rep(c("Separable\n(Correct)", "Nonseparable\n(Misspecified)"),
              c(length(obs_sep), length(obs_nonsep)))
)

p2 <- ggplot(comparison_df, aes(x = Observed, y = Predicted, color = Model)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~Model, ncol = 2) +
  labs(title = "Observed vs Predicted: Covariance Misspecification",
       x = "Observed", y = "Predicted") +
  theme_minimal() +
  theme(legend.position = "none")

print(p2)
ggsave("nonsep_obs_vs_pred.png", plot = p2, width = 10, height = 5, dpi = 300)

write.csv(results_df, "nonsep_test_results.csv", row.names = FALSE)

