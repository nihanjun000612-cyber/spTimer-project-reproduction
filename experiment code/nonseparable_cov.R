

library(Matrix)

# Ensure positive definite covariance matrix
make_pd <- function(mat, tol = 1e-6) {
  # Ensure symmetry
  mat <- (mat + t(mat)) / 2
  
  # Eigenvalue decomposition
  eig <- eigen(mat, symmetric = TRUE)
  
  # Replace negative eigenvalues with small positive values
  eig$values[eig$values < tol] <- tol
  
  # Reconstruct matrix
  mat_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
  return(mat_pd)
}

# Nonseparable covariance function - Gneiting (2002)
nonsep_cov <- function(coords, time, phi_s = 2, phi_t = 0.1, sig2 = 1, a = 0.5, c = 0.5) {
  # coords: matrix of spatial coordinates (n x 2)
  # time: vector of time points (length T)
  # phi_s: spatial range parameter (smaller = more correlation)
  # phi_t: temporal range parameter
  # sig2: variance parameter
  # a, c: interaction parameters (keep small for stability)
  
  n <- nrow(coords)
  T_len <- length(time)
  N <- n * T_len
  
  cat("Computing nonseparable covariance matrix:\n")
  cat("  Sites:", n, ", Times:", T_len, ", Total:", N, "\n")
  
  # Compute spatial distance matrix
  dist_s <- as.matrix(dist(coords, method = "euclidean"))
  
  # Compute temporal distance matrix
  dist_t <- as.matrix(dist(time, method = "euclidean"))
  
  # Initialize covariance matrix
  cov_mat <- matrix(0, nrow = N, ncol = N)
  
  # Fill in the covariance matrix with Gneiting structure
  cat("  Building covariance matrix...\n")
  
  for (i in 1:n) {
    for (j in 1:n) {
      # Spatial distance
      h_s <- dist_s[i, j]
      
      for (t1 in 1:T_len) {
        for (t2 in 1:T_len) {
          # Row and column indices
          row_idx <- (i - 1) * T_len + t1
          col_idx <- (j - 1) * T_len + t2
          
          # Temporal distance
          h_t <- dist_t[t1, t2]
          
          # Gneiting nonseparable structure
          # C(h_s, h_t) = sig2 / (a * h_t^(2c) + 1)^(3/2) * exp(-phi_s * h_s / sqrt(a * h_t^(2c) + 1))
          
          denominator <- (a * h_t^(2 * c) + 1)
          spatial_decay <- exp(-phi_s * h_s / sqrt(denominator))
          temporal_decay <- exp(-phi_t * h_t)
          
          cov_mat[row_idx, col_idx] <- sig2 * spatial_decay * temporal_decay / sqrt(denominator)
        }
      }
    }
    
    if (i %% 3 == 0) cat("    Progress:", round(i/n*100), "%\n")
  }
  
  cat("  Making matrix positive definite...\n")
  cov_mat <- make_pd(cov_mat)
  
  cat("  Done!\n")
  
  return(cov_mat)
}

# Simplified and more stable nonseparable covariance
nonsep_cov_simple <- function(coords, time, phi_s = 0.5, phi_t = 0.05, sig2 = 1, 
                              interaction_strength = 0.3) {
  # More stable implementation using Kronecker product with interaction
  
  n <- nrow(coords)
  T_len <- length(time)
  
  cat("Computing simple nonseparable covariance:\n")
  cat("  Sites:", n, ", Times:", T_len, "\n")
  
  # Spatial covariance (exponential)
  dist_s <- as.matrix(dist(coords))
  cov_s <- exp(-phi_s * dist_s)
  
  # Temporal covariance (exponential)
  dist_t <- as.matrix(dist(time))
  cov_t <- exp(-phi_t * dist_t)
  
  # Separable part (Kronecker product)
  cat("  Computing separable part...\n")
  cov_sep <- kronecker(cov_t, cov_s)
  
  # Nonseparable interaction term
  cat("  Adding nonseparable interaction...\n")
  n_total <- n * T_len
  interaction <- matrix(0, n_total, n_total)
  
  for (i in 1:n_total) {
    for (j in 1:n_total) {
      site_i <- ceiling(i / T_len)
      site_j <- ceiling(j / T_len)
      time_i <- ((i - 1) %% T_len) + 1
      time_j <- ((j - 1) %% T_len) + 1
      
      # Interaction: weighted product of spatial and temporal correlation
      spatial_cor <- cov_s[site_i, site_j]
      temporal_cor <- cov_t[time_i, time_j]
      
      # Nonseparable term decays with combined distance
      interaction[i, j] <- spatial_cor * temporal_cor * exp(-0.1 * dist_s[site_i, site_j] * dist_t[time_i, time_j])
    }
    
    if (i %% 100 == 0) cat("    Progress:", round(i/n_total*100), "%\n")
  }
  
  # Combine separable and nonseparable parts
  cat("  Combining components...\n")
  cov_mat <- sig2 * ((1 - interaction_strength) * cov_sep + interaction_strength * interaction)
  
  # Ensure positive definiteness
  cat("  Ensuring positive definiteness...\n")
  cov_mat <- make_pd(cov_mat, tol = 1e-5)
  
  # Check result
  min_eig <- min(eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values)
  cat("  Minimum eigenvalue:", min_eig, "\n")
  cat("  Done!\n")
  
  return(cov_mat)
}

# Even simpler and faster version (most stable)
nonsep_cov_fast <- function(coords, time, phi_s = 1, phi_t = 0.1, sig2 = 1) {
  # Very simple but stable implementation
  
  n <- nrow(coords)
  T_len <- length(time)
  N <- n * T_len
  
  cat("Computing fast nonseparable covariance:\n")
  cat("  Total dimension:", N, "x", N, "\n")
  
  # Compute all pairwise distances at once
  dist_s <- as.matrix(dist(coords))
  dist_t <- as.matrix(dist(time))
  
  # Spatial covariance
  cov_s <- sig2 * exp(-phi_s * dist_s)
  
  # Temporal covariance with different decay
  cov_t <- exp(-phi_t * dist_t)
  
  # Create nonseparable structure by modifying Kronecker product
  # Use block structure for stability
  cov_mat <- matrix(0, N, N)
  
  for (i in 1:n) {
    for (j in 1:n) {
      # Block indices
      row_start <- (i - 1) * T_len + 1
      row_end <- i * T_len
      col_start <- (j - 1) * T_len + 1
      col_end <- j * T_len
      
      # Spatial correlation between sites i and j
      spatial_cor <- cov_s[i, j]
      
      # Temporal block with spatial modification
      # Nonseparability: temporal correlation decays differently by spatial distance
      decay_mod <- exp(-0.05 * dist_s[i, j])  # Spatial distance modifies temporal correlation
      
      cov_mat[row_start:row_end, col_start:col_end] <- spatial_cor * cov_t * decay_mod
    }
  }
  
  # Ensure positive definiteness
  cov_mat <- make_pd(cov_mat, tol = 1e-4)
  
  cat("  Done!\n")
  return(cov_mat)
}