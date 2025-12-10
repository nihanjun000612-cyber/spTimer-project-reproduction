# Model fitting function - spTimer (GP, AR, GPP)

library(spTimer)

# General model fitting function
fit_model <- function(model_type = "GP", data.sim, nItr = 200, nBurn = 100, n_knots = 4) {
  coords <- data.sim$coords
  time_length <- ncol(data.sim$Y)
  n <- nrow(data.sim$Y)
  
  # Ensure coords are in the correct format
  if (!is.matrix(coords)) {
    coords <- as.matrix(coords)
  }
  if (is.null(colnames(coords))) {
    colnames(coords) <- c("coords.Var1", "coords.Var2")
  }
  
  # Build data frame
  df <- data.frame(
    site = rep(1:n, each = time_length),
    time = rep(1:time_length, times = n),
    Y = as.vector(t(data.sim$Y)),
    x1 = rep(data.sim$X[, 2], times = n),
    x2 = rep(data.sim$X[, 3], times = n),
    x3 = rep(data.sim$X[, 4], times = n),
    lon = rep(coords[, 1], each = time_length),
    lat = rep(coords[, 2], each = time_length)
  )
  
  if (model_type == "GPP") {
    # Method: Generate uniform grid knots within the observation range, 
    # slightly shrunken inward to ensure they are inside the convex hull
    lon_range <- range(coords[, 1])
    lat_range <- range(coords[, 2])
    
    # Shrink inward by 15% margin
    lon_margin <- diff(lon_range) * 0.15
    lat_margin <- diff(lat_range) * 0.15
    
    # Calculate grid dimension
    grid_dim <- ceiling(sqrt(n_knots))
    
    # Generate uniform grid
    lon_seq <- seq(lon_range[1] + lon_margin, 
                   lon_range[2] - lon_margin, 
                   length.out = grid_dim)
    lat_seq <- seq(lat_range[1] + lat_margin, 
                   lat_range[2] - lat_margin, 
                   length.out = grid_dim)
    
    knots_grid <- expand.grid(lon = lon_seq, lat = lat_seq)
    knots_coords <- as.matrix(knots_grid)
    colnames(knots_coords) <- c("coords.Var1", "coords.Var2")
    
    # Limit to the specified number of knots
    if (nrow(knots_coords) > n_knots) {
      knots_coords <- knots_coords[1:n_knots, ]
    }
    
    cat("\n=== GPP Model Knots Configuration ===\n")
    cat("Number of knots:", nrow(knots_coords), "\n")
    cat("Sites range: lon [", lon_range[1], ",", lon_range[2], 
        "], lat [", lat_range[1], ",", lat_range[2], "]\n")
    cat("Knots range: lon [", min(knots_coords[,1]), ",", max(knots_coords[,1]), 
        "], lat [", min(knots_coords[,2]), ",", max(knots_coords[,2]), "]\n")
    cat("Knots coordinates:\n")
    print(knots_coords)
    
    # Visualize knots and sites
    tryCatch({
      plot(coords, pch = 19, col = "blue", cex = 1.5,
           main = "Sites (blue) and Knots (red)",
           xlab = "Longitude", ylab = "Latitude",
           xlim = range(c(coords[,1], knots_coords[,1])),
           ylim = range(c(coords[,2], knots_coords[,2])))
      points(knots_coords, pch = 17, col = "red", cex = 2)
      legend("topright", legend = c("Sites", "Knots"), 
             pch = c(19, 17), col = c("blue", "red"), cex = 1.2)
    }, error = function(e) {
      cat("Could not generate visualization\n")
    })
    
    # Fit model
    model_out <- spT.Gibbs(
      formula = Y ~ x1 + x2 + x3,
      data = df,
      model = model_type,
      coords = ~lon + lat,
      knots.coords = knots_coords,
      nItr = nItr,
      nBurn = nBurn,
      report = 10,
      tol.dist = 0.1  # Increase distance tolerance
    )
    
  } else {
    # For GP or AR models
    model_out <- spT.Gibbs(
      formula = Y ~ x1 + x2 + x3,
      data = df,
      model = model_type,
      coords = ~lon + lat,
      nItr = nItr,
      nBurn = nBurn,
      report = 10
    )
  }
  
  return(model_out)
}