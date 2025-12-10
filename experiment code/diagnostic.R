# diagnostic.R
setwd('D:/python project')
source("D:/python project/data_simulation.R")
library('spTimer')

# 模拟数据
sim_data <- data.sim.GP(sig2eps = 1, sig2eta = 1)

# 手动构建数据框
n <- nrow(sim_data$Y)
time_length <- ncol(sim_data$Y)

df <- data.frame(
  site = rep(1:n, each = time_length),
  time = rep(1:time_length, times = n),
  Y = as.vector(t(sim_data$Y)),
  x1 = rep(sim_data$X[, 2], times = n),
  x2 = rep(sim_data$X[, 3], times = n),
  x3 = rep(sim_data$X[, 4], times = n),
  lon = rep(sim_data$coords[, 1], each = time_length),
  lat = rep(sim_data$coords[, 2], each = time_length)
)

cat("Data frame created successfully\n")
print(str(df))
print(summary(df))

# 检查是否有 NA 或 Inf
cat("\nChecking for NA/Inf values:\n")
print(sapply(df, function(x) sum(is.na(x))))
print(sapply(df, function(x) sum(is.infinite(x))))

# 尝试拟合
cat("\nAttempting to fit model...\n")
tryCatch({
  model_out <- spT.Gibbs(
    formula = Y ~ x1 + x2 + x3,
    data = df,
    model = "GP",
    coords = ~lon + lat,
    nItr = 100,
    nBurn = 50,
    report = 10
  )
  cat("SUCCESS!\n")
  print(summary(model_out))
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  print(traceback())
})



setwd('D:/python project')
source("D:/python project/data_simulation.R")

set.seed(123)
sim_data <- data.sim.GP(sig2eps = 1, sig2eta = 1)

coords <- sim_data$coords
cat("=== Coordinate Information ===\n")
cat("Class:", class(coords), "\n")
cat("Dimensions:", dim(coords), "\n")
cat("Column names:", colnames(coords), "\n\n")

print("Coordinates:")
print(coords)

cat("\nSummary:\n")
print(summary(coords))

# 检查是否有重复点
cat("\nDuplicate points:", any(duplicated(coords)), "\n")

# 检查坐标范围
cat("\nCoordinate ranges:\n")
cat("  X:", range(coords[,1]), "\n")
cat("  Y:", range(coords[,2]), "\n")

# 可视化
plot(coords, pch = 19, col = "blue", 
     main = "Simulated Site Locations",
     xlab = "Longitude", ylab = "Latitude")
text(coords, labels = 1:nrow(coords), pos = 3, cex = 0.7)

model_test <- fit_model("GP", data.sim = sim_data, nItr = 100, nBurn = 50)

cat("Model slots:\n")
print(names(model_test))

cat("\nChecking prediction components:\n")
if (!is.null(model_test$fitted.values)) {
  cat("  fitted.values: length =", length(model_test$fitted.values), "\n")
  cat("  range:", range(model_test$fitted.values), "\n")
}

if (!is.null(model_test$prediction)) {
  cat("  prediction structure:\n")
  print(str(model_test$prediction))
}

if (!is.null(model_test$mcmc.beta)) {
  cat("  mcmc.beta: dim =", dim(model_test$mcmc.beta), "\n")
}