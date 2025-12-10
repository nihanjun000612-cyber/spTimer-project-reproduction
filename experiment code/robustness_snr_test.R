
# 鲁棒性分析：在不同 SNR（信噪比）下模拟数据并拟合 GP 模型

source("R/functions/data_simulation.R")
source("R/functions/model_fitting.R")

# 设置信噪比测试值
snr_values <- c(1, 10, 15)
results <- list()

for (snr in snr_values) {
  cat("SNR =", snr, "\n")

  # 设置信号与噪声的方差
  sig2eta <- 1
  sig2eps <- sig2eta / snr

  # 模拟数据
  sim_data <- data.sim.GP(sig2eps = sig2eps, sig2eta = sig2eta)

  # 模型拟合
  model_result <- fit_model(model_type = "GP", data.sim = sim_data, nItr = 200, nBurn = 100)

  # 保存结果
  results[[paste0("SNR_", snr)]] <- model_result
}
# 提示：可以用 spT.validation() 对 results 进行性能评估
