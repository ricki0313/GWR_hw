# =====================================================
# cv vs. gcv 
# run basic_estimates_maual.R before running these codes
# =====================================================
library(bench)

# visualization
cv_result <- select_bw_cv(
  X = X,
  y = y,
  obs = obs_xy,
  bw_candidates = bw_candidates
)

gcv_result <- select_bw_gcv(
  X = X,
  y = y,
  obs = obs_xy,
  bw_candidates = bw_candidates
)

par(mfrow = c(1, 2))

# CV
plot(
  plot_data$bw,
  plot_data$CV,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "Bandwidth",
  ylab = "CV value",
  main = "CV Criterion"
)

points(
  cv_result$best$bw,
  cv_result$best$cv,
  col = "blue",
  pch = 16
)

text(
  cv_result$best$bw,
  cv_result$best$cv,
  labels = paste0(
    "min = ", round(cv_result$best$cv, 2),
    "\nbw = ", cv_result$best$bw
  ),
  pos = 3,
  col = "blue",
  cex = 0.8
)

# GCV
plot(
  plot_data$bw,
  plot_data$GCV,
  type = "l",
  col = "red",
  lwd = 2,
  xlab = "Bandwidth",
  ylab = "GCV value",
  main = "GCV Criterion"
)

points(
  gcv_result$best$bw,
  gcv_result$best$gcv,
  col = "red",
  pch = 16
)

text(
  gcv_result$best$bw,
  gcv_result$best$gcv,
  labels = paste0(
    "min = ", round(gcv_result$best$gcv, 2),
    "\nbw = ", gcv_result$best$bw
  ),
  pos = 3,
  col = "red",
  cex = 0.8
)


# efficacy
benchmark_result <- bench::mark(
  CV = select_bw_cv(
    X = X,
    y = y,
    obs = obs_xy,
    bw_candidates = bw_candidates
  ),
  
  GCV = select_bw_gcv(
    X = X,
    y = y,
    obs = obs_xy,
    bw_candidates = bw_candidates
  ),
  
  iterations = 20,
  check = FALSE
)
benchmark_result