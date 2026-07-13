# ============================================================
# Plot 6 GWR local coefficient surfaces on Georgia map
# 2 x 3 layout, masked outside Georgia, with legends
# ============================================================
library(GWmodel)
library(sp)

data(GeorgiaCounties)

coef_names <- c(
  "Intercept",
  "PctRural",
  "PctEld",
  "PctFB",
  "PctPov",
  "PctBlack"
)

# 取得網格點所有不重複的XY座標並排序，image()需要規則的XY座標與對應矩陣
x_unique <- sort(unique(gwr_result$coef_estimates$x_coord))
y_unique <- sort(unique(gwr_result$coef_estimates$y_coord))

# 將座標轉為SpatialPoints空間點物件，指定與Georgia相同的座標參考系統
grid_points <- SpatialPoints(
  coords = gwr_result$coef_estimates[, c("x_coord", "y_coord")],
  proj4string = CRS(proj4string(Gedu.counties))
)

# 點落在州界內，使用該郡的屬性資料，點落在州界外，結果是NA
inside <- !is.na(over(grid_points, Gedu.counties)[, 1])

# 建立含有100種顏色的色階
pal <- heat.colors(6)

# 儲存原繪圖參數，避免程式執行後永久改變R的繪圖設定
old_par <- par(no.readonly = TRUE)

# 2 x 3 layout
par(
  mfrow = c(2, 3),
  mar = c(3, 3, 3, 5)
)

for (coef_name in coef_names) {
  # Mask coefficient values outside Georgia
  coef_value <- gwr_result$coef_estimates[[coef_name]]
  coef_value_masked <- coef_value
  coef_value_masked[!inside] <- NA
  
  # 向量轉矩陣 供image()使用
  # coef 原始資料排序需和網格座標排列一致
  z_masked <- matrix(
    coef_value_masked,
    nrow = length(x_unique),
    ncol = length(y_unique)
  )
  
  # 取得該系數最小值與最大值，切成100顏色區間
  z_range <- range(z_masked, na.rm = TRUE)
  breaks <- seq(z_range[1], z_range[2], length.out = length(pal) + 1)
  
  # Plot coefficient surface
  image(
    x_unique,
    y_unique,
    z_masked,
    col = pal,
    breaks = breaks,
    xlab = "X",
    ylab = "Y",
    main = coef_name,
    axes = TRUE
  )
  
  # Add Georgia boundaries
  plot(
    Gedu.counties,
    add = TRUE,
    border = "black",
    lwd = 0.5
  )
  
  # Add observation points
  points(
    obs_xy[, 1],
    obs_xy[, 2],
    pch = 16,
    cex = 0.15
  )
  
  # Add legend / color scale
  legend_labels <- paste0(
    round(breaks[-length(breaks)], 3),
    " – ",
    round(breaks[-1], 3)
  )
  
  legend(
    "right",
    inset = c(-0.25, 0),
    legend = legend_labels,
    fill = pal,
    title = "Coef.",
    cex = 0.65,
    bty = "n",
    xpd = TRUE
  )
}

# Restore original plotting parameters
par(old_par)
