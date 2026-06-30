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

# Unique grid coordinates
x_unique <- sort(unique(gwr_result$x_coord))
y_unique <- sort(unique(gwr_result$y_coord))

# Build SpatialPoints for grid
grid_points <- SpatialPoints(
  coords = gwr_result[, c("x_coord", "y_coord")],
  proj4string = CRS(proj4string(Gedu.counties))
)

# Check which grid points are inside Georgia
inside <- !is.na(over(grid_points, Gedu.counties)[, 1])

# Color palette
pal <- heat.colors(100)

# Save original plotting parameters
old_par <- par(no.readonly = TRUE)

# 2 x 3 layout
par(
  mfrow = c(2, 3),
  mar = c(3, 3, 3, 5)
)

for (coef_name in coef_names) {
  
  # Mask coefficient values outside Georgia
  coef_value <- gwr_result[[coef_name]]
  coef_value_masked <- coef_value
  coef_value_masked[!inside] <- NA
  
  # Convert to matrix for image()
  z_masked <- matrix(
    coef_value_masked,
    nrow = length(x_unique),
    ncol = length(y_unique)
  )
  
  # Breaks for color scale
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
  legend_values <- pretty(z_range, n = 5)
  legend_colors <- pal[
    findInterval(
      legend_values,
      breaks,
      all.inside = TRUE
    )
  ]
  
  legend(
    "right",
    inset = c(-0.25, 0),
    legend = round(legend_values, 3),
    fill = legend_colors,
    title = "Coef.",
    cex = 0.65,
    bty = "n",
    xpd = TRUE
  )
}

# Restore original plotting parameters
par(old_par)
