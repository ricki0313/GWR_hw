# install.packages("GWmodel")
# install.packages("sp")

library(GWmodel)
library(sp)

# ----- init. -----
# Georgia census data set (GWmodel doc. p.16) 
# transform data set to spdf
data(Georgia) 
coords <- cbind(Gedu.df$X, Gedu.df$Y)
georgia_spdf <- SpatialPointsDataFrame(
  coords = coords,
  data = Gedu.df
)

names(georgia_spdf)
head(georgia_spdf)
# ----- variable & bw selection -----
# choose variable (GWR ref. p.97 table 4.5)
formula_gwr <- PctBach ~ PctRural + PctEld + PctFB + PctPov + PctBlack

# find the best bandwidth (CV, AICc)
# bw_cv <- bw.gwr(
#   formula = formula_gwr,
#   data = georgia_spdf,
#   approach = "CV",
#   kernel = "bisquare",
#   adaptive = TRUE
# )
# bw_cv

bw_aicc <- bw.gwr(
  formula = formula_gwr,
  data = georgia_spdf,
  approach = "AICc",
  kernel = "bisquare",
  adaptive = TRUE
)
bw_aicc

# ----- fit model with the best bandwidth (regression points = data points) -----
gwr_basic <- gwr.basic(
  formula = formula_gwr,
  data = georgia_spdf,
  regression.points = georgia_spdf,
  bw = bw_aicc,
  kernel = "bisquare",
  adaptive = TRUE
)
names(gwr_basic$SDF@data)
head(gwr_basic$SDF@data)

# Extract local coefficient estimates
basic_estimates <- data.frame(
  grid_id = 1:nrow(gwr_basic$SDF@data),
  X = coordinates(gwr_basic$SDF)[, 1],
  Y = coordinates(gwr_basic$SDF)[, 2],
  gwr_basic$SDF@data
)

head(basic_estimates)
nrow(basic_estimates)

# ----- fit model with the best bandwidth (40*40 grid) -----
# x_seq <- seq(
#   from = min(georgia_spdf$X),
#   to = max(georgia_spdf$X),
#   length.out = 40
# )
# 
# y_seq <- seq(
#   from = min(georgia_spdf$Y),
#   to = max(georgia_spdf$Y),
#   length.out = 40
# )
# 
# grid_df <- expand.grid(
#   X = x_seq,
#   Y = y_seq
# )
# 
# grid_df$grid_id <- 1:nrow(grid_df)
# 
# grid_coords <- cbind(grid_df$X, grid_df$Y)
# 
# grid.spdf <- SpatialPointsDataFrame(
#   coords = grid_coords,
#   data = grid_df
# )
# 
# nrow(grid_df)
# 
# gwr_basic_grid <- gwr.basic(
#   formula = formula_gwr,
#   data = georgia_spdf,
#   regression.points = grid.spdf,
#   bw = bw_aicc,
#   kernel = "bisquare",
#   adaptive = TRUE
# )
# names(gwr_basic_grid$SDF@data)
# head(gwr_basic_grid$SDF@data)
# 
# # Extract local coefficient estimates
# basic_estimates_grid <- data.frame(
#   grid_id = 1:nrow(gwr_basic_grid$SDF@data),
#   X = coordinates(gwr_basic_grid$SDF)[, 1],
#   Y = coordinates(gwr_basic_grid$SDF)[, 2],
#   gwr_basic_grid$SDF@data
# )
# 
# head(basic_estimates_grid)
# nrow(basic_estimates_grid)

# ----- Monte Carlo permutation test -----
set.seed(123)
mc <- gwr.montecarlo(formula=formula_gwr, 
                     data=georgia_spdf, 
                     nsims=1000, 
                     kernel="bisquare",
                     adaptive=TRUE,
                     bw=bw_aicc)

# ----- Leung f123 test -----
ftest <- gwr.basic(formula=formula_gwr,
                   data=georgia_spdf,
                   bw=bw_aicc,
                   kernel="bisquare",
                   adaptive=TRUE,
                   F123.test=TRUE)
