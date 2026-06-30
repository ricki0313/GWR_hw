library(GWmodel)

data(Georgia)
names(Gedu.df)
head(Gedu.df)

# ----- function -----
grid <- function(x, nrow, y, ncol){
  x_seq <- seq(
    from = min(x),
    to = max(x),
    length.out = nrow
  )
  
  y_seq <- seq(
    from = min(y),
    to = max(y),
    length.out = ncol
  )
  
  df <- expand.grid(
    X = x_seq,
    Y = y_seq
  )
  
  return(df)
}

calc_distance_matrix <- function(grid, obs){
  m <- nrow(grid)
  n <- nrow(obs)
  D <- matrix(NA, nrow = m, ncol = n)
  for (i in 1:m){
    dx <- grid[i, 1] - obs[, 1]
    dy <- grid[i, 2] - obs[, 2]
    
    D[i, ] <- sqrt(dx^2 + dy^2)
  }
  
  return(D)
}

bisquare_adaptive_weight <- function(d, bw){
  b <- sort(d)[bw]
  
  w <- numeric(length(d))
  for(j in 1:length(d)){
    if(d[j] < b){
      w[j] <- (1 - (d[j] / b)^2)^2
    }
    else{
      w[j] <- 0
    }
  }
  
  return(w)
}

local_lwr <- function(X, y, w){
  W <- diag(as.vector(w))
  beta_hat <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  
  return (beta_hat)
}

gwr_grid <- function(X, y, obs, grid, bw){
  m <- nrow(grid)
  k <- ncol(X)
  
  # distance matrix, coefficient matrix
  D <- calc_distance_matrix(grid, obs)
  beta_mat <- matrix(NA, nrow = m, ncol = k)
  colnames(beta_mat) <- colnames(X)
  
  # loop over each regression point
  for (i in 1:m){
    w_i <- bisquare_adaptive_weight(D[i,], bw)
    beta_i <- local_lwr(X, y, w_i)
    
    beta_mat[i, ] <- as.vector(beta_i)
  }
  
  # combine grid & coefficient
  result <- data.frame(
    id = 1:m,
    x_coord = grid[, 1],
    y_coord = grid[, 2],
    beta_mat
  )
  
  return(result)
}

# ----- main -----
# prepare variable: n, x, y, obs_xy, grid_xy, bw, distance matrix
n <- nrow(Gedu.df)
y <- Gedu.df$PctBach
X <- cbind(
    Intercept = 1,
    PctRural = Gedu.df$PctRural,
    PctEld = Gedu.df$PctEld,
    PctFB = Gedu.df$PctFB,
    PctPov = Gedu.df$PctPov,
    PctBlack = Gedu.df$PctBlack
)
obs_xy <- cbind(Gedu.df$X, Gedu.df$Y)
grid_df <- grid(Gedu.df$X, 40, Gedu.df$Y, 40)
grid_xy <- matrix(
  c(grid_df$X, grid_df$Y),
  ncol = 2
)
bw <- 136

# run
gwr_result <- gwr_grid(X=X, y=y, obs=obs_xy, grid=grid_xy, bw=bw)
names(gwr_result)
head(gwr_result)

# ----- compare with GWmodel -----
coef_names <- c(
  "Intercept",
  "PctRural",
  "PctEld",
  "PctFB",
  "PctPov",
  "PctBlack"
)

diff_result <- gwr_result[, coef_names] - grid_basic_estimates[, coef_names]
head(diff_result)
summary(diff_result)
max(abs(as.matrix(diff_result)))
