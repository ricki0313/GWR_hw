library(GWmodel)

data(Georgia)
names(Gedu.df)
head(Gedu.df)

# ========== function ==========
# ----- basic estimate -----
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

# ----- AICc comparison -----
calc_gwr_aicc <- function(X, y, obs, bw){
  D <- calc_distance_matrix(obs, obs)
  
  n <- nrow(X)
  y_hat <- numeric(n)
  tr_S <- 0
  for (i in 1:n){
    # y_hat
    w_i <- bisquare_adaptive_weight(D[i, ], bw)
    beta_i <- local_lwr(X, y, w_i)
    y_hat[i] <- X[i, ,drop=FALSE] %*%beta_i
    
    # trace(S)
    W_i <- diag(as.vector(w_i))
    s_i <- X[i, ,drop=FALSE] %*% solve(t(X) %*% W_i %*% X) %*% t(X) %*% W_i # row i of S
    tr_S <- tr_S + s_i[1, i]
  }
  
  # sigma_hat for AICc
  res <- as.vector(y) - y_hat
  RSS <- sum(res^2)
  sigma_hat <- sqrt(RSS / n)
  
  # AICc
  aicc <- 2*n*log(sigma_hat) + n*log(2 * pi) + n*((n + tr_S)/(n - 2 - tr_S))
  return(aicc)
}

select_bw_aicc <- function(X, y, obs, bw_canditates){
  aicc_table <- data.frame(
    bw = bw_canditates,
    aicc = NA
  )
  
  for (k in 1:length(bw_canditates)){
    bw_k <- bw_canditates[k]
    result_k <- calc_gwr_aicc(X=X, y=y, obs=obs, bw=bw_k)
    
    aicc_table$aicc[k] <- result_k
  }
  
  min_aicc <- aicc_table[which.min(aicc_table$aicc), ]
  
  return(list(
    all_results = aicc_table,
    best = min_aicc
  ))
}

# ========== main ==========
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

bw_canditates <- seq(130, 140, by=1)
bw <- select_bw_aicc(X=X, y=y, obs=obs_xy, bw_canditates=bw_canditates)$best$bw

# run
gwr_result <- gwr_grid(X=X, y=y, obs=obs_xy, grid=grid_xy, bw=bw)
names(gwr_result)
head(gwr_result)

# ----- compare with GWmodel (must run package version first) -----
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
