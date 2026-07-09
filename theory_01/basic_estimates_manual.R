library(GWmodel)

data(Georgia)
names(Gedu.df)
head(Gedu.df)

# ============================================================
# Basic estimate 
# ============================================================
# ---------- helper functions ----------
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

# ---------- main GWR estimation functions ----------
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
  
  # standard deviation of each coefficient
  coef_sd <- numeric(k)
  for (i in 1:k){
    beta_bar <- mean(beta_mat[, i])
    
    ss <- 0
    for (j in 1:m){
      ss <- ss + (beta_mat[j, i] - beta_bar)^2
    }
    coef_sd[i] <- sqrt(ss / (m - 1))
  }
  names(coef_sd) <- colnames(X)
  
  # combine grid & coefficient
  result <- data.frame(
    id = 1:m,
    x_coord = grid[, 1],
    y_coord = grid[, 2],
    beta_mat
  )
  
  return(list(
    coef_estimates = result,
    coef_sd = coef_sd
  ))
}

# ============================================================
# bandwidth selection
# ============================================================
# ---------- aicc ----------
calc_gwr_aicc <- function(X, y, obs, bw){
  D <- calc_distance_matrix(obs, obs)
  
  n <- nrow(X)
  y_hat <- numeric(n)
  tr_S <- 0
  for (i in 1:n){
    # y_hat
    w_i <- bisquare_adaptive_weight(D[i, ], bw)
    beta_i <- local_lwr(X, y, w_i)
    y_hat[i] <- X[i, ,drop=FALSE] %*% beta_i
    
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

select_bw_aicc <- function(X, y, obs, bw_candidates){
  aicc_table <- data.frame(
    bw = bw_candidates,
    aicc = NA
  )
  
  for (k in 1:length(bw_candidates)){
    bw_k <- bw_candidates[k]
    result_k <- calc_gwr_aicc(X=X, y=y, obs=obs, bw=bw_k)
    
    aicc_table$aicc[k] <- result_k
  }
  
  min_aicc <- aicc_table[which.min(aicc_table$aicc), ]
  
  return(list(
    all_results = aicc_table,
    best = min_aicc
  ))
}

# ---------- cv ----------
calc_gwr_cv <- function(X, y, obs, bw){
  D <- calc_distance_matrix(obs, obs)
  n <- nrow(X)
  y_hat <- numeric(n)
  for(i in 1:n){
    w_i <- bisquare_adaptive_weight(D[i, ], bw)
    
    # leave one out
    w_i[i] <- 0
    
    # fit
    beta_i <- local_lwr(X, y, w_i)
    y_hat[i] <- X[i, ,drop=FALSE] %*% beta_i
  }
  cv <- sum((as.vector(y) - y_hat)^2)
  
  return(cv)
}

select_bw_cv <- function(X, y, obs, bw_candidates){
  cv_table <- data.frame(
    bw = bw_candidates,
    cv = NA
  )
  for (k in seq_along(bw_candidates)){
    bw_k <- bw_candidates[k]
    cv_k <- calc_gwr_cv(X=X, y=y, obs=obs, bw=bw_k)
    cv_table$cv[k] <- cv_k
  }
  
  min_cv <- cv_table[which.min(cv_table$cv), ]
  
  return(list(
    all_results = cv_table,
    best = min_cv
  ))
}

# ---------- gcv ----------
calc_gwr_gcv <- function(X, y, obs, bw){
  D <- calc_distance_matrix(obs, obs)
  
  n <- nrow(X)
  y_hat <- numeric(n)
  tr_S <- 0
  for (i in 1:n){
    # y_hat
    w_i <- bisquare_adaptive_weight(D[i, ], bw)
    beta_i <- local_lwr(X, y, w_i)
    y_hat[i] <- X[i, ,drop=FALSE] %*% beta_i
    
    # trace(S)
    W_i <- diag(as.vector(w_i))
    s_i <- X[i, ,drop=FALSE] %*% solve(t(X) %*% W_i %*% X) %*% t(X) %*% W_i # row i of S
    tr_S <- tr_S + s_i[1, i]
  }
  gcv <- n * (as.vector(y) - y_hat)^2 / (n - tr_S)^2
  
  return(gcv)
}

# ============================================================
# Monte Carlo permutation test 
# ============================================================
mc_per <- function(X, y, obs, bw, sims, seed=123){
  set.seed(seed)
  sd_ori <- gwr_grid(X=X, y=y, obs=obs, grid=obs, bw)$coef_sd
  k <- length(sd_ori)
  
  # whether simulation std. is larger than original std.
  larger_than_ori <- matrix(NA, nrow=sims, ncol=k)
  colnames(larger_than_ori) <- names(sd_ori) 
  
  # std. of random location beta
  sd_rand_mat <- matrix(NA, nrow=sims, ncol=k)
  colnames(sd_rand_mat) <- names(sd_ori)
  
  for (i in 1:sims){
    # permutation
    obs_rand <- obs[sample(1:nrow(obs)), ]
    
    # std. of randomized beta
    sd_rand <- gwr_grid(X=X, y=y, obs=obs_rand, grid=obs_rand, bw)$coef_sd
    sd_rand_mat[i, ] <- sd_rand
  }
  
  # calculate p-value
  pv <- colMeans(sd_rand_mat >= matrix(sd_ori, nrow=sims, ncol=k, byrow=TRUE))
  
  return(pv)
}

# ============================================================
# main 
# ============================================================
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

bw_candidates <- seq(130, 145, by=1)
bw_aicc <- select_bw_aicc(X=X, y=y, obs=obs_xy, bw_candidates=bw_candidates)$best$bw
bw_cv <- select_bw_cv(X=X, y=y, obs=obs_xy, bw_candidates=bw_candidates)$best$bw


# run
gwr_result <- gwr_grid(X=X, y=y, obs=obs_xy, grid=grid_xy, bw=bw)
names(gwr_result$coef_estimates)
head(gwr_result$coef_estimates)

p_value_mc <- mc_per(X=X, y=y, obs=obs_xy, bw=bw, sims=1000)
p_value_mc

# ========== compare with GWmodel (must run package version first) ==========
coef_names <- c(
  "Intercept",
  "PctRural",
  "PctEld",
  "PctFB",
  "PctPov",
  "PctBlack"
)

diff_result <- gwr_result$coef_estimates[, coef_names] - grid_basic_estimates[, coef_names]
head(diff_result)
summary(diff_result)
max(abs(as.matrix(diff_result)))