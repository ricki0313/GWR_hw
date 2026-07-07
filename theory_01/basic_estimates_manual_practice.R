  library(GWmodel)
  
  data(Georgia)
  names(Gedu.df)
  
  # ========== function ==========
  # ----- basic estimate -----
  grid <- function(x, nrow, y, ncol){
    seq_x <- seq(from=min(x), to=max(x), length.out=nrow)
    seq_y <- seq(from=min(y), to=max(y), length.out=ncol)
    df <- expand.grid(X=seq_x, Y=seq_y)
    
    return(df)
  }
  
  cal_distance_matrix <- function(grid, obs){
    m <- nrow(grid)
    n <- nrow(obs)
    D <- matrix(data=NA, nrow=m, ncol=n)
    for(i in 1:m){
      dx <- grid[i, 1] - obs[, 1]
      dy <- grid[i, 2] - obs[, 2]
      D[i, ] <- sqrt(dx^2 + dy^2) 
    }
    
    return(D)
  }
  
  bisquare_adaptive_weight <- function(d, bw){
    b <- sort(d)[bw]
    
    w <- numeric(length(d))
    for(i in 1:length(d)){
      if(d[i] < b){
        w[i] <- (1 - (d[i] / b)^2)^2      
      }else{
        w[i] <- 0
      }
    }
    
    return(w)
  }
  
  local_lwr <- function(X, y, w){
    W <- diag(as.vector(w))
    beta_hat <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
    
    return(beta_hat)
  }
  
  gwr_grid <- function(X, y, obs, grid, bw){
    m <- nrow(grid)
    k <- ncol(X)
    
    D <- cal_distance_matrix(grid=grid, obs=obs)
    
    # loop over each regression point
    beta_mat <- matrix(NA, nrow=m, ncol=k)
    colnames(beta_mat) <- colnames(X)
    for(i in 1:m){
      w_i <- bisquare_adaptive_weight(d=D[i, ], bw=bw)
      beta_i <- local_lwr(X=X, y=y, w=w_i)
      beta_mat[i, ] <- as.vector(beta_i)
    }
    
    # result
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
    D <- cal_distance_matrix(obs, obs)
    n <- nrow(X)
    
    y_hat <- numeric(n)
    tr_S <- 0
    for(i in 1:n){
      # y_hat
      w_i <- bisquare_adaptive_weight(d=D[i, ], bw=bw)
      beta_i <- local_lwr(X=X, y=y, w=w_i)
      y_hat[i] <- X[i, , drop=FALSE] %*% beta_i
      
      # trace(S)
      W_i <- diag(as.vector(w_i))
      s_i <- X[i, , drop=FALSE] %*% solve(t(X) %*% W_i %*% X) %*% t(X) %*% W_i
      tr_S <- tr_S + s_i[1, i]
    }
  
    # sigma_hat
    res <- as.vector(y_hat) - y
    RSS <- sum(res^2)
    sigma_hat <- sqrt(RSS / n)
    
    # aicc
    aicc <- 2*n*log(sigma_hat) + n*log(2 * pi) + n*((n + tr_S)/(n - 2 - tr_S))
    
    return(aicc)
  }
  
  select_bw_aicc <- function(X, y, obs, bw_candidates){
    aicc_table <- data.frame(
      bw = bw_candidates,
      aicc = NA
    )
    
    for(k in 1:length(bw_candidates)){
      bw_k <- bw_candidates[k]
      aicc_k <- calc_gwr_aicc(X=X, y=y, obs=obs, bw=bw_k)
      
      aicc_table$aicc[k] <- aicc_k
    }
    
    min_aicc <- aicc_table[which.min(aicc_table$aicc),]
    
    return(list(
      all_results = aicc_table,
      best = min_aicc
    ))
  }
  
