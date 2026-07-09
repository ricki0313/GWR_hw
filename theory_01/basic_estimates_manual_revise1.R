gwr_grid <- function(X, y, obs, grid, bw){
  m <- nrow(grid)
  n <- nrow(obs)
  k <- ncol(X)
  
  # distance matrix, coefficient matrix
  D <- matrix(NA, nrow = m, ncol = n)
  for (i in 1:m){
    dx <- grid[i, 1] - obs[, 1]
    dy <- grid[i, 2] - obs[, 2]
    
    D[i, ] <- sqrt(dx^2 + dy^2)
  }
  beta_mat <- matrix(NA, nrow = m, ncol = k)
  colnames(beta_mat) <- colnames(X)
  
  # loop over each regression point, calculate coefficient matrix
  for (i in 1:m){
    d <- D[i, ]
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
    W <- diag(as.vector(w))
    beta_hat <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
    beta_mat[i, ] <- as.vector(beta_hat)
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