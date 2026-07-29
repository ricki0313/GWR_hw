golden_search_ori <- function(f, a, c , error=sqrt(.Machine$double.eps), max_iter=1000){
  if(a >= c){
    stop("a must be smaller than c")
  }
  
  rho <- (sqrt(5) - 1) / 2
  x1 <- c - (c - a)*rho
  x2 <- a + (c - a)*rho
  
  f1 <- f(x1)
  f2 <- f(x2)
  
  iter <- 0
  while((c-a) > error && iter <= max_iter){
    if(f1 < f2){
      c <- x2
      
      x2 <- x1
      f2 <- f1
      
      x1 <- c - (c - a)*rho
      f1 <- f(x1)
    }else{
      a <- x1
      
      x1 <- x2
      f1 <- f2
      
      x2 <- a + (c-a)*rho
      f2 <- f(x2)
    }
    
    iter <- iter + 1
  }
  
  xmin <- (a + c) / 2
  
  list(
    min = xmin,
    value = f(xmin),
    interval = c(a,c),
    iter = iter,
    converged = iter < max_iter
  )
}

golden_search <- function(f, a, c , adaptive=FALSE,
                          error=1e-4, 
                          max_iter=1000){
  if(a >= c){
    stop("a must be smaller than c")
  }
  
  rho <- (sqrt(5) - 1) / 2
  
  # initial golden section points
  d <- (c - a)*rho
  
  if(adaptive){
    x1 <- round(c - d)
    x2 <- round(a + d)
  }else{
    x1 <- c - d
    x2 <- a + d
  }
  
  f1 <- f(x1)
  f2 <- f(x2)
  
  # difference between 2 objective points
  d1 <- f2 - f1
  
  iter <- 0
  while(abs(d) > error &&
        abs(d1)> error &&
        iter < max_iter){
    d <- rho * d
    
    if(f1 < f2){
      c <- x2
      
      x2 <- x1
      f2 <- f1
      
      if(adaptive){x1 <- round(c - d)}else{x1 <- c - d}
      f1 <- f(x1)
      
    }else{
      a <- x1
      
      x1 <- x2
      f1 <- f2
      
      if(adaptive){x2 <- round(a + d)}else{x2 <- a + d}
      f2 <- f(x2)
    }
    
    d1 <- f2 - f1
    iter <- iter + 1
  }
  
  if(adaptive){
    bw_canditate <- seq(from=ceiling(a), to=floor(c), by=1)
    value <- sapply(bw_canditate, f)
    min_id <- which.min(value)
    
    xmin <- bw_canditate[min_id]
    fmin <- value[min_id]
  }else{
    if(f1 < f2){
      xmin <- x1
      fmin <- f1
    }else{
      xmin <- x2
      fmin <- f2
    }
  }
  
  list(
    min = xmin,
    value = fmin,
    interval = c(a,c),
    iter = iter
  )
}
