# ----- function -----
global_parameter_config <- function(grid_x=10L,
                                    grid_y=10L,
                                    n_predictors=2L,
                                    raw_dim=10L,
                                    sims=200L,
                                    error_var=c(0.001, 1),
                                    confidence_level=0.95,
                                    seed=123){
  # derived parameters
  n <- grid_x * grid_y
  error_sd <- sqrt(error_var)
  alpha <- 1 - confidence_level
  critical_value <- qnorm(1 - alpha/2)
  
  # GWR set
  include_intercept <- FALSE
  kernel <- "exponential"
  bw_type <- "fixed"
  bw_method <- "CV"
  
  parameter <- list(
    grid_x = grid_x,
    grid_y = grid_y,
    n = n,
    n_predictors = n_predictors,
    raw_dim = raw_dim,
    sims = sims,
    error_var = error_var,
    error_sd = error_sd,
    confidence_level = confidence_level,
    alpha = alpha,
    critical_value = critical_value,
    seed = seed,
    include_intercept = include_intercept,
    kernel = kernel,
    bw_type = bw_type,
    bw_method = bw_method
  )
  
  return(parameter)
}

grid_setup <- function(parameter){
  # x, y values
  x_values <- 1:parameter$grid_x
  y_values <- 1:parameter$grid_y
  
  # set the grid & extract all x,y
  coordinates <- expand.grid(x=x_values, y=y_values)
  x_coord <- coordinates$x
  y_coord <- coordinates$y
  
  # calculate distance matrix
  distance_matrix <- as.matrix(
    dist(coordinates[, c("x", "y")])
  )
  
  spatial_grid <- list(
    coordinates = coordinates,
    x_coord = x_coord,
    y_coord = y_coord,
    distance_matrix = distance_matrix
  )
  
  return(spatial_grid)
}

true_coef_setup <- function(grid){
  beta1_true <- grid$x_coord
  beta2_true <- grid$y_coord
  beta_true <- cbind(beta1_true=beta1_true, beta2_true=beta2_true)
  
  true_coef <- list(
    beta1_true=beta1_true,
    beta2_true=beta2_true,
    beta_true=beta_true
  )
  
  return(true_coef)
}

run_sim <- function(parameter, grid, true_coef){
  set.seed(parameter$seed)
  
  # save result of all scenarios
  all_results <- vector(mode="list", length=length(parameter$error_var))
  
  # Outer loops (different error variance scenario)
  for(i in seq_along(parameter$error_var)){
    # current var. & std.
    current_var <- parameter$error_var[i]
    current_sd <- parameter$error_sd[i]
    
    cat("Scenario", i, "error var = ", current_var, "\n")
    
    # initialize realization results (one realization for each row)
    realizations <- data.frame(
      realization = seq_len(parameter$sims),
      error_var = rep(current_var, parameter$sims),
      error_sd = rep(current_sd, parameter$sims),
      x_corr = rep(NA_real_, parameter$sims),
      optim_bw = rep(NA_real_, parameter$sims),
      rmse = rep(NA_real_, parameter$sims),
      rmse_beta1 = rep(NA_real_, parameter$sims),
      rmse_beta2 = rep(NA_real_, parameter$sims)
    )
    
    # initialize estimated coefficients for each location
    # row:`n` locations, col: `sims` realization
    beta1_est <- matrix(NA_real_, nrow=parameter$n, ncol=parameter$sims)
    beta2_est <- matrix(NA_real_, nrow=parameter$n, ncol=parameter$sims)
    beta1_coverage <- matrix(NA, nrow=parameter$n, ncol=parameter$sims)
    beta2_coverage <- matrix(NA, nrow=parameter$n, ncol=parameter$sims)
    
    # Inner loops (each simulation)
    for (r in seq_len(parameter$sims)){
      cat("\r", "scenario", i, "-realization", r, "of", parameter$sims)
      
      # 1. generate 10-dimensional N(0,I) & 100*10 matrix
      raw_data <- matrix(
        rnorm(parameter$n * parameter$raw_dim, mean=0, sd=1),
        nrow = parameter$n, 
        ncol = parameter$raw_dim
      )

      # 2. PCA (100*10 -> 100*2)
      pca_result <- prcomp(raw_data, center=TRUE, scale.=FALSE)
      X <- pca_result$x[, 1:parameter$n_predictors]
      colnames(X) <- c("x1", "x2")
      
      realizations$x_corr[r] <- cor(X[, 1], X[, 2])
      
      # 3. generate errors
      epsilon <- rnorm(n=parameter$n, mean=0, sd=current_sd)
      
      # 4. generate true y
      y <- true_coef$beta1_true*X[, 1] + true_coef$beta2_true*X[, 2] + epsilon
        
      # 5. the best bandwidth
      gwr_data <- data.frame(
        y = y,
        x1 = X[, 1],
        x2 = X[, 2],
        coord_x = grid$x_coord,
        coord_y = grid$y_coord
      )
      
      gwr_data <- sf::st_as_sf(gwr_data, coords=c("coord_x", "coord_y"), remove=FALSE)
      
      invisible(
        capture.output(
          optim_bw <- GWmodel::bw.gwr(
            formula = y ~ 0 + x1 + x2,
            data = gwr_data,
            approach = parameter$bw_method,
            kernel = parameter$kernel,
            adaptive = FALSE,
            dMat = grid$distance_matrix
          )
        )
      )
      realizations$optim_bw[r] <- optim_bw
      
      # 6. fit GWR & store local coefficient estimates
      gwr_fit <- GWmodel::gwr.basic(
        formula = y ~ 0 + x1 + x2,
        data = gwr_data,
        bw = optim_bw,
        kernel = parameter$kernel,
        adaptive = FALSE,
        dMat = grid$distance_matrix
      )
      
      gwr_result <- as.data.frame(gwr_fit$SDF)
      beta1_hat <- gwr_result$x1
      beta2_hat <- gwr_result$x2
      
      beta1_est[, r] <- beta1_hat
      beta2_est[, r] <- beta2_hat
      
      # 7. calculate standard errors & CI
      se_beta1 <- gwr_result$x1_SE
      se_beta2 <- gwr_result$x2_SE
      
      beta1_lower <- beta1_hat - parameter$critical_value*se_beta1
      beta1_upper <- beta1_hat + parameter$critical_value*se_beta1
      beta2_lower <- beta2_hat - parameter$critical_value*se_beta2
      beta2_upper <- beta2_hat + parameter$critical_value*se_beta2
      
      # 8. calculate coverage & RMSE
      # coverage
      coverage_beta1 <- (
        beta1_lower <= true_coef$beta1_true & true_coef$beta1_true <= beta1_upper
      )
      coverage_beta2 <- (
        beta2_lower <= true_coef$beta2_true & true_coef$beta2_true <= beta2_upper
      )
      
      beta1_coverage[, r] <- coverage_beta1
      beta2_coverage[, r] <- coverage_beta2
      
      # local RMSE
      rmse_beta1_r <- sqrt(mean((beta1_hat - true_coef$beta1_true)^2))
      rmse_beta2_r <- sqrt(mean((beta2_hat - true_coef$beta2_true)^2))
      
      # overall RMSE
      rmse_r <- sqrt(mean(c(
        (beta1_hat - true_coef$beta1_true)^2,
        (beta2_hat - true_coef$beta2_true)^2
      )))
      
      # store realization RMSE results
      realizations$rmse[r] <- rmse_r
      realizations$rmse_beta1[r] <- rmse_beta1_r
      realizations$rmse_beta2[r] <- rmse_beta2_r
    }
    cat("\n")
    
    # summarize current scenario
    scenario_summary <- list(
      error_var = current_var,
      error_sd = current_sd,
      mean_rmse = mean(realizations$rmse, na.rm=TRUE),
      mean_rmse_beta1 = mean(realizations$rmse_beta1, na.rm=TRUE),
      mean_rmse_beta2 = mean(realizations$rmse_beta2, na.rm=TRUE),
      mean_cp_beta1 = mean(beta1_coverage, na.rm=TRUE),
      mean_cp_beta2 = mean(beta2_coverage, na.rm=TRUE),
      mean_cp_beta1_loc = rowMeans(beta1_coverage, na.rm=TRUE),
      mean_cp_beta2_loc = rowMeans(beta2_coverage, na.rm=TRUE),
      mean_bw = mean(realizations$optim_bw, na.rm=TRUE),
      median_bw = median(realizations$optim_bw, na.rm=TRUE),
      sd_bw = sd(realizations$optim_bw, na.rm=TRUE),
      min_bw = min(realizations$optim_bw, na.rm=TRUE),
      max_bw = max(realizations$optim_bw, na.rm=TRUE)
    )
    
    # summarize all results of current scenario
    all_results[[i]] <- list(
      scenario = i,
      error_var = current_var,
      error_sd = current_sd,
      realizations = realizations,
      beta1_est = beta1_est,
      beta2_est = beta2_est,
      beta1_coverage = beta1_coverage,
      beta2_coverage = beta2_coverage,
      summary = scenario_summary
    )
  }
  names(all_results) <- paste0("error_var_", parameter$error_var)
  
  return(all_results)
}

# ----- main ----
parameter_test <- global_parameter_config(sims=200L, seed=123)
grid_test <- grid_setup(parameter_test)
true_coef_test <- true_coef_setup(grid_test)
result_test <- run_sim(parameter=parameter_test, 
                       grid=grid_test, 
                       true_coef=true_coef_test)

summary_table <- data.frame(
  error_var = c(
    result_test$error_var_0.001$summary$error_var,
    result_test$error_var_1$summary$error_var
  ),
  mean_rmse = c(
    result_test$error_var_0.001$summary$mean_rmse,
    result_test$error_var_1$summary$mean_rmse
  ),
  mean_cp_beta1 = c(
    result_test$error_var_0.001$summary$mean_cp_beta1,
    result_test$error_var_1$summary$mean_cp_beta1
  ),
  mean_cp_beta2 = c(
    result_test$error_var_0.001$summary$mean_cp_beta2,
    result_test$error_var_1$summary$mean_cp_beta2
  ),
  mean_bw = c(
    result_test$error_var_0.001$summary$mean_bw,
    result_test$error_var_1$summary$mean_bw
  )
)

summary_table
result_test$error_var_0.001$summary
result_test$error_var_1$summary