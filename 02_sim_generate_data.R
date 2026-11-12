# Project: MFPCA digital outcomes 
# Script Purpose: functions to generate ECG curves for simulation study   
# Author: Mia Tackney 
# Date: 2026-04-15
# GitHub: https://github.com/mst1g15/MFPCA_digital_outcomes


generate_healthy_group <- function(res, all_dat_name, increments){
  #Purpose: Given an mfpca object from a reference dataset, use the estimated
  #components to generate a dataset which is similar to the reference dataset.
  #input: res - refund mfpca object 
  #all_dat_name: subject identified for each curves 
  #increments: time increments
  #output: a list containing (1) identifiers for each curve, and two datasets 
  #representing ECG curves from two measurements occasions. 
  
  mu  <- res$mu                        # Mean function
  T_len    <- length(mu)                  # Grid length (e.g., 308)
  phi_b <- res$efunctions$level1       # Between eigenfunctions (T x Kb)   #308 x 3
  phi_w <- res$efunctions$level2       # Within eigenfunctions (T x Kw)   #308 x 3
  scores_b <- res$scores$level1        # Between scores (N_subjects x Kb)   #64 * 3
  scores_w <- res$scores$level2        # Within scores (N_obs x Kw)   #1392 x 3
  
  # Estimated measurement error SD from refund object
  sigma_eps <- sqrt(max(res$sigma2, 0))
  
  all_ids <- all_dat_name
  id_unique <- unique(all_dat_name)    # Which subject each observation belongs to
  time <- increments                   # Evaluation grid
  score_b_mean <- colMeans(scores_b)
  scores_b_sigma <- cov(scores_b)
  score_w_mean <- colMeans(scores_w)
  scores_w_sigma <- cov(scores_w)
  ids <- c()
  
  #dataset for two measurement occasions 
  X_hat_day1 <- X_hat_day2 <- c()
  
  for (ii in 1:length(id_unique)) {
    subj <- id_unique[ii]
    ind <- which(all_ids == subj)
    n_i <- length(ind)
    
    # Non-parametric bootstrap sampling for between-subject scores:
    sampled_row <- sample(1:nrow(scores_b), size = 1, replace = TRUE)
    sb <- as.numeric(scores_b[sampled_row, ])
    
    #within subject scores obtained via multivarite normal distributino 
    sw1 <- mvrnorm(n = n_i, mu = score_w_mean, Sigma = scores_w_sigma) 
    sw2 <- mvrnorm(n = n_i, mu = score_w_mean, Sigma = scores_w_sigma) 
    
    # Independent measurement error matrices (T_len x n_i)
    err1 <- matrix(rnorm(T_len * n_i, mean = 0, sd = sigma_eps), nrow = T_len, ncol = n_i)
    err2 <- matrix(rnorm(T_len * n_i, mean = 0, sd = sigma_eps), nrow = T_len, ncol = n_i)
    
    # make matrices compatible for addition
    mu_mat <- matrix(mu, nrow = T_len, ncol = n_i)
    
    between_mat <- matrix(phi_b %*% sb, nrow = T_len, ncol = n_i)
    
    # reconstruct
    X_subj_day1 <- mu_mat + between_mat + phi_w %*% t(sw1) +err1
    X_subj_day2 <- mu_mat + between_mat + phi_w %*% t(sw2) +err2
    
    X_hat_day1 <- rbind(X_hat_day1, t(X_subj_day1)) 
    X_hat_day2 <- rbind(X_hat_day2, t(X_subj_day2)) 
    
    ids <- c(ids, rep(subj, n_i))
    
  }
  return(list(ids, X_hat_day1, X_hat_day2))
  
  
}



generate_changes <- function(X_hat, change, ids, ...){
  #purpose: given ECG curves from a healthy group, induce a specific change 
  #input: X_hat - a matrix containing curves from the healthy group 
  #       change - a character for the change, either "ST", "changing_T", "flat_T_high", reduce_all_high"
  #       ids - an identifier for the subject
  #output: a matrix containing updated curves (with changes induced)
  
  #ST elevation-----------------------------------------------------------------
  
  if(change=="ST"){
    
    Y_mat_new <- X_hat
    
    A <- 120          # height of elevation in mV
    t_J <- 0.32      # J-point (s)
    t_STend <- 0.42   # end of ST segment (s)
    mu <- (t_J + t_STend)/2
    sigma <- (t_STend - t_J)/2  # controls width of bump
    
    for (i in 1:nrow(Y_mat_new)){
      pert <- A*exp(-0.5*((increments-mu)/sigma)^2)
      Y_mat_new[i, ] <-   X_hat[i, ] + pert
    }
    
  #changing T-wave amplitude----------------------------------------------------
    
  }else if(change=="changing_T"){
    
    Y_mat_new <- X_hat
    
    alpha_T <- ids %>%
      tibble(val = .) %>%
      group_by(val) %>%
      mutate(sampled = sample(c(0.7, 0.2), 1)) %>%  #changing amplitudes for T 
      pull(sampled)
    
    t_T <- 0.508; sigma_T <- 0.1; 
    
    for (i in 1:nrow(Y_mat_new)){
      
      T_reduction <- alpha_T[i]* exp(-0.5 * ((increments - t_T) / sigma_T)^2)
      
      Y_mat_new[i, ] <-   X_hat[i, ] * (1- T_reduction)
      
    }
    
  #flattening of T wave to a high degree-----------------------------------------
    
  }else if(change=="flat_T_high"){
    
    Y_mat_new <- X_hat
    
    #reduction in T
    t_T <- 0.508; sigma_T <- 0.1; alpha_T <- 0.5  #strong amount of suppression 
    
    for (i in 1:nrow(Y_mat_new)){
      
      T_reduction <- alpha_T* exp(-0.5 * ((increments - t_T) / sigma_T)^2)
      
      Y_mat_new[i, ] <-   X_hat[i, ] * (1-T_reduction)
      
    }
    
  #reduce all amplitudes to a high degree---------------------------------------
  }else if(change=="reduce_all_high"){
    
    Y_mat_new <- X_hat
    
    t_T <- 0.508; sigma_T <- 0.1; alpha_T <- 0.5 
    t_R <- 0.25; sigma_R <- 0.1; alpha_R <- 0.5  
    t_P <- 0.104; sigma_P<- 0.1; alpha_P <- 0.5 
    
    for (i in 1:nrow(Y_mat_new)){
      
      T_reduction <- alpha_T* exp(-0.5 * ((increments - t_T) / sigma_T)^2)
      R_reduction <- alpha_R* exp(-0.5 * ((increments - t_R) / sigma_R)^2)
      P_reduction <- alpha_P* exp(-0.5 * ((increments - t_P) / sigma_P)^2)
      
      Y_mat_new[i, ] <-   X_hat[i, ] * (1-R_reduction) * (1-T_reduction)  * (1-P_reduction) 
      
    }
  }
  
  return(Y_mat_new)
  
}
