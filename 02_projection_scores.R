
# Project: MFPCA digital outcomes 
# Script Purpose: functions to obtain MFPCA projection scores  
# Author: Mia Tackney 
# Date: 2026-04-15
# GitHub: https://github.com/mst1g15/MFPCA_digital_outcomes



predict_fpca_same_grid <- function(fpca_obj, newdata) {
  #Purpose: obtain fpca projection scores for a new dataset, given an fpca object 
  #Input: fpca_obj - is an fpca object from the refund package; 
  #       this has been fit on a reference dataset
  #       Note this is in the univariate fpca setting.
  #       newdata - is a dataset from a new group
  #Output: new_scores - a vector of scores for each individual
  #       in the newdata dataset. 
  
  
  #ensure newdata is a matrix
  newdata <- as.matrix(newdata)
  
  # extract components from the fpca object
  mu  <- fpca_obj$mu           # The estimated mean function
  phi <- fpca_obj$efunctions  # The eigenfunctions (basis)
  
  #center the new data: subtract the mean from every row
  newdata_centered <- sweep(newdata, 2, mu, "-")
  
  #calculate Scores: Projection onto the eigenfunctions
  # Scores = (Y - mu) %*% phi
  new_scores <- newdata_centered %*% phi
  
  
  return(new_scores)
}



predict_mfpca_scores_wide <- function(object,
                                      Y_new_matrix,
                                      id_new_vector,
                                      sigma_e = object$sigma2) {
  # Purpose: Obtain MFPCA projection scores for a new dataset given a fitted refund MFPCA object
  # Inputs:
  #   object        - Two-level MFPCA object from the refund package
  #   Y_new_matrix  - Matrix of new curves (N_curves x L)
  #   id_new_vector - Subject identifier vector of length N_curves
  #   sigma_e       - Estimated pointwise error variance (defaults to object$sigma2)
  
  if (is.null(sigma_e) || missing(sigma_e))
    stop("Must supply measurement error variance sigma_e (e.g., object$sigma2).")
  
  Y_new_matrix <- as.matrix(Y_new_matrix)
  N <- nrow(Y_new_matrix)
  L <- ncol(Y_new_matrix)
  if (length(id_new_vector) != N)
    stop("id_new_vector length must match number of curves (N).")
  
  # Extract components from reference MFPCA object
  mu_hat <- as.numeric(object$mu)
  if (length(mu_hat) != L)
    stop("Length of object$mu must match # of columns in Y_new_matrix (L).")
  
  Phi1 <- as.matrix(object$efunctions$level1)    # L x K1
  Phi2 <- as.matrix(object$efunctions$level2)    # L x K2
  
  K1 <- ncol(Phi1)
  K2 <- ncol(Phi2)
  
  lambda1 <- as.numeric(object$evalues$level1)   # length K1
  lambda2 <- as.numeric(object$evalues$level2)   # length K2
  
  if (length(lambda1) != K1)
    stop("object$evalues$level1 length must equal ncol(Phi1).")
  if (length(lambda2) != K2)
    stop("object$evalues$level2 length must equal ncol(Phi2).")
  
  # Uniform quadrature weights on [0, 1]
  weights <- rep(1 / L, L)
  
  # Center new curves using reference mean function
  Y_centered <- sweep(Y_new_matrix, 2, mu_hat, "-")
  
  # Raw Level-1 projections (c_{ijk1})
  raw_L1 <- Y_centered %*% (Phi1 * weights)      # N x K1
  
  # Subject-average Level-1 projections (\bar{c}_{ik1})
  df <- data.frame(ID = id_new_vector, raw_L1)
  mean_by_id <- aggregate(df[, -1, drop = FALSE],
                          by = list(ID = df$ID),
                          FUN = mean)
  
  subject_ids <- mean_by_id$ID
  cbar <- as.matrix(mean_by_id[, -1, drop = FALSE])   # n_subjects x K1
  
  # Number of curves per subject (n_i)
  n_tbl <- table(id_new_vector)
  n_vec <- as.numeric(n_tbl[match(subject_ids, names(n_tbl))])
  
  # EBLUP Level-1 score estimation (\hat{\eta}_{ik1})
  n_subjects <- length(subject_ids)
  b_hat <- matrix(0, n_subjects, K1)
  
  for (i in seq_len(n_subjects)) {
    shrink <- lambda1 / (lambda1 + sigma_e / (L * n_vec[i]))
    b_hat[i, ] <- shrink * cbar[i, ]
  }
  rownames(b_hat) <- subject_ids
  
  # Expand back to N curves according to original ID order
  b_hat_expanded <- b_hat[match(id_new_vector, subject_ids), , drop = FALSE]
  
  # Level-1 fitted functions (\hat{f}^{(1)}_i) and subject-detrended residuals (r_{ij})
  F1 <- b_hat_expanded %*% t(Phi1)               # N x L
  residuals <- Y_centered - F1
  
  # Raw Level-2 projections (\tilde{c}_{ijk2})
  raw_L2_tilde <- residuals %*% (Phi2 * weights) # N x K2
  
  # EBLUP Level-2 score estimation (\hat{\xi}_{ijk2})
  L2_scores <- matrix(0, N, K2)
  for (k in seq_len(K2)) {
    shrink <- lambda2[k] / (lambda2[k] + sigma_e / L)
    L2_scores[, k] <- shrink * raw_L2_tilde[, k]
  }
  
  final_scores <- cbind(b_hat_expanded, L2_scores)
  colnames(final_scores) <- c(
    paste0("L1_PC", 1:K1),
    paste0("L2_PC", 1:K2)
  )
  
  return(data.frame(ID = id_new_vector, final_scores))
}
