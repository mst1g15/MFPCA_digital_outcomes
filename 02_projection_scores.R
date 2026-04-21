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
                                      sigma_e) {
  #Purpose: obtain mfpca projection scores for a new dataset, given an mfpca object 
  #Input: object - is an mfpca object from the refund package; 
  #       this has been fit on a reference dataset
  #       Note this is in the two-level fpca setting.
  #       Y_new_matrix - is a dataset from a new group
  #       id_new_vector - this is an identifier for curves from 
  #       the same individual (since we are in the repeated measures setting)
  #       sigma_e - this is the error vector from the mfcpa
  #Output: new_scores - a matrix of scores for each individual and each curves
  #       in the new dataset.
  
  #Check input------------------------------------------------------------------
  
  if (missing(sigma_e))
    stop("Must supply measurement error variance sigma_e.")
  
  Y_new_matrix <- as.matrix(Y_new_matrix)
  N <- nrow(Y_new_matrix)
  L <- ncol(Y_new_matrix)
  if (length(id_new_vector) != N)
    stop("id_new_vector length must match number of curves (N).")
  
  #extract components from object-----------------------------------------------
  mu_hat <- as.numeric(object$mu)
  if (length(mu_hat) != L)
    stop("Length of object$mu must match # of columns in Y_new_matrix (L).")
  
  Phi1 <- as.matrix(object$efunctions$level1)    # L × K1
  Phi2 <- as.matrix(object$efunctions$level2)    # L × K2
  
  K1 <- ncol(Phi1)
  K2 <- ncol(Phi2)
  
  lambda1 <- as.numeric(object$evalues$level1)  # length K1
  lambda2 <- as.numeric(object$evalues$level2)  # <--- NEW: Level 2 eigenvalues
  
  if (length(lambda1) != K1)
    stop("object$evalues$level1 length must equal ncol(Phi1).")
  if (length(lambda2) != K2)
    stop("object$evalues$level2 length must equal ncol(Phi2).")
  
  # uniform grid weights
  weights <- rep(1 / L, L)
  
  #Center data------------------------------------------------------------------
  
  Y_centered <- sweep(Y_new_matrix, 2, mu_hat, "-")
  
  #Obtain raw Level 1 Projections-----------------------------------------------
  # Raw projection of each curve onto L1 basis (c_ij in the notation)
  raw_L1 <- Y_centered %*% (Phi1 * weights)    # N × K1
  
  #Compute subject-average Level 1 projections
  df <- data.frame(ID = id_new_vector, raw_L1)
  
  #Calculate mean projection per subject (c̄_i in the notation)
  mean_by_id <- aggregate(df[ , -1, drop = FALSE],
                          by = list(ID = df$ID),
                          FUN = mean)
  
  subject_ids <- mean_by_id$ID
  cbar <- as.matrix(mean_by_id[ , -1, drop = FALSE])    # n_subject × K1
  
  #Number of curves (visits) per subject (m_i)
  m_tbl <- table(id_new_vector)
  m_vec <- as.numeric(m_tbl[match(subject_ids, names(m_tbl))])
  
  #EBLUP Level-1 score estimation (between-subect Scores)-------------------------
  # Formula: b_i = (λ1 / (λ1 + σ_e/m_i)) * c̄_i
  n_subjects <- length(subject_ids)
  b_hat <- matrix(0, n_subjects, K1)
  
  for (i in seq_len(n_subjects)) {
    # Calculate shrinkage factor for each PC, based on the subject's m_i
    shrink <- lambda1 / (lambda1 + sigma_e / m_vec[i])
    b_hat[i, ] <- shrink * cbar[i, ]
  }
  
  rownames(b_hat) <- subject_ids
  
  #Expand back to N curves according to ID order
  b_hat_expanded <- b_hat[match(id_new_vector, subject_ids), , drop = FALSE]
  
  # Level-1 fitted functions------------------------------------------------
  #F1: Fitted function based on L1 scores: mu_hat + Phi1 %*% b_hat
  F1 <- b_hat_expanded %*% t(Phi1)    # N × L
  
  #EBLUP Level-2 scores (within-subject scores)---------------------------------
  # Calculate residuals after removing the L1 fitted part
  residuals <- Y_centered - F1
  
  # Raw projection of residuals onto L2 basis (̃c_ij in the notation)
  raw_L2_tilde <- residuals %*% (Phi2 * weights)  # N × K2
  
  # Apply shrinkage for Level 2 scores <--- THE CRUCIAL CORRECTION
  # Formula: a_ij = (λ2 / (λ2 + σ_e)) * ̃c_ij
  L2_scores <- matrix(0, N, K2)
  for (k in seq_len(K2)) {
    # Calculate shrinkage factor (does not depend on m_i or ID)
    shrink <- lambda2[k] / (lambda2[k] + sigma_e)
    L2_scores[, k] <- shrink * raw_L2_tilde[, k]
  }
  
  # Final combined scores
  
  
  final_scores <- cbind(b_hat_expanded, L2_scores)
  colnames(final_scores) <- c(
    paste0("L1_PC", 1:K1),
    paste0("L2_PC", 1:K2)
  )
  
  return(data.frame(ID = id_new_vector, final_scores))
}

