# Project: MFPCA digital outcomes 
# Script Purpose: additional helper functions for running simulation study  
# Author: Mia Tackney 
# Date: 2026-04-15
# GitHub: https://github.com/mst1g15/MFPCA_digital_outcomes






compute_ICC <- function(dataset, outcome){
  #Purpose: Compute ICC(C,1) (consistency) and ICC(A,1) (agreement), given summary scores  
  #Input: dataset - a dataframe assumed to contain an outcome column, an ids column and days column
  #       outcome - a character for the column name of the ids column 
  #Output: the estimated ICC(C,1) and ICC(A,1) values 
  
  lmm.1 = lmer(paste0(outcome,"~ 1 + (1|ids)"), data=dataset, REML=TRUE) 
  lmm.1.Person = as.data.frame(VarCorr(lmm.1))$vcov[1]
  lmm.1.res = as.data.frame(VarCorr(lmm.1))$vcov[2]
  ICC.consistency = lmm.1.Person/(lmm.1.Person+lmm.1.res)
  
  lmm.1 = lmer(paste0(outcome,"~ 1 + (1|ids) + (1|day)"), data=dataset, REML=TRUE) 
  lmm.1.Person = as.data.frame(VarCorr(lmm.1))$vcov[1]
  lmm.1.Days = as.data.frame(VarCorr(lmm.1))$vcov[2]
  lmm.1.res = as.data.frame(VarCorr(lmm.1))$vcov[3]
  ICC.agreement = lmm.1.Person/(lmm.1.Person+lmm.1.Days+lmm.1.res)
    
    
  return(c(ICC.consistency, ICC.agreement))

}

