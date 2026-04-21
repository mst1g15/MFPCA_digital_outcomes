# Project: MFPCA digital outcomes 
# Script Purpose: functions to obtain MFPCA projection scores  
# Author: Mia Tackney 
# Date: 2026-04-15
# GitHub: https://github.com/mst1g15/MFPCA_digital_outcomes


#load packages and functions------------------------------------------------
source("00_init.R")
source("02_projection_scores.R")
source("02_sim_generate_data.R")
source("02_sim_functions.R")

#Load healthy landmarked ECGs which is the reference group
#- note this data is not made available on Github
load("Analysis_Outputs/healthy_ECG_landmarked.RData")


#mfpca with refund package----------------------------------------------------
Y_mat <- t(eval.fd(increments, regfd))
all_dat_name <- ecg_dat$Name
k=5
res <- mfpca.face(Y = Y_mat, id = all_dat_name, visit=NULL, npc = c(k,k))


#obtain average curve for each subject and run fpca----------------------------
df_curves <- as.data.frame(Y_mat)
df_curves <- df_curves %>%
  mutate(ID = all_dat_name)

# summarize by subject
mean_curves_df <- df_curves %>%
  group_by(ID) %>%
  summarise(across(.cols = everything(), .fns = ~ mean(.x)), .groups = "drop")

fpca_single <- refund::fpca.face(
  Y = as.matrix(mean_curves_df[,-1]),
  argvals = NULL,         # or your time grid
  npc = 5
)


#Simulation study-------------------------------------------------------------- 

nsim=5000
change_types <- c("none", "ST", "changing_T", "flat_T_high", "reduce_all_high")

results_all <- c()

for(change_type in change_types){
  
  for(i in 1:nsim){
    
    #generate healthy dataset--------------------------------------------------
    healthy <- generate_healthy_group(res, all_dat_name, increments)
    
    #id name 
    healthy_id <- healthy[[1]] 
    
    #healthy data from day 1 
    healthy_dat <- healthy[[2]]
    
    #healthy data from day 2
    healthy_dat_day2 <- healthy[[3]]
    
    #generate another healthy dataset ------------------------------------------
    healthy2<- generate_healthy_group(res, all_dat_name, increments)
    
    #ids
    healthy2_id <- healthy2[[1]]
    
    #healthy data from day 1
    healthy2_dat <- healthy2[[2]]
    
    #healthy data from day 2
    healthy2_dat_day2 <- healthy2[[3]]
    
    #induce changes -------------------------------------------------------------------
    if(change_type=="none"){
      dat_disease <- healthy2_dat; dat_disease_day2 <- healthy2_dat_day2
      dis_gp <- "none"
    } else {
      dat_disease <- generate_changes(healthy2_dat, change=change_type, ids=healthy_id)
      dat_disease_day2 <- generate_changes(healthy2_dat_day2, change=change_type, ids = healthy2_id)
      dis_gp <- "disease"
    }
    
    
    #Obtain amplitudes------------------------------------------------------------
    #Create one dataframe for amplitudes
    
    all_results <- rbind(
      data.frame(ids=healthy_id, group="healthy", day=1, P=healthy_dat[,54], R=healthy_dat[,129], T=healthy_dat[,261]),
      data.frame(ids=healthy_id, group="healthy", day=2, P=healthy_dat_day2[,54], R=healthy_dat_day2[,129], T=healthy_dat_day2[,261]),
      data.frame(ids=healthy2_id, group=dis_gp, day=1, P=dat_disease[,54], R=dat_disease[,129], T=dat_disease[,261]),
      data.frame(ids=healthy2_id, group=dis_gp, day=2, P=dat_disease_day2[,54], R=dat_disease_day2[,129], T=dat_disease_day2[,261])
    )
    all_results$ids <- paste0(all_results$ids, all_results$group)
    
    #obtain median amplitudes 
    results_medians <- all_results %>% 
      group_by(ids, group, day) %>% 
      summarise(across(c(P, R, T), median, na.rm=TRUE), .groups = "drop")
    
    day1_meds <- filter(results_medians, day==1)
    
    # Compute ROC, p-value of non-parametric Wilcoxon test and ICC(A,1) and ICC(C,1) for P, R, T-amplitudes 
    amp_stats <- do.call(rbind, lapply(c("P", "R", "T"), function(p) {
      col_name <- paste0(p, "_amplitude")
      rbind(
        data.frame(change=change_type, sim=i, eval="ROC", parameter=col_name, value=roc(day1_meds$group, day1_meds[[p]], quiet=T)$auc[1]),
        data.frame(change=change_type, sim=i, eval="p-value", parameter=col_name, value=wilcox.test(day1_meds[[p]] ~ day1_meds$group)$p.value),
        data.frame(change=change_type, sim=i, eval="ICC_consistency", parameter=col_name, value=compute_ICC(results_medians, p)[1]),
        data.frame(change=change_type, sim=i, eval="ICC_agreement", parameter=col_name, value=compute_ICC(results_medians, p)[2])
      )
    }))
    
    
    #MFPCA scores ---------------------------------------------------------------
    predict_healthy <-  predict_mfpca_scores_wide(res, healthy_dat, paste0(healthy_id, "_healthy"), res$sigma2)
    predict_healthy_day2 <- predict_mfpca_scores_wide(res, healthy_dat_day2, paste0(healthy2_id, "_healthy"), res$sigma2)
    predict_disease <- predict_mfpca_scores_wide(res, dat_disease, paste0(healthy2_id, "_disease"), res$sigma2)
    predict_disease_day2 <- predict_mfpca_scores_wide(res, dat_disease_day2, paste0(healthy2_id, "_disease"), res$sigma2)
    
    day1_scores <- rbind(as.matrix(predict_healthy %>% dplyr::select(L1_PC1, L1_PC2, L1_PC3, L1_PC4, L1_PC5) %>% distinct()), 
                         as.matrix(predict_disease %>% dplyr::select(L1_PC1, L1_PC2, L1_PC3, L1_PC4, L1_PC5) %>% distinct()))
    
    # Extract subject-level scores for day 2-------------------------------------
    healthy_day2_L1 <- predict_healthy_day2 %>% dplyr::select(ID, L1_PC1, L1_PC2, L1_PC3, L1_PC4, L1_PC5) %>% distinct()
    healthy_day2_L1$day=2
    healthy_day2_L1$group="healthy"
    
    disease_day2_L1 <- predict_disease_day2 %>% dplyr::select(ID, L1_PC1, L1_PC2, L1_PC3, L1_PC4, L1_PC5) %>% distinct()
    disease_day2_L1$day=2
    disease_day2_L1$group="disease"
    
    healthy_day1_L1 <- predict_healthy %>% dplyr::select(ID, L1_PC1, L1_PC2, L1_PC3, L1_PC4, L1_PC5) %>% distinct()
    healthy_day1_L1$day=1
    healthy_day1_L1$group="healthy"
    
    disease_day1_L1 <-  predict_disease %>% dplyr::select(ID, L1_PC1, L1_PC2, L1_PC3, L1_PC4, L1_PC5) %>% distinct()
    disease_day1_L1$day=1
    disease_day1_L1$group="disease"
    
    # Combine in long format for ICC
    all_scores <-rbind(healthy_day2_L1, healthy_day1_L1, disease_day2_L1, disease_day1_L1)
    all_scores$ids <- all_scores$ID
    
    mfpca_stats <- do.call(rbind, lapply(1:5, function(pc) {
      p_name <- paste0("L1_PC", pc)
      label <- paste0("MFPCA_Score", pc)
      # Extract relevant day 1 data for ROC/Wilcox
      d1_scores <- filter(all_scores, day==1) 
      rbind(
        data.frame(change=change_type, sim=i, eval="ROC", parameter=label, value=roc(d1_scores$group, d1_scores[[p_name]], quiet=T)$auc[1]),
        data.frame(change=change_type, sim=i, eval="p-value", parameter=label, value=wilcox.test(d1_scores[[p_name]] ~ d1_scores$group)$p.value),
        data.frame(change=change_type, sim=i, eval="ICC_consistency", parameter=label, value=compute_ICC(all_scores, p_name)[1]),
        data.frame(change=change_type, sim=i, eval="ICC_agreement", parameter=label, value=compute_ICC(all_scores, p_name)[2])
      )
    }))
    
    
    #single fpca-------------------------------------------------------------
    #average healthy curves
    av_healthy_dat <- as.data.frame(healthy_dat)
    av_healthy_dat <- av_healthy_dat %>%
      mutate(ID = healthy_id)%>%
      group_by(ID) %>%
      summarise(across(.cols = everything(), .fns = ~ mean(.x)), .groups = "drop")
    
    #average healthy curves - day2
    av_healthy_dat_day2 <- as.data.frame(healthy_dat_day2)
    av_healthy_dat_day2 <- av_healthy_dat_day2 %>%
      mutate(ID = paste0(healthy_id, "_healthy")) %>%
      group_by(ID) %>%
      summarise(across(.cols = everything(), .fns = ~ mean(.x)), .groups = "drop")
    
    #average disease curves - day1
    av_dat_disease <- as.data.frame(dat_disease)
    av_dat_disease <- av_dat_disease %>%
      mutate(ID = paste0(healthy_id, "_disease")) %>%
      group_by(ID) %>%
      summarise(across(.cols = everything(), .fns = ~ mean(.x)), .groups = "drop")
    
    #average disease curves - day2
    av_dat_disease_day2 <- as.data.frame(dat_disease_day2)
    av_dat_disease_day2 <- av_dat_disease_day2 %>%
      mutate(ID = paste0(healthy_id, "_disease")) %>%
      group_by(ID) %>%
      summarise(across(.cols = everything(), .fns = ~ mean(.x)), .groups = "drop")
    
    
    predict_healthy <- predict_fpca_same_grid(fpca_single, av_healthy_dat[,-1])
    predict_healthy_day2 <- predict_fpca_same_grid(fpca_single, av_healthy_dat_day2[,-1])
    
    predict_disease <- predict_fpca_same_grid(fpca_single, av_dat_disease[,-1])
    predict_disease_day2 <-predict_fpca_same_grid(fpca_single, av_dat_disease_day2[,-1])
    
    
    fpca_scores <- rbind(predict_healthy, predict_disease)
    fpca_all_scores <- as.data.frame(fpca_scores)
    colnames(fpca_all_scores) <- c("L1_PC1", "L1_PC2", "L1_PC3",  "L1_PC4",  "L1_PC5")
    fpca_all_scores$day=1
    fpca_all_scores$ids  <- c(paste0(unique(healthy_id), "_healthy"), paste0(unique(healthy_id), "_disease"))
    
    
    
    predict_all_scores_day2 <- data.frame(rbind(predict_healthy_day2, predict_disease_day2))
    colnames(predict_all_scores_day2) <- c("L1_PC1", "L1_PC2", "L1_PC3",  "L1_PC4",  "L1_PC5")
    predict_all_scores_day2$day=2
    predict_all_scores_day2$ids  <- c(paste0(unique(healthy_id), "_healthy"), paste0(unique(healthy_id), "_disease"))
    
    all_scores <- rbind(fpca_all_scores, predict_all_scores_day2)
    all_scores$group= sub(".*_", "", all_scores$ids)
    
    fpca_stats <- do.call(rbind, lapply(1:5, function(pc) {
      p_name <- paste0("L1_PC", pc)
      label <- paste0("FPCA_Score", pc)
      # Extract relevant day 1 data for ROC/Wilcox
      d1_scores <- filter(all_scores, day==1) 
      rbind(
        data.frame(change=change_type, sim=i, eval="ROC", parameter=label, value=roc(d1_scores$group, d1_scores[[p_name]], quiet=T)$auc[1]),
        data.frame(change=change_type, sim=i, eval="p-value", parameter=label, value=wilcox.test(d1_scores[[p_name]] ~ d1_scores$group)$p.value),
        data.frame(change=change_type, sim=i, eval="ICC_consistency", parameter=label, value=compute_ICC(all_scores, p_name)[1]),
        data.frame(change=change_type, sim=i, eval="ICC_agreement", parameter=label, value=compute_ICC(all_scores, p_name)[2])
      )
    }))
    
    
    results_all <- rbind(results_all, amp_stats, mfpca_stats, fpca_stats) 
    
  }
  
}


saveRDS(results_all, "results_all.RDS")
