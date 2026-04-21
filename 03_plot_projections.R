#do fPCA on single ECGs and see if there are any groupings 
source("00_init.R")
source("02_sim_generate_data.R")
source("02_sim_functions.R")

load( "Analysis_Outputs/healthy_ECG_landmarked.RData")

Y_mat <- t(eval.fd(increments, regfd))
all_dat_name <- ecg_dat$Name

k=10
res <- mfpca.face(Y = Y_mat, id = all_dat_name, visit=NULL, npc = c(k,k))

mu  <- res$mu                        # Mean function
phi_b <- res$efunctions$level1          # Between eigenfunctions (T x Kb)   #308 x 3
phi_w <- res$efunctions$level2           # Within eigenfunctions (T x Kw)   #308 x 3
scores_b <- res$scores$level1            # Between scores (N_subjects x Kb)   #64 * 3
scores_w <- res$scores$level2            # Within scores (N_obs x Kw)   #1392 x 3
all_ids <- all_dat_name
id_unique <- unique(all_dat_name)             # Which subject each observation belongs to
time <- increments    


dat_disease <- generate_changes(Y_mat, change="flat_T_high", ids=healthy_id)
predict_disease <- predict_mfpca_scores_wide(res, dat_disease, paste0(all_dat_name, "_disease"), res$sigma2)
predict_disease <- predict_disease %>% dplyr::select(L1_PC1, L1_PC2)%>% distinct()



dat_disease2 <- generate_changes(Y_mat, change="reduce_all_high", ids=healthy_id)
predict_disease2 <- predict_mfpca_scores_wide(res, dat_disease2, paste0(all_dat_name, "_disease"), res$sigma2)
predict_disease2 <- predict_disease2 %>% dplyr::select(L1_PC1, L1_PC2)%>% distinct()


plot_scores <- ggplot(data = data.frame(score1=scores_b[,1], 
                                        score2=scores_b[,2]), 
                      aes(x=score1, y=score2)) + theme_bw(base_size = 11)+geom_point(aes(color="Healthy"))+ 
  geom_point(aes(x=predict_disease$L1_PC1, y= predict_disease$L1_PC2, color="Flattened T-wave"))+
  geom_point(aes(x=predict_disease2$L1_PC1, y= predict_disease2$L1_PC2, color="All Waves Flattened"))+
  xlab("Subject-level FPC1 score") + ylab("Subject-level FPC2 score") + 
  ggtitle("Subject-level scores and Projections") + scale_color_manual(values = c(
    "Flattened T-wave" = "darkorange", 
    "All Waves Flattened" = "#b5367a",
    "Healthy" = "black"
  ))+labs(color="Group")+  guides(color = guide_legend(ncol=2))+
  theme(legend.position="bottom")


ggsave("Simulation_plots/plot_projections.png", plot = plot_scores, width = 4, height = 4, dpi = 200)
