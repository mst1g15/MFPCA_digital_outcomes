# Project: MFPCA digital outcomes 
# Script Purpose: plot results of simulations   
# Author: Mia Tackney 
# Date: 2026-04-15
# GitHub: https://github.com/mst1g15/MFPCA_digital_outcomes


source("00_init.R")

results_all <- readRDS("results_all_14Apr2026.RDS")

nsim=5000

#tidy up column names for manuscript figure 

results_all$value <- as.numeric(results_all$value)
results_all$change <- as.factor(results_all$change)

levels(results_all$change) <- c("Changing T-wave", "Flattened T-wave", "No change", "All Waves Flattened", "ST elevation")
results_all$change <- factor(results_all$change, levels=c("No change", "Flattened T-wave","Changing T-wave",  "All Waves Flattened", "ST elevation"))
results_all$parameter <- as.factor(results_all$parameter)
levels(results_all$parameter) <- c( "P amplitude",  "R amplitude",  "T amplitude", 
                                    "FPCA Score1",  "FPCA Score2",  "FPCA Score3",  "FPCA Score4",  "FPCA Score5",
                                    "MFPCA Score1", "MFPCA Score2", "MFPCA Score3", "MFPCA Score4", "MFPCA Score5")

results_all$eval <- as.factor(results_all$eval)
levels(results_all$eval) <-  c("ICC(A,1)", "ICC(C,1)", "p-value", "ROC")

results_all <- results_all %>% 
  mutate(Metric = case_when(
    parameter %in% c("P amplitude",  "R amplitude",  "T amplitude")~ "scalar summary",
    parameter %in% c("FPCA Score1",  "FPCA Score2",  "FPCA Score3",  "FPCA Score4",  "FPCA Score5") ~ "FPCA Score",
    parameter %in% c("MFPCA Score1", "MFPCA Score2", "MFPCA Score3", "MFPCA Score4", "MFPCA Score5") ~ "MFPCA Score"
  ))

results_all_plot <- results_all %>% group_by(eval, parameter, change) %>% mutate(result_mean= mean(value, na.rm=TRUE))%>%
  dplyr::select(-c(value, sim)) %>% distinct()

results_all_plot <- results_all_plot %>% mutate(result_sd = sqrt(result_mean*(1-result_mean)))

#do not need to present 5th MFPA score as first four eigenfunctions explain over 95% of variability 
results_all_plot <- results_all_plot %>% dplyr::filter(!(parameter %in% c("MFPCA Score5", "FPCA Score5")))


p <- ggplot(results_all_plot, aes(x=parameter, y=result_mean, color=Metric)) + 
  geom_point(size=0.8) + 
  geom_errorbar(aes(ymin=result_mean-1.96*result_sd/sqrt(nsim), 
                    ymax=result_mean+1.96*result_sd/sqrt(nsim)),  width=0.5)+
  facet_grid(eval~change, scales="free_y")+theme_bw() + 
  xlab("Summary Metric") + ylab("Mean of metric across simulations")+
  scale_color_manual(values = c( "#d95f02", "#b5367a", "#000000"))+
  theme(
    strip.text.y = element_text(angle = 0), # 0 makes text horizontal
    strip.background = element_rect(fill = "white", color = "black"), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 50, hjust = 1, size=8,  color = "black")
  )

ggsave("Simulation_plots/initial_results.png", plot = p, width = 8, height = 5.5, dpi = 200)



