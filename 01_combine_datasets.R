# Project: MFPCA digital outcomes 
# Script Purpose: extract left knee flexion/extension data from open source gait data
# Author: Mia Tackney 
# Date: 2026-04-15
# GitHub: https://github.com/mst1g15/MFPCA_digital_outcomes

#load packages 
source("00_init.R")

#Parkinson's Disease Dataset----------------------------------------------------

#This is obtained from: 
#Boari, Daniel (2021). A dataset of overground walking full-body kinematics and 
#kinetics in individuals with Parkinson’s disease. figshare. Dataset. 
#https://doi.org/10.6084/m9.figshare.14896881.v4

#Obtain left Knee flexion/extension data
files <- list.files(
  path = "Datasets/Gait cycles/",
  pattern = "Left_Knee.*\\.csv$")  

all_data <- c()

for (i in files){
  
  data <- read_excel(paste0("Datasets/Gait cycles/", i), col_names=TRUE)
  data <- data[-1,]
  
  #select knee flx/extension 
  data <- data %>% dplyr::select(`Gait cycle [%]`, contains("Knee Flx/Extension"))
  
  #obtain additional variables 
  meta_data <- data.frame(filename=i)%>%
    mutate(filename = sub("\\.csv$", "", filename)) %>%
    separate(
      col = filename,
      into = c("Subject", "Medication", "Task", "Joint", "Type", "Angle"),
      sep = "_"
    )
  
  data_longer <- data %>% pivot_longer(cols=contains("Knee Flx/Extension"), 
                                       names_to="Stride", 
                                       values_to="value")
  
  data_longer$Stride <- str_remove(data_longer$Stride , "Knee Flx/Extension\\.*\\s*")
  data_longer$value <- as.numeric(data_longer$value)
  data_longer$Stride <- as.factor(data_longer$Stride)
  
  data_longer <- cbind(meta_data, data_longer)
  
  all_data <- rbind(all_data, data_longer)
  
}

saveRDS(all_data, "Datasets/left_knee_data.RDS")

#Healthy Dataset ---------------------------------------------------------------

#This is obtained from: 
#Helwig, N. & Hsiao-Wecksler, E. (2016). Multivariate Gait Data. Dataset. 
#UCI Machine Learning Repository. https://doi.org/10.24432/C5861T.

#obtain knee flexion/extension data 
gait <- read.csv("Datasets_healthy_UCI/gait.csv")

knee <- gait %>% filter(condition==1 &  #unbraced setting 
                          leg==1 &        #left side 
                          joint==2)

colnames(knee) <- c("Subject", "Condition", "Stride",  "Leg", "Joint", "Gait cycle [%]", "Angle")

#Gait cycle on this dataset starts at zero instead of at one (which is the case for the PD dataset)
knee$`Gait cycle [%]` <- knee$`Gait cycle [%]` +1
saveRDS(knee, "Datasets/healthy_knee_data.RDS")
