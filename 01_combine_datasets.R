# Project: MFPCA digital outcomes 
# Script Purpose: extract left knee flexion/extension data from open source gait data
# Author: Mia Tackney 
# Date: 2026-11-17
# GitHub: https://github.com/mst1g15/MFPCA_digital_outcomes

#load packages 
source("00_init.R")

#Parkinson's Disease Dataset----------------------------------------------------

#This is obtained from: 
#Daniel Boari(2021). A dataset of overground walking full-body kinematics and 
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
#Claudiane Fukuchi, Reginaldo Fukuchi, and Marcos Duarte. A public data set of overground and
#treadmill walking kinematics and kinetics of healthy individuals, 2018. URL 
#https://figshare.com/articles/dataset/A_public_data_set_of_overground_and_treadmill_walking_kinematics_and_kinetics_of_healthy_individuals/5722711


files <- list.files(
  path = "Datasets_Healthy_Fukuchi/51subjs/",
  pattern = "ang.*\\.txt$")  

all_data <- c()

for (i in files){
  
  data <- read.delim(paste0("Datasets_Healthy_Fukuchi/51subjs/", i))

  #select knee flx/extension 
  data <- data %>% dplyr::select(Time, LKneeAngleZ)
  data$Subject <- as.integer(str_extract(i, "(?<=WBDS)\\d+"))
  data$Stride <- str_extract(i, "(?<=walk)[A-Za-z0-9]+(?=ang)")
  
  all_data <- rbind(all_data, data)
  
}

colnames(all_data) <- c("Gait cycle [%]", "Angle", "Subject", "Stride")
saveRDS(all_data, "Datasets/Fukuchi_left_knee.RDS")

#plot results 
ggplot(all_data, aes(x=`Gait cycle [%]`, y=Angle, group=Stride)) + geom_line() + facet_wrap(.~Subject)
