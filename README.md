# Multilevel Functional Principal Components Analysis (MFPCA) for digital outcome measures

Welcome! This repository contains R scripts to run MFPCA to obtain low-dimensional summary scalars, which may be used as digital outcome measure in a trial setting. 


## Datasets

Three real-world datasets are part of the project. 

### 1. ECG data from Apple Watch Electrocardiograms 
Apple Watch (AW) ECG data from healthy participants from a registry study aimed to assess the validity of smartwatch-based ECG recordings are used. This data is not publicly available. 

### 2. Gait data from individuals with Parkinson's disease
This is an open source dataset obtained from: 
Boari, Daniel (2021). A dataset of overground walking full-body kinematics and kinetics in individuals with Parkinson’s disease. figshare. Dataset. https://doi.org/10.6084/m9.figshare.14896881.v4

### 3. Gait data from healthy individuals 
This is an open source dataset obtained from:
Helwig, N. & Hsiao-Wecksler, E. (2016). Multivariate Gait Data. Dataset. UCI Machine Learning Repository. https://doi.org/10.24432/C5861T.

## R code for Simulation study based on Apple Watch ECG data 

TO DO 

## R code for analysis of gait data
The repository contains the following R scripts: 

- `00_init.R` loads packages.
- `01_combine_datasets.R` extracts left knee flexion/extension data from open source gait data
- `02_projection_scores.R` contains functions to compute MFPCA projection scores

## Vignette 
A vignette is provided for the analysis of knee flexion-extension data from individuals with PD.

- `PD MFPCA vignette.qmd` is a Quarto markdown file for the vignette 
- `PD MFPCA vignette.qmd` is an html file for the vignette



