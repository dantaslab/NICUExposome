###### Packages ######
library(ggplot2)
library(reshape)
library(labdsv)
library(vegan)
library(ggpubr)
library(permute)
library(dplyr)
library(rsample)
library(purrr)
library(tidyr)
library(rowr)
library(rcompanion)
library(ggpubr)
library("BiodiversityR")
library("Maaslin2")
library("randomForest")

###### Load and edit DF #######
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/211006_rel_pathabundances.txt", header = TRUE,row.names = 1,check.names=FALSE)
SM_DF2<-t(SM_DF)
metadata<-read.delim('Desktop/Projects/NEC/final_analysis/taxonomy/210429_nec_metadata.txt', header = TRUE)
SM_DF_selection <- SM_DF2[row.names(SM_DF2) %in% metadata$SampleID, ]
SM_DF_BC <-vegdist(SM_DF_selection, method = "bray")

metadata<-read.table('Desktop/Projects/NEC/final_analysis/taxonomy/PERMANOVA/Patient_IDs.txt', header = TRUE, row.names = 1)

adonis(SM_DF_BC~Patient_ID,data=metadata,method="bray",permutations=999)
#######if this fails it is due to loaded package incomptability 

subject_2<-read.table('Desktop/Projects/NEC/final_analysis/taxonomy/PERMANOVA/Patients_BMI.txt', header = FALSE)
subject<-subject_2$V1
subject_data<-read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/PERMANOVA/subject/subject_matPPBMI.txt", header = TRUE,check.names = FALSE, row.names=1)
sample_data<-read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/PERMANOVA/sample/sample_recent_pen.txt", header = TRUE,check.names = FALSE, row.names=1)

#Subject_data
PERMANOVA_repeat_measures(
  SM_DF_BC,
  subject, subject_data = subject_data,
  metadata_order = c(names(subject_data)),
  permutations=999, ncores=1)

#Sample_data
PERMANOVA_repeat_measures(
  SM_DF_BC,
  subject,
  sample_data = sample_data,
  metadata_order = c(names(sample_data)),
  permutations=999, ncores=1)

#ALL
PERMANOVA_repeat_measures(
  SM_DF_BC,
  subject, subject_data = subject_data,
  sample_data = sample_data,
  metadata_order = c(names(subject_data), names(sample_data)),
  permutations=999, ncores=1)

