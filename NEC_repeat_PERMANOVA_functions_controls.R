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
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_PERMANOVA/ALL_control_human_community.txt", header = TRUE,row.names = 1)

SM_DF_BC <-vegdist(SM_DF, method = "bray")

#metadata<-read.table('Desktop/Projects/NEC/final_analysis/control_analysis/tax_PERMANOVA/Patient_IDs.txt', header = TRUE, row.names = 1)

#adonis(SM_DF_BC~Patient_ID,data=metadata,method="bray",permutations=999)
#######if this fails it is due to loaded package incomptability 


study_2<-read.table('Desktop/Projects/NEC/final_analysis/control_analysis/DNA_PERMANOVA/ALL_study.txt', header = FALSE)
study_2<-study_2$V1
study_longitudinal_2<-read.table('Desktop/Projects/NEC/final_analysis/control_analysis/DNA_PERMANOVA/ALL_study_longitudinal.txt', header = FALSE)
study_longitudinal<-as.logical(study_longitudinal_2$V2)
names(study_longitudinal)<-study_longitudinal_2$V1
subject_2<-read.table('Desktop/Projects/NEC/final_analysis/control_analysis/DNA_PERMANOVA/ALL_patients.txt', header = FALSE)
subject<-subject_2$V1
subject_data<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/", header = TRUE,check.names = FALSE, row.names=1,na.strings=c(""))
sample_data<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_PERMANOVA/ALL_sample/sample_DOL_ALL.txt", header = TRUE,check.names = FALSE, row.names=1)


#Sample_data
PERMANOVA_repeat_measures(
  SM_DF_BC,
  subject,
  sample_data = sample_data,
  metadata_order = c(names(sample_data)),
  permutations=999, ncores=1)

#Subject_data
PERMANOVA_repeat_measures(
  SM_DF_BC,
  subject, subject_data = subject_data,
  metadata_order = c(names(subject_data)),
  permutations=999, ncores=1)




PERMANOVA_repeat_measures_meta(
  SM_DF_BC,
  study_2, study_longitudinal,
  subject,
  sample_data = sample_data,
  metadata_order = c(names(sample_data)),
  permutations=999, ncores=1)


PERMANOVA_repeat_measures_meta(
  SM_DF_BC,
  study_2, study_longitudinal,
  subject,
  subject_data = subject_data,
  metadata_order = c(names(subject_data)),
  permutations=999, ncores=1)












library(corrplot)
library(RColorBrewer)
library(ggplot2)


###ALL
medications_all<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/metadata/Meds_correlations/medications_correlations.txt", header = TRUE,check.names = FALSE, row.names=1)
medications_all_BC <-vegdist(medications_all, method = "bray")
medications_all_BC<-as.matrix(medications_all_BC)
write.csv(medications_10_BC,"./Desktop/Projects/NEC/final_analysis/control_analysis/metadata/medications_jaccard.csv")

corrplot(medications_all_BC,method="color",type="lower",cl.lim = c(0, 1),col=brewer.pal(n=8, name="RdYlBu"),
         is.corr = FALSE,diag=FALSE,outline = TRUE,
         tl.pos="ld",tl.col="black")

####DOL 10
medications_10<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/metadata/Meds_correlations/medications_correlations_10_DOL.txt", header = TRUE,check.names = FALSE, row.names=1)
medications_10_BC <-vegdist(medications_10, method = "bray")
medications_10_BC<-as.matrix(medications_10_BC)

corrplot(medications_10_BC,method="color",type="lower",cl.lim = c(0, 1),col=brewer.pal(n=8, name="RdYlBu"),
         is.corr = FALSE,diag=FALSE,outline = TRUE,
         tl.pos="ld",tl.col="black",na.label="square",na.label.col="#818181")


####DOL 20
medications_20<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/metadata/Meds_correlations/medications_correlations_20_DOL.txt", header = TRUE,check.names = FALSE, row.names=1)
medications_20_BC <-vegdist(medications_20, method = "bray")
medications_20_BC<-as.matrix(medications_20_BC)

corrplot(medications_20_BC,method="color",type="lower",cl.lim = c(0, 1),col=brewer.pal(n=8, name="RdYlBu"),
         is.corr = FALSE,diag=FALSE,outline = TRUE,
         tl.pos="ld",tl.col="black",na.label="square",na.label.col="#818181")

####DOL 30
medications_30<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/metadata/Meds_correlations/medications_correlations_30_DOL.txt", header = TRUE,check.names = FALSE, row.names=1)
medications_30_BC <-vegdist(medications_30, method = "bray")
medications_30_BC<-as.matrix(medications_30_BC)

corrplot(medications_30_BC,method="color",type="lower",cl.lim = c(0, 1),col=brewer.pal(n=8, name="RdYlBu"),
         is.corr = FALSE,diag=FALSE,outline = TRUE,
         tl.pos="ld",tl.col="black",na.label="square",na.label.col="#818181")


####DOL 40
medications_40<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/metadata/Meds_correlations/medications_correlations_40_DOL.txt", header = TRUE,check.names = FALSE, row.names=1)
medications_40_BC <-vegdist(medications_40, method = "bray")
medications_40_BC<-as.matrix(medications_40_BC)

corrplot(medications_40_BC,method="color",type="lower",cl.lim = c(0, 1),col=brewer.pal(n=8, name="RdYlBu"),
         is.corr = FALSE,diag=FALSE,outline = TRUE,
         tl.pos="ld",tl.col="black",na.label="square",na.label.col="#818181")


####DOL 50
medications_50<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/metadata/Meds_correlations/medications_correlations_50_DOL.txt", header = TRUE,check.names = FALSE, row.names=1)
medications_50_BC <-vegdist(medications_50, method = "bray")
medications_50_BC<-as.matrix(medications_50_BC)

corrplot(medications_50_BC,method="color",type="lower",cl.lim = c(0, 1),col=brewer.pal(n=8, name="RdYlBu"),
         is.corr = FALSE,diag=FALSE,outline = TRUE,
         tl.pos="ld",tl.col="black",na.label="square",na.label.col="#818181")


####DOL 60
medications_60<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/metadata/Meds_correlations/medications_correlations_60_DOL.txt", header = TRUE,check.names = FALSE, row.names=1)
medications_60_BC <-vegdist(medications_60, method = "bray")
medications_60_BC<-as.matrix(medications_60_BC)

corrplot(medications_60_BC,method="color",type="lower",cl.lim = c(0, 1),col=brewer.pal(n=8, name="RdYlBu"),
         is.corr = FALSE,diag=FALSE,outline = TRUE,
         tl.pos="ld",tl.col="black",na.label="square",na.label.col="#818181")



