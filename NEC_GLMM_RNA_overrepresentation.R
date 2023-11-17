


RPKM_args_LA_merge$Cohort<-as.factor(RPKM_args_LA_merge$Cohort)

model.b = lme(Abundance ~ Cohort*DOL, 
              random = ~1|Patient, data=RPKM_args_LA_merge)
summary(model.b)

model.null = lme(Abundance ~ DOL, 
                 random = ~1|Patient, data=RPKM_args_LA_merge)

anova(update(model.b, . ~ ., method = 'ML'),
      update(model.null, . ~ ., method = 'ML'))

library (multcomp)
summary(glht(model.b, mcp(Cohort='Tukey')))


library('lme4')
library('lmerTest')
library('MuMIn')
library('multcomp')
library(dplyr)

GLMM_RNA_DNA<-read.delim("./Desktop/Projects/NEC/final_analysis/RNA/taxonomy/GLMM_RNA_DNA/211109_RNADNAmetaphlan_latecasevscontrol_forRTGLMM.txt", header = TRUE,check.names = FALSE)

ggplot(GLMM_RNA_DNA, aes(x=dayspreonset, y=log(RNA_DNA_FC), fill=NEC_status)) +geom_point(size=5, pch=21)+
  theme_classic()  +geom_smooth(method='lm')+facet_wrap(~DNA_variable, scales="free") 

GLMM_RNA_DNA<-read.delim("./Desktop/Projects/NEC/final_analysis/RNA/taxonomy/GLMM_RNA_DNA/211109_RNADNAmetaphlan_earlycasevslatecase_forRTGLMM_noZeros.txt", header = TRUE,check.names = FALSE)

ggplot(GLMM_RNA_DNA, aes(x=dayspreonset_log, y=log(RNA_DNA_FC), fill=Onset)) +geom_point(size=5, pch=21)+
  theme_classic()  +geom_smooth(method='lm')+facet_wrap(~DNA_variable, scales="free") 

GLMM_RNA_DNA<-read.delim("./Desktop/Projects/NEC/final_analysis/RNA/taxonomy/GLMM_RNA_DNA/211110_GLMM_RNA_DNA_no_zeros.txt", header = TRUE,check.names = FALSE)

ggplot(GLMM_RNA_DNA, aes(x=dayspreonset, y=log(RNA_DNA_FC), fill=NEC_status)) +geom_point(size=5, pch=21)+
  theme_classic()  +geom_smooth(method='lm')+facet_wrap(~DNA_variable, scales="free") 

#Staph epidermidis
GLMM_RNA_DNA_Staphylococcus_epidermidis<-subset(GLMM_RNA_DNA,GLMM_RNA_DNA$DNA_variable=="Staphylococcus_epidermidis")
#GLMM_RNA_DNA_Staphylococcus_epidermidis[sapply(GLMM_RNA_DNA_Staphylococcus_epidermidis, is.character)] <- lapply(GLMM_RNA_DNA_Staphylococcus_epidermidis[sapply(GLMM_RNA_DNA_Staphylococcus_epidermidis, is.character)], 
                                                                                                                 as.factor)
#cat <- sapply(GLMM_RNA_DNA_Staphylococcus_epidermidis, is.factor)
#GLMM_RNA_DNA_Staphylococcus_epidermidis2<-Filter(function(x) nlevels(x)>1, GLMM_RNA_DNA_Staphylococcus_epidermidis[,cat])

#cat <- sapply(GLMM_RNA_DNA_Staphylococcus_epidermidis, is.numeric)
#GLMM_RNA_DNA_Staphylococcus_epidermidis_numeric<-GLMM_RNA_DNA_Staphylococcus_epidermidis[,cat]
#GLMM_RNA_DNA_Staphylococcus_epidermidis_numeric$Sample_ID<-GLMM_RNA_DNA_Staphylococcus_epidermidis$Sample_ID
#GLMM_RNA_DNA_Staphylococcus_epidermidis_filtered<-merge(GLMM_RNA_DNA_Staphylococcus_epidermidis_numeric, GLMM_RNA_DNA_Staphylococcus_epidermidis2,by.x = 'Sample_ID', by.y = 'Sample_ID')

#GLMM_RNA_DNA_Staphylococcus_epidermidis[,colnames(GLMM_RNA_DNA_Staphylococcus_epidermidis)%ni%colnames(GLMM_RNA_DNA_Staphylococcus_epidermidis_filtered)]

str(GLMM_SD_acute)

FIRSTPASSMODEL.glm <-lmer(log(RNA_DNA_FC)~ (1|Patient)+(0 + DOL_log | Patient) +NEC_status+dayspreonset+GA_birth_log+ABX_recent+Recent_penicillin+Recent_carbapenem+
                             Recent_amingoglycoside+Recent_macrolide+Recent_glycopeptide+Recent_1st_gen_cephalosporin+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
                             intraventricular_hemorrhage+Abx_cumulative_log+ANTIM_cumulative_log+mat_ppbmi_log+
                             mat_inpatient_meds+milk_cumulative_log+milk_current+formula_cumulative+inf+mat_alcohol,
                           REML=FALSE,
                           data = GLMM_RNA_DNA_Staphylococcus_epidermidis)

#2B. run Step function- initial backfitting
lmerTest::step(FIRSTPASSMODEL.glm,
               type = 3,
               alpha.random = 0.1,
               alpha.fixed = 0.05,
               reduce.fixed = TRUE,
               keep = c("dayspreonset","NEC_status"),
               fixed.calc = TRUE,
               lsmeans.calc = TRUE,
               difflsmeans.calc = TRUE)


BACKFITTEDMODEL.glm <-lmer(RNA_DNA_FC_final~ (1|Patient)+(0 + DOL_log | Patient) +NEC_status+dayspreonset+ABX_recent+Recent_amingoglycoside+
                             inf+milk_current+mat_inpatient_meds,
                           REML=FALSE,
                           data = GLMM_RNA_DNA_Staphylococcus_epidermidis)

#pseudo-R2
r.squaredGLMM(BACKFITTEDMODEL.glm)
#summary. Look at correlation matrix. If any significant interactions with time (<-0.1 or >0.1), add to the list of variables in the first past model in part A, rerun Step. 
summary(BACKFITTEDMODEL.glm)

#2D once model has been selected, correct for multiple comparisons
g <- glht(BACKFITTEDMODEL.glm, lincfit = mcp(tension = "Tukey"))
summary(g)

ggplot(GLMM_RNA_DNA_Escherichia_coli, aes(x=dayspreonset, y=log(RNA_DNA_FC), fill=NEC_status)) +geom_point(size=5, pch=21)+
  theme_classic()  +geom_smooth(method='lm')

#Enterococcus_faecalis
GLMM_RNA_DNA_Enterococcus_faecalis<-subset(GLMM_RNA_DNA,GLMM_RNA_DNA$DNA_variable=="Enterococcus_faecalis")
#GLMM_RNA_DNA_Staphylococcus_epidermidis[sapply(GLMM_RNA_DNA_Staphylococcus_epidermidis, is.character)] <- lapply(GLMM_RNA_DNA_Staphylococcus_epidermidis[sapply(GLMM_RNA_DNA_Staphylococcus_epidermidis, is.character)], 
as.factor)
#cat <- sapply(GLMM_RNA_DNA_Staphylococcus_epidermidis, is.factor)
#GLMM_RNA_DNA_Staphylococcus_epidermidis2<-Filter(function(x) nlevels(x)>1, GLMM_RNA_DNA_Staphylococcus_epidermidis[,cat])

#cat <- sapply(GLMM_RNA_DNA_Staphylococcus_epidermidis, is.numeric)
#GLMM_RNA_DNA_Staphylococcus_epidermidis_numeric<-GLMM_RNA_DNA_Staphylococcus_epidermidis[,cat]
#GLMM_RNA_DNA_Staphylococcus_epidermidis_numeric$Sample_ID<-GLMM_RNA_DNA_Staphylococcus_epidermidis$Sample_ID
#GLMM_RNA_DNA_Staphylococcus_epidermidis_filtered<-merge(GLMM_RNA_DNA_Staphylococcus_epidermidis_numeric, GLMM_RNA_DNA_Staphylococcus_epidermidis2,by.x = 'Sample_ID', by.y = 'Sample_ID')

#GLMM_RNA_DNA_Staphylococcus_epidermidis[,colnames(GLMM_RNA_DNA_Staphylococcus_epidermidis)%ni%colnames(GLMM_RNA_DNA_Staphylococcus_epidermidis_filtered)]

str(GLMM_SD_acute)

BACKFITTEDMODEL.glm <-lmer(RNA_DNA_FC_final~ (1|Patient)+(0 + DOL_log | Patient) +NEC_status+dayspreonset+GA_birth_log+ABX_recent+Recent_penicillin+Recent_carbapenem+
                             Recent_amingoglycoside+Recent_macrolide+Recent_glycopeptide+Recent_1st_gen_cephalosporin+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
                             intraventricular_hemorrhage+Abx_cumulative_log+ANTIM_cumulative_log+mat_ppbmi_log+
                             mat_inpatient_meds+milk_cumulative_log+milk_current+formula_cumulative+inf+mat_alcohol,
                           REML=FALSE,
                           data = GLMM_RNA_DNA_Enterococcus_faecalis)

#pseudo-R2
r.squaredGLMM(BACKFITTEDMODEL.glm)
#summary. Look at correlation matrix. If any significant interactions with time (<-0.1 or >0.1), add to the list of variables in the first past model in part A, rerun Step. 
summary(BACKFITTEDMODEL.glm)

#2D once model has been selected, correct for multiple comparisons
g <- glht(BACKFITTEDMODEL.glm, lincfit = mcp(tension = "Tukey"))
summary(g)


#Escherichia_coli
GLMM_RNA_DNA_Escherichia_coli<-subset(GLMM_RNA_DNA,GLMM_RNA_DNA$DNA_variable=="Escherichia_coli")
GLMM_RNA_DNA_Escherichia_coli[sapply(GLMM_RNA_DNA_Escherichia_coli, is.character)] <- lapply(GLMM_RNA_DNA_Escherichia_coli[sapply(GLMM_RNA_DNA_Escherichia_coli, is.character)], 
as.factor)
cat <- sapply(GLMM_RNA_DNA_Escherichia_coli, is.factor)
GLMM_RNA_DNA_Escherichia_coli2<-Filter(function(x) nlevels(x)>1, GLMM_RNA_DNA_Escherichia_coli[,cat])

cat <- sapply(GLMM_RNA_DNA_Escherichia_coli, is.numeric)
GLMM_RNA_DNA_Escherichia_coli_numeric<-GLMM_RNA_DNA_Escherichia_coli[,cat]
GLMM_RNA_DNA_Escherichia_coli_numeric$Sample_ID<-GLMM_RNA_DNA_Escherichia_coli$Sample_ID
GLMM_RNA_DNA_Escherichia_coli_filtered<-merge(GLMM_RNA_DNA_Escherichia_coli_numeric, GLMM_RNA_DNA_Escherichia_coli2,by.x = 'Sample_ID', by.y = 'Sample_ID')
colnames(GLMM_RNA_DNA_Escherichia_coli_filtered)

#GLMM_RNA_DNA_Staphylococcus_epidermidis[,colnames(GLMM_RNA_DNA_Staphylococcus_epidermidis)%ni%colnames(GLMM_RNA_DNA_Staphylococcus_epidermidis_filtered)]

str(GLMM_SD_acute)

#2A. Make First Past Model
FIRSTPASSMODEL.glm <- lmer(log(RNA_DNA_FC)~ (1|Patient)+NEC_status+dayspreonset+GA_birth_log+ABX_recent+Recent_penicillin+Recent_carbapenem+
                             Recent_amingoglycoside+Recent_glycopeptide+Recent_1st_gen_cephalosporin+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
                             intraventricular_hemorrhage+Abx_cumulative_log+ANTIM_cumulative_log+mat_ppbmi_log+
                             milk_cumulative_log+milk_current+formula_cumulative+inf,
                           REML=FALSE,
                           data = GLMM_RNA_DNA_Escherichia_coli)

#2B. run Step function- initial backfitting
lmerTest::step(FIRSTPASSMODEL.glm,
               type = 3,
               alpha.random = 0.1,
               alpha.fixed = 0.05,
               reduce.fixed = TRUE,
               keep = c("dayspreonset","NEC_status"),
               fixed.calc = TRUE,
               lsmeans.calc = TRUE,
               difflsmeans.calc = TRUE)


BACKFITTEDMODEL.glm <-lmer(RNA_DNA_FC_final~ (1|Patient)+(0 + DOL_log | Patient) +NEC_status+dayspreonset+intraventricular_hemorrhage+Recent_carbapenem+
                             Recent_3rd_gen_cephalosporin,
                           REML=FALSE,
                           data = GLMM_RNA_DNA_Escherichia_coli)

#pseudo-R2
r.squaredGLMM(BACKFITTEDMODEL.glm)
#summary. Look at correlation matrix. If any significant interactions with time (<-0.1 or >0.1), add to the list of variables in the first past model in part A, rerun Step. 
summary(BACKFITTEDMODEL.glm)

#2D once model has been selected, correct for multiple comparisons
g <- glht(BACKFITTEDMODEL.glm, lincfit = mcp(tension = "Tukey"))
summary(g)


#Klebsiella_pneumoniae
GLMM_RNA_DNA_Klebsiella_pneumoniae<-subset(GLMM_RNA_DNA,GLMM_RNA_DNA$DNA_variable=="Klebsiella_pneumoniae")
GLMM_RNA_DNA_Klebsiella_pneumoniae[sapply(GLMM_RNA_DNA_Klebsiella_pneumoniae, is.character)] <- lapply(GLMM_RNA_DNA_Klebsiella_pneumoniae[sapply(GLMM_RNA_DNA_Klebsiella_pneumoniae, is.character)], 
                                                                                             as.factor)
cat <- sapply(GLMM_RNA_DNA_Klebsiella_pneumoniae, is.factor)
GLMM_RNA_DNA_Klebsiella_pneumoniae2<-Filter(function(x) nlevels(x)>1, GLMM_RNA_DNA_Klebsiella_pneumoniae[,cat])

cat <- sapply(GLMM_RNA_DNA_Escherichia_coli, is.numeric)
GLMM_RNA_DNA_Escherichia_coli_numeric<-GLMM_RNA_DNA_Escherichia_coli[,cat]
GLMM_RNA_DNA_Escherichia_coli_numeric$Sample_ID<-GLMM_RNA_DNA_Escherichia_coli$Sample_ID
GLMM_RNA_DNA_Escherichia_coli_filtered<-merge(GLMM_RNA_DNA_Escherichia_coli_numeric, GLMM_RNA_DNA_Escherichia_coli2,by.x = 'Sample_ID', by.y = 'Sample_ID')
colnames(GLMM_RNA_DNA_Escherichia_coli_filtered)

#GLMM_RNA_DNA_Staphylococcus_epidermidis[,colnames(GLMM_RNA_DNA_Staphylococcus_epidermidis)%ni%colnames(GLMM_RNA_DNA_Staphylococcus_epidermidis_filtered)]

str(GLMM_SD_acute)

BACKFITTEDMODEL.glm <-lmer(RNA_DNA_FC_final~ (1|Patient)+(0 + DOL_log | Patient) +NEC_status+dayspreonset+GA_birth_log+ABX_recent+Recent_penicillin+Recent_carbapenem+
                             Recent_amingoglycoside+Recent_macrolide+Recent_glycopeptide+Recent_1st_gen_cephalosporin+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
                             intraventricular_hemorrhage+Abx_cumulative_log+ANTIM_cumulative_log+mat_ppbmi_log+
                             mat_inpatient_meds+milk_cumulative_log+milk_current+formula_cumulative+inf+mat_alcohol,
                           REML=FALSE,
                           data = GLMM_RNA_DNA_Klebsiella_pneumoniae)

#pseudo-R2
r.squaredGLMM(BACKFITTEDMODEL.glm)
#summary. Look at correlation matrix. If any significant interactions with time (<-0.1 or >0.1), add to the list of variables in the first past model in part A, rerun Step. 
summary(BACKFITTEDMODEL.glm)

#2D once model has been selected, correct for multiple comparisons
g <- glht(BACKFITTEDMODEL.glm, lincfit = mcp(tension = "Tukey"))
summary(g)
