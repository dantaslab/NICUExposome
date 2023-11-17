#####################################
####### NEC GLMM SD fitting #########
#####################################

library('lme4')
library('lmerTest')
library('MuMIn')
library('multcomp')

SM_DF <- read.delim('./Desktop/Projects/NEC/final_analysis/humann/211006_rel_pathabundances.txt', header = TRUE, row.names = 1,check.names = FALSE)
metadata<-read.delim('Desktop/Projects/NEC/final_analysis/humann/210920_NEC_metadata_reduced_acute.txt', header = TRUE)
SM_DF2<-t(SM_DF)
SM_DF_lifetime <- SM_DF2[row.names(SM_DF2) %in% metadata$Sample_ID, ]
shannDiv_humann<-as.data.frame(diversity(SM_DF_lifetime, index='shannon'))
shannDiv_humann$Sample_ID<-row.names(shannDiv_humann)
shannDivMetaData<-merge(shannDiv_humann, metadata,by.x = 'Sample_ID', by.y = 'Sample_ID')
names(shannDivMetaData)[2] <- 'ShannonDiversity'
row.names(shannDivMetaData)<-shannDivMetaData[,1]
shannDivMetaData<-shannDivMetaData[,-1]
write.table(shannDivMetaData, './Desktop/Projects/NEC/final_analysis/humann/210921_SD_Meta_acute.txt', sep='\t')


GLMM_SD_acute<-read.delim("./Desktop/Projects/NEC/final_analysis/humann/210921_SD_Meta_acute.txt", header = TRUE,check.names = FALSE, row.names=1)
GLMM_SD_acute[sapply(GLMM_SD_acute, is.character)] <- lapply(GLMM_SD_acute[sapply(GLMM_SD_acute, is.character)], 
                                       as.factor)
str(GLMM_SD_acute)
#2A. Make First Past Model
FIRSTPASSMODEL.glm <- lmer(ShannonDiversity~ (1|Patient_ID) +DOL_log + GA_birth_log+
                             ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_carbapenem+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
                             Other_ANTIM_recent+mat_ppbmi_log+mat_smoking+mat_diabetes+mat_preeclampsia+mat_chorio+mat_steroids+mat_inpt_abx_recent+
                             route+small_for_GA+hifi+pda+inf+chronic_lung_dis+intraventricular_hemorrhage+feeding_current,
                           REML=FALSE,
                           data = GLMM_SD_acute)

#2B. run Step function- initial backfitting
lmerTest::step(FIRSTPASSMODEL.glm,
     type = 3,
     alpha.random = 0.1,
     alpha.fixed = 0.05,
     reduce.fixed = TRUE,
     keep = c("Patient_ID","DOL_log"),
     fixed.calc = TRUE,
     lsmeans.calc = TRUE,
     difflsmeans.calc = TRUE)

#2C. inspect model
BACKFITTEDMODEL.glm <-lmer(ShannonDiversity~ (1|Patient_ID)+(0 + DOL_log | Patient_ID) +DOL_log+GA_birth_log+ABX_recent+Recent_3rd_gen_cephalosporin+
                                   Recent_4th_gen_cephalosporin+Other_ANTIM_recent,
                           REML=FALSE,
                           data = GLMM_SD_acute)

#pseudo-R2
r.squaredGLMM(BACKFITTEDMODEL.glm)
#summary. Look at correlation matrix. If any significant interactions with time (<-0.1 or >0.1), add to the list of variables in the first past model in part A, rerun Step. 
summary(BACKFITTEDMODEL.glm)

#2D once model has been selected, correct for multiple comparisons
g <- glht(BACKFITTEDMODEL.glm, lincfit = mcp(tension = "Tukey"))
summary(g)



####################
##### Lifetime #####
####################

GLMM_SD_lifetime<-read.delim("./Desktop/Projects/NEC/final_analysis/humann/210921_SD_Meta_lifetime.txt", header = TRUE,check.names = FALSE, row.names=1)
GLMM_SD_lifetime[sapply(GLMM_SD_lifetime, is.character)] <- lapply(GLMM_SD_lifetime[sapply(GLMM_SD_lifetime, is.character)], 
                                                             as.factor)
str(GLMM_SD_lifetime)
#2A. Make First Past Model
FIRSTPASSMODEL_lifetime.glm <- lmer(ShannonDiversity~ (1|Patient_ID) +DOL_log + GA_birth_log+
                                      mat_ppbmi_log+mat_smoking+mat_diabetes+mat_preeclampsia+mat_chorio+mat_steroids+mat_inpt_abx_recent+
                             route+small_for_GA+hifi+pda+inf+intraventricular_hemorrhage+milk_cumulative_log+formula_cumulative_log,
                           REML=FALSE,
                           data = GLMM_SD_lifetime)

#2B. run Step function- initial backfitting
lmerTest::step(FIRSTPASSMODEL_lifetime.glm,
               type = 3,
               alpha.random = 0.1,
               alpha.fixed = 0.05,
               reduce.fixed = TRUE,
               keep = c("(0+DOL_log | Patient_ID )","DOL_log"),
               fixed.calc = TRUE,
               lsmeans.calc = TRUE,
               difflsmeans.calc = TRUE)

#2C. inspect model
BACKFITTEDMODEL_lifetime.glm <-lmer(ShannonDiversity~ (1|Patient_ID)+(0 + DOL_log | Patient_ID) +DOL_log+GA_birth_log+milk_cumulative_log+ formula_cumulative_log,
                           REML=FALSE,
                           data = GLMM_SD_lifetime)

#pseudo-R2
r.squaredGLMM(BACKFITTEDMODEL_lifetime.glm)
#summary. Look at correlation matrix. If any significant interactions with time (<-0.1 or >0.1), add to the list of variables in the first past model in part A, rerun Step. 
summary(BACKFITTEDMODEL_lifetime.glm)

#2D once model has been selected, correct for multiple comparisons
g <- glht(BACKFITTEDMODEL_lifetime.glm, lincfit = mcp(tension = "Tukey"))
summary(g)

