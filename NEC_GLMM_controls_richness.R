#####################################
####### NEC GLMM SD fitting #########
#####################################

library('lme4')
library('lmerTest')
library('MuMIn')
library('multcomp')



GLMM_SD_acute<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/tax_richness_GLMM/211116_SpeciesRichness_metadata_selection_lifetime.txt", header = TRUE,check.names = FALSE, row.names=1)
GLMM_SD_acute[sapply(GLMM_SD_acute, is.character)] <- lapply(GLMM_SD_acute[sapply(GLMM_SD_acute, is.character)], 
                                       as.factor)
str(GLMM_SD_acute)
#2A. Make First Past Model
FIRSTPASSMODEL.glm <- lmer(Richness~ (1|Patient_ID) +DOL_log + GA_birth_log+Birthweight_log+mat_age_log+Abx_cumulative_log+milk_cumulative_log+formula_cumulative_log+
                                   mat_ppbmi_log+gravida_log+para_log+enteral_first_DOL_log+payer_source+Maternal_race+mat_smoking+
                                   mat_preeclampsia+mat_chorio+mat_steroids+mat_inpt_abx_recent+mat_inpatient_meds+tracheal_intubate+surfactant+hifi+
                                   pda+inf+chronic_lung_dis+intraventricular_hemorrhage,
                           REML=FALSE,
                           data = GLMM_SD_acute)

#2B. run Step function- initial backfitting
lmerTest::step(FIRSTPASSMODEL.glm,
     type = 3,
     alpha.random = 0.1,
     alpha.fixed = 0.05,
     reduce.fixed = TRUE,
     keep = c("(0+DOL_log | Patient_ID )","DOL_log"),
     fixed.calc = TRUE,
     lsmeans.calc = TRUE,
     difflsmeans.calc = TRUE)

#2C. inspect model
BACKFITTEDMODEL.glm <-lmer(Richness~ (1|Patient_ID)+(0 + DOL_log | Patient_ID) +DOL_log+inf+mat_age_log+Abx_cumulative_log+milk_cumulative_log+
                                   formula_cumulative_log+enteral_first_DOL_log+payer_source+mat_inpatient_meds+hifi,
                           REML=FALSE,
                           data = GLMM_SD_acute)

#pseudo-R2
r.squaredGLMM(BACKFITTEDMODEL.glm)
#summary. Look at correlation matrix. If any significant interactions with time (<-0.1 or >0.1), add to the list of variables in the first past model in part A, rerun Step. 
summary(BACKFITTEDMODEL.glm)

#2D once model has been selected, correct for multiple comparisons
g <- glht(BACKFITTEDMODEL.glm, lincfit = mcp(tension = "Tukey"))
summary(g)


