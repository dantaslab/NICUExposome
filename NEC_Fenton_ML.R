############################
##### Machine learning #####
#####  Weight/length   #####
############################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SIAMCAT",force = TRUE)
library("curatedMetagenomicData")
library("dplyr")
library("tidyverse")
library("SIAMCAT")
library("ranger")

ML_microbiome <- read.delim("./Desktop/Projects/NEC/final_analysis/machine_learning/Combined_data_transposed_onset.txt", header = TRUE, row.names = 1,check.names = FALSE)
ML_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/machine_learning/Combined_metadata_onset.txt", header = TRUE, row.names = 1)
x<-colnames(ML_meta)
colnames_meta<-x[2:33]

####create label object
label.col <- create.label(meta=ML_meta,
                          label='NEC_status', case='1')
#####Create SIAMCAT object
sc.obj <- siamcat(feat=ML_microbiome,
                  label=label.col,
                  meta=ML_meta)

###Filter/Normalize features
sc.obj <- filter.features(sc.obj, filter.method = 'prevalence',
                          cutoff = 0.01)
sc.obj <- normalize.features(sc.obj,
                             norm.method = 'log.std',
                             norm.param=list(log.n0=1e-05,
                                             sd.min.q=0))
sc.obj <- normalize.features(sc.obj,
                             norm.method = 'pass')

############################
########## LASSO ###########
############################

sc.obj.block_25_2 <- create.data.split(sc.obj, num.folds = 10, num.resample = 10,
                                       inseparable = 'Patient')

####Full model####

#######Train LASSO with nested feature selection during cross-validation
sc.obj.block_NEC_25_enet <- train.model(sc.obj.block_25_2, method = 'enet',
                                 perform.fs = TRUE,
                                 param.fs = list(no_features = 25,
                                                 method = "AUC",
                                                 direction='absolute'))

sc.obj.block_NEC_50_enet <- train.model(sc.obj.block_25_2, method = 'enet',
                                           perform.fs = TRUE,
                                           param.fs = list(no_features = 50,
                                                           method = "AUC",
                                                           direction='absolute'))

sc.obj.block_NEC_100_enet <- train.model(sc.obj.block_25_2, method = 'enet',
                                           perform.fs = TRUE,
                                           param.fs = list(no_features = 100,
                                                           method = "AUC",
                                                           direction='absolute'))

#### Onset model #####

sc.obj.block_NEC_25_onset_2_enet <- train.model(sc.obj.block_25_2, method = 'enet',
                                     perform.fs = TRUE,
                                     param.fs = list(no_features = 25,
                                                     method = "AUC",
                                                     direction='absolute'))

sc.obj.block_NEC_50_onset_2_enet <- train.model(sc.obj.block_25_2, method = 'enet',
                                     perform.fs = TRUE,
                                     param.fs = list(no_features = 50,
                                                     method = "AUC",
                                                     direction='absolute'))

sc.obj.block_NEC_100_onset_2_enet <- train.model(sc.obj.block_25_2, method = 'enet',
                                      perform.fs = TRUE,
                                      param.fs = list(no_features = 100,
                                                      method = "AUC",
                                                      direction='absolute'))
#### 3days model #####

sc.obj.block_NEC_25_3days_2_enet <- train.model(sc.obj.block_25_2, method = 'randomForest',
                                           perform.fs = TRUE,param.fs =list(no_features = 25, method = "AUC", direction='absolute'))

sc.obj.block_NEC_50_3days_2_enet <- train.model(sc.obj.block_25_2, method = 'randomForest',
                                           perform.fs = TRUE,
                                           param.fs = list(no_features = 50,
                                                           method = "AUC",
                                                           direction='absolute'))

sc.obj.block_NEC_100_3days_2_enet <- train.model(sc.obj.block_25_2, method = 'randomForest',
                                            perform.fs = TRUE,
                                            param.fs = list(no_features = 100,
                                                            method = "AUC",
                                                            direction='absolute'))

#### 7days model #####

sc.obj.block_NEC_25_7days_2_enet <- train.model(sc.obj.block_25_2, method = 'enet',
                                           perform.fs = TRUE,
                                           param.fs = list(no_features = 25,
                                                           method = "AUC",
                                                           direction='absolute'))

sc.obj.block_NEC_50_7days_2_enet <- train.model(sc.obj.block_25_2, method = 'enet',
                                           perform.fs = TRUE,
                                           param.fs = list(no_features = 50,
                                                           method = "AUC",
                                                           direction='absolute'))

sc.obj.block_NEC_100_7days_2_enet <- train.model(sc.obj.block_25_2, method = 'enet',
                                            perform.fs = TRUE,
                                            param.fs = list(no_features = 100,
                                                            method = "AUC",
                                                            direction='absolute'))

##########Predictions
sc.obj.block_NEC_100_onset_2_enet <- make.predictions(sc.obj.block_NEC_100_onset_2_enet)
sc.obj.block_NEC_100_onset_2_enet <- evaluate.predictions(sc.obj.block_NEC_100_onset_2_enet)
model.evaluation.plot(sc.obj.block_NEC_100_enet,sc.obj.block_NEC_50_enet,sc.obj.block_NEC_25_enet,
                      colours=c('#aed6dc', '#6883bc','#1e2761'))
model.evaluation.plot(sc.obj.block_NEC_100_7days_2_enet,sc.obj.block_NEC_50_7days_2_enet,sc.obj.block_NEC_25_7days_2_enet,
                      colours=c('#aed6dc', '#6883bc','#1e2761'))
model.evaluation.plot(sc.obj.block_NEC_100_3days_2_enet,sc.obj.block_NEC_50_3days_2_enet,sc.obj.block_NEC_25_3days_2_enet,
                      colours=c('#aed6dc', '#6883bc','#1e2761'))
model.evaluation.plot(sc.obj.block_NEC_100_onset_2_enet,sc.obj.block_NEC_50_onset_2_enet,sc.obj.block_NEC_25_onset_2_enet,
                      colours=c('#aed6dc', '#6883bc','#1e2761'))

model.evaluation.plot(sc.obj.block_NEC_100_2,sc.obj.block_NEC_50_2,sc.obj.block_NEC_25_2,
                      sc.obj.block_NEC_100_RF,sc.obj.block_NEC_50_RF,sc.obj.block_NEC_25_RF,
                      colours=c('#aed6dc', '#6883bc','#1e2761',
                                '#ff9a8d','#f47a60','#ff6e40'))

model.evaluation.plot(sc.obj.block_NEC_25_2,sc.obj.block_NEC_50_2,sc.obj.block_NEC_100_2, colours=c('red', 'blue','black'))
model.interpretation.plot(sc.obj.block_30_2, consens.thres = 0.8,
                          fn.plot = './Desktop/interpret_LASSO_30.pdf')


############################
####### Elastic Net ########
############################

sc.obj.block_25_3 <- create.data.split(sc.obj, num.folds = 10, num.resample = 10,
                                       inseparable = 'Patient')
#######Train ELASTIC NET with nested feature selection during cross-validation
sc.obj.block_NEC_50_3 <- train.model(sc.obj.block_25_3, method = 'enet',
                                 perform.fs = TRUE,
                                 param.fs = list(thres.fs = 50,
                                                 method.fs = "AUC",
                                                 direction='absolute'))
##########Predictions
sc.obj.block_30_3 <- make.predictions(sc.obj.block_30_3)
sc.obj.block_30_3 <-  evaluate.predictions(sc.obj.block_30_3)
model.evaluation.plot(sc.obj.block_30_3)



############################
###### Random Forest #######
############################

######Block Cross-validation
sc.obj.block <- create.data.split(sc.obj, num.folds = 10, num.resample = 10,
                                  inseparable = 'Patient')
#######Train RF with nested feature selection during cross-validation
sc.obj.block <- train.model(sc.obj.block, method = 'randomForest', 
                            perform.fs = TRUE,
                            param.fs = list(thres.fs = 25,
                                            method.fs = "AUC",
                                            direction='absolute'))
##########Predictions
sc.obj.block <- make.predictions(sc.obj.block)
sc.obj.block <-  evaluate.predictions(sc.obj.block)
model.evaluation.plot(sc.obj.block)
model.interpretation.plot(sc.obj.block, consens.thres = 0.8,
                          fn.plot = './Desktop/interpret_nielsen.pdf')




