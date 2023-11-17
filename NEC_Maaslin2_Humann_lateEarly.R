####################################
############ Maaslin2 ##############
####################################

library("Maaslin2")

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/Maaslin/210923_maaslin_metadataEarly_log.txt", header = TRUE,row.names = 1, na.strings = "NA")

SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/211006_rel_pathabundances.txt", header = TRUE,row.names = 1,check.names=FALSE)
SM_DF2<-t(SM_DF)

SM_DF_selection <- SM_DF2[row.names(SM_DF2) %in% row.names(SM_meta), ]

fit_data = Maaslin2(
  input_data = SM_DF_selection, 
  input_metadata = SM_meta,
  plot_scatter=FALSE,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  output = "./Desktop/Projects/NEC/final_analysis/humann/Maaslin/Early_log", 
  fixed_effects = c("GA_birth","NEC_status","DOL","ABX_recent","Recent_penicillin","Recent_amingoglycoside",
                    "Recent_carbapenem","Recent_macrolide","Recent_1st_gen_cephalosporin","Recent_3rd_gen_cephalosporin",
                    "Recent_4th_gen_cephalosporin","Other_ANTIM_recent","Abx_cumulative","ANTIM_cumulative","mat_preeclampsia","inf",
                    "feeding_current","milk_current","milk_cumulative","formula_cumulative"),
  random_effects = c("Patient_ID"))



##################################
######### Community PCoA #########
##################################

library(ape)
library(vegan)
library(ggplot2)
library(labdsv)
library(dplyr)
library(compositions)
library(remotes)
library(colorspace)
library("pals")

HM_LF <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/211025_rel_pathabundances_late_community.txt", header = TRUE,check.names=FALSE)
HM_matrix<-cast(HM_LF,Sample~Pathway, value='Abundance',fun.aggregate =sum)
write.table(Humann_matrix_arcsinetransformed, "Desktop/Projects/NEC/final_analysis/humann/211025_test.txt", sep="\t")

Humann_matrix_arcsinetransformed <- asin(sqrt(HM_DF))
Humann_matrix_arcsinetransformed <- as.data.frame(Humann_matrix_arcsinetransformed)

HM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/211025_rel_pathabundances_earlyCommunity.txt", header = TRUE,row.names=1)


Humann_arcsine_matrix <-vegdist(HM_DF, method = "bray")
pco_bray_humann<-pco(Humann_arcsine_matrix,k=4)
df_out_humann<-as.data.frame(pco_bray_humann$points)

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/Maaslin/210923_maaslin_metadataEarly_log.txt", header = TRUE,row.names = 1, na.strings = "NA")

df_out_humann$NEC_status<-SM_meta$NEC_status
cut<-row.names(df_out_humann)

#Calculate variance explained by axis 1 
pco_bray_humann$eig[1]/sum(pco_bray_humann$eig)
#Calculate variance explained by axis 2 
pco_bray_humann$eig[2]/sum(pco_bray_humann$eig)

ggplot(df_out_humann, aes(x=V1, y=V2, fill=as.factor(NEC_status))) +
  geom_point(size=4, shape=21) + 
  labs(fill = "NEC status") +
  scale_fill_manual(labels = c("Control", "Case"), values = c("#0D1137", "#e52165")) + 
  labs(y = "PCo2 (15.3%)", x ="PCo1 (22.5%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="bottom",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        #axis.text.x= element_text(size=12,face="bold",colour="black"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),axis.line = element_line(size=1),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y= element_blank()) +
  coord_fixed((pco_bray_humann$eig[2]/sum(pco_bray_humann$eig))/(pco_bray_humann$eig[1]/sum(pco_bray_humann$eig)))


subject_2<-read.table('Desktop/Projects/NEC/final_analysis/humann/PERMANOVA_late/Patients.txt', header = FALSE)
subject<-subject_2$V1
subject_data<-read.delim("./Desktop/Projects/NEC/final_analysis/humann/PERMANOVA_late/subject_NEC.txt", header = TRUE,check.names = FALSE, row.names=1)

PERMANOVA_repeat_measures(
  Humann_arcsine_matrix,
  subject, subject_data = subject_data,
  metadata_order = c(names(subject_data)),
  permutations=999, ncores=1)



#abx CAP with transformed data##
arcvare.cap <- capscale(HM_DF ~ NEC_status, SM_meta, dist="bray")
#vare.cap
#anova(vare.cap)
#plot(vare.cap)
arcsmry_cap <- summary(arcvare.cap)
arcCAP_for_plot <- data.frame(arcsmry_cap$sites[,1:4])
arcCAP_for_plot$Sample_ID<-row.names(arcCAP_for_plot)
arcMerged_CAP<-merge(arcCAP_for_plot, SM_meta,by.x = "row.names", by.y = "row.names")
cap_early_functions<-ggplot(arcMerged_CAP, aes(x=CAP1, y=MDS1, fill=as.factor(NEC_status))) +
  geom_point(size=5, shape=21) +
  labs(color = "NEC status") +
  scale_fill_manual(labels = c("Control", "Case"), values = c("#0D1137", "#e52165")) +
  labs(y = "", x ="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="right",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        #axis.text.x= element_text(size=12,face="bold",colour="black"),
        #axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),axis.line = element_line(size=1),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y= element_blank())
cap_early_functions
ggsave('./Desktop/Projects/NEC/final_analysis/humann/PERMANOVA_late/functions_early.png', cap_late_functions, width=8, height=4.5, bg = 'transparent')


########################################
############## Metaphlan ###############
########################################
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/PERMANOVA/PERMANOVA_late/210201_metaphlan_table_early.txt", header = TRUE,row.names = 1)

SM_DF_BC <-vegdist(SM_DF, method = "bray")

subject_2<-read.table('Desktop/Projects/NEC/final_analysis/taxonomy/PERMANOVA/PERMANOVA_late/Patients_early.txt', header = FALSE)
subject<-subject_2$V1
subject_data<-read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/PERMANOVA/PERMANOVA_late/subject_NEC_early.txt", header = TRUE,check.names = FALSE, row.names=1)

PERMANOVA_repeat_measures(
  Humann_arcsine_matrix,
  subject, subject_data = subject_data,
  metadata_order = c(names(subject_data)),
  permutations=999, ncores=1)


arcvare.cap <- capscale(SM_DF ~ NEC_status, SM_meta, dist="bray")
#vare.cap
#anova(vare.cap)
#plot(vare.cap)
arcsmry_cap <- summary(arcvare.cap)
arcCAP_for_plot <- data.frame(arcsmry_cap$sites[,1:4])
arcCAP_for_plot$Sample_ID<-row.names(arcCAP_for_plot)
arcMerged_CAP<-merge(arcCAP_for_plot, SM_meta,by.x = "row.names", by.y = "row.names")
cap_late_taxa<-ggplot(arcMerged_CAP, aes(x=CAP1, y=MDS1, fill=as.factor(NEC_status))) +
  geom_point(size=5, shape=21) +
  labs(color = "NEC status") +
  scale_fill_manual(labels = c("Control", "Case"), values = c("#0D1137", "#e52165")) +
  labs(y = "", x ="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="right",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        #axis.text.x= element_text(size=12,face="bold",colour="black"),
        #axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),axis.line = element_line(size=1),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y= element_blank())
cap_late_taxa
ggsave('./Desktop/Projects/NEC/final_analysis/taxonomy/PERMANOVA/PERMANOVA_late/PERMANOVA_taxa_early.png', cap_late_taxa, width=8, height=4.5, bg = 'transparent')



