


library(NBZIMM)

fit1=lme.zig(Abunda_asin ~ Type+DOL_log,
             random=~1|Patient/Type, data=Table_RNA_DNA_Staph_epi, zi_fixed=~1,zi_random=NULL)

summary(fit1)


Table_RNA_DNA_Kleb_pneu<-Table_RNA_DNA[which(Table_RNA_DNA$Species=="Klebsiella_pneumoniae"),]
Table_RNA_DNA_Kleb_pneu$Abunda_asin<-asin(sqrt(Table_RNA_DNA_Kleb_pneu$Abunda_bino))

fit1=lme.zig(Abunda_asin ~ Type+log(DOL)+Type:log(DOL),
             random=~1|Patient/Type, data=Table_RNA_DNA_Kleb_pneu,
             correlation=corAR1(), zi_fixed=~1,zi_random=NULL)

summary(fit1)


Table_RNA_DNAClostr_perfri<-Table_RNA_DNA[which(Table_RNA_DNA$Species=="Clostridium_perfringens"),]
Table_RNA_DNAClostr_perfri$Abunda_asin<-asin(sqrt(Table_RNA_DNAClostr_perfri$Abunda_bino))

fit1=lme.zig(Abunda_asin ~ Type+log(DOL)+Type:log(DOL),
             random=~1|Patient/Type, data=Table_RNA_DNAClostr_perfri,
             correlation=corAR1(), zi_fixed=~1,zi_random=NULL)

summary(fit1)


library("Maaslin2")

SM_meta <- read.delim("Desktop/Projects/NEC/final_analysis/control_analysis/RNA_overrepresentation/211122_S_mitis_in.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("Desktop/Projects/NEC/final_analysis/control_analysis/RNA_overrepresentation/211122_S_mitis_data.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  output = "Desktop/Projects/NEC/final_analysis/control_analysis/RNA_overrepresentation/S_mitis_test", 
  fixed_effects = c("DOL_log","Type"),
  random_effects = c("Patient"))





library("R2admb")
library(glmmADMB)

Table_RNA <- read.table("Desktop/Projects/NEC/final_analysis/control_analysis/RNA_overrepresentation/211121_RNA_input_RNA.txt", header = TRUE, row.names = 1, check.names=FALSE)
Table_DNA <- read.table("Desktop/Projects/NEC/final_analysis/control_analysis/RNA_overrepresentation/211121_RNA_input_DNA.txt", header = TRUE, row.names = 1, check.names=FALSE)

##########################################
#######Staphylococcus epidermidis#########

Table_RNA_Staph_epi<-Table_RNA[which(Table_RNA$Species=="Staphylococcus_epidermidis"),]
Table_RNA_Staph_epi$Abunda_asin<-asin(sqrt(Table_RNA_Staph_epi$Abunda_bino))
Table_DNA_Staph_epi<-Table_DNA[which(Table_DNA$Species=="Staphylococcus_epidermidis"),]
Table_DNA_Staph_epi$Abunda_asin<-asin(sqrt(Table_DNA_Staph_epi$Abunda_bino))

Table_RNA_DNA_Staph_epi<-rbind(Table_RNA_Staph_epi,Table_DNA_Staph_epi)

Table_RNA_DNA_Staph_epi$Patient<-as.factor(Table_RNA_DNA_Staph_epi$Patient)
Table_RNA_DNA_Staph_epi$Type<-as.factor(Table_RNA_DNA_Staph_epi$Type)

fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Staph_epi, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)




#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Staph_epi, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


#######################
#######C. diff#########

Table_RNA_C_diff<-Table_RNA[which(Table_RNA$Species=="Clostridioides_difficile"),]
Table_RNA_C_diff$Abunda_asin<-asin(sqrt(Table_RNA_C_diff$Abunda_bino))
Table_DNA_C_diff<-Table_DNA[which(Table_DNA$Species=="Clostridioides_difficile"),]
Table_DNA_C_diff$Abunda_asin<-asin(sqrt(Table_DNA_C_diff$Abunda_bino))

Table_RNA_DNA_C_diff<-rbind(Table_RNA_C_diff,Table_DNA_C_diff)

Table_RNA_DNA_C_diff$Patient<-as.factor(Table_RNA_DNA_C_diff$Patient)
Table_RNA_DNA_C_diff$Type<-as.factor(Table_RNA_DNA_C_diff$Type)

fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_C_diff, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_C_diff, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


##########################################
#######Clostridium_paraputrificum#########

Table_RNA_C_paraput<-Table_RNA[which(Table_RNA$Species=="Clostridium_paraputrificum"),]
Table_RNA_C_paraput$Abunda_asin<-asin(sqrt(Table_RNA_C_paraput$Abunda_bino))
Table_DNA_C_paraput<-Table_DNA[which(Table_DNA$Species=="Clostridium_paraputrificum"),]
Table_DNA_C_paraput$Abunda_asin<-asin(sqrt(Table_DNA_C_paraput$Abunda_bino))

Table_RNA_DNA_C_paraput<-rbind(Table_RNA_C_paraput,Table_DNA_C_paraput)

Table_RNA_DNA_C_paraput$Patient<-as.factor(Table_RNA_DNA_C_paraput$Patient)
Table_RNA_DNA_C_paraput$Type<-as.factor(Table_RNA_DNA_C_paraput$Type)

fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+ (0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_C_paraput, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_C_paraput, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


##############################
#######C. perfringens#########

Table_RNA_C_perfringens<-Table_RNA[which(Table_RNA$Species=="Clostridium_perfringens"),]
Table_RNA_C_perfringens$Abunda_asin<-asin(sqrt(Table_RNA_C_perfringens$Abunda_bino))
Table_DNA_C_perfringens<-Table_DNA[which(Table_DNA$Species=="Clostridium_perfringens"),]
Table_DNA_C_perfringens$Abunda_asin<-asin(sqrt(Table_DNA_C_perfringens$Abunda_bino))

Table_RNA_DNA_C_perfringens<-rbind(Table_RNA_C_perfringens,Table_DNA_C_perfringens)

Table_RNA_DNA_C_perfringens$Patient<-as.factor(Table_RNA_DNA_C_perfringens$Patient)
Table_RNA_DNA_C_perfringens$Type<-as.factor(Table_RNA_DNA_C_perfringens$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_C_perfringens, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)

#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_C_perfringens, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


########################################
#######Clostridium_sp_7_2_43FAA#########

Table_RNA_Clostridium_sp_7_2_43FAA<-Table_RNA[which(Table_RNA$Species=="Clostridium_sp_7_2_43FAA"),]
Table_RNA_Clostridium_sp_7_2_43FAA$Abunda_asin<-asin(sqrt(Table_RNA_Clostridium_sp_7_2_43FAA$Abunda_bino))
Table_DNA_Clostridium_sp_7_2_43FAA<-Table_DNA[which(Table_DNA$Species=="Clostridium_sp_7_2_43FAA"),]
Table_DNA_Clostridium_sp_7_2_43FAA$Abunda_asin<-asin(sqrt(Table_DNA_Clostridium_sp_7_2_43FAA$Abunda_bino))

Table_RNA_DNA_Clostridium_sp_7_2_43FAA<-rbind(Table_RNA_Clostridium_sp_7_2_43FAA,Table_DNA_Clostridium_sp_7_2_43FAA)

Table_RNA_DNA_Clostridium_sp_7_2_43FAA$Patient<-as.factor(Table_RNA_DNA_Clostridium_sp_7_2_43FAA$Patient)
Table_RNA_DNA_Clostridium_sp_7_2_43FAA$Type<-as.factor(Table_RNA_DNA_Clostridium_sp_7_2_43FAA$Type)

fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+ +(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Clostridium_sp_7_2_43FAA, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Clostridium_sp_7_2_43FAA, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


###################################
#######Cutibacterium_acnes#########

Table_RNA_Cutibacterium_acnes<-Table_RNA[which(Table_RNA$Species=="Cutibacterium_acnes"),]
Table_RNA_Cutibacterium_acnes$Abunda_asin<-asin(sqrt(Table_RNA_Cutibacterium_acnes$Abunda_bino))
Table_DNA_Cutibacterium_acnes<-Table_DNA[which(Table_DNA$Species=="Cutibacterium_acnes"),]
Table_DNA_Cutibacterium_acnes$Abunda_asin<-asin(sqrt(Table_DNA_Cutibacterium_acnes$Abunda_bino))

Table_RNA_DNA_Cutibacterium_acnes<-rbind(Table_RNA_Cutibacterium_acnes,Table_DNA_Cutibacterium_acnes)

Table_RNA_DNA_Cutibacterium_acnes$Patient<-as.factor(Table_RNA_DNA_Cutibacterium_acnes$Patient)
Table_RNA_DNA_Cutibacterium_acnes$Type<-as.factor(Table_RNA_DNA_Cutibacterium_acnes$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Cutibacterium_acnes, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)

#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Cutibacterium_acnes, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


####################################
#######Cutibacterium_avidum#########

Table_RNA_Cutibacterium_avidum<-Table_RNA[which(Table_RNA$Species=="Cutibacterium_avidum"),]
Table_RNA_Cutibacterium_avidum$Abunda_asin<-asin(sqrt(Table_RNA_Cutibacterium_avidum$Abunda_bino))
Table_DNA_Cutibacterium_avidum<-Table_DNA[which(Table_DNA$Species=="Cutibacterium_avidum"),]
Table_DNA_Cutibacterium_avidum$Abunda_asin<-asin(sqrt(Table_DNA_Cutibacterium_avidum$Abunda_bino))

Table_RNA_DNA_Cutibacterium_avidum<-rbind(Table_RNA_Cutibacterium_avidum,Table_DNA_Cutibacterium_avidum)

Table_RNA_DNA_Cutibacterium_avidum$Patient<-as.factor(Table_RNA_DNA_Cutibacterium_avidum$Patient)
Table_RNA_DNA_Cutibacterium_avidum$Type<-as.factor(Table_RNA_DNA_Cutibacterium_avidum$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Cutibacterium_avidum, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)

#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Cutibacterium_avidum, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


############################################
#######Enterobacter_cloacae_complex#########

Table_RNA_Enterobacter_cloacae_complex<-Table_RNA[which(Table_RNA$Species=="Enterobacter_cloacae_complex"),]
Table_RNA_Enterobacter_cloacae_complex$Abunda_asin<-asin(sqrt(Table_RNA_Enterobacter_cloacae_complex$Abunda_bino))
Table_DNA_Enterobacter_cloacae_complex<-Table_DNA[which(Table_DNA$Species=="Enterobacter_cloacae_complex"),]
Table_DNA_Enterobacter_cloacae_complex$Abunda_asin<-asin(sqrt(Table_DNA_Enterobacter_cloacae_complex$Abunda_bino))

Table_RNA_DNA_Enterobacter_cloacae_complex<-rbind(Table_RNA_Enterobacter_cloacae_complex,Table_DNA_Enterobacter_cloacae_complex)

Table_RNA_DNA_Enterobacter_cloacae_complex$Patient<-as.factor(Table_RNA_DNA_Enterobacter_cloacae_complex$Patient)
Table_RNA_DNA_Enterobacter_cloacae_complex$Type<-as.factor(Table_RNA_DNA_Enterobacter_cloacae_complex$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Enterobacter_cloacae_complex, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Enterobacter_cloacae_complex, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


#####################################
#######Enterococcus_faecalis#########

Table_RNA_Enterococcus_faecalis<-Table_RNA[which(Table_RNA$Species=="Enterococcus_faecalis"),]
Table_RNA_Enterococcus_faecalis$Abunda_asin<-asin(sqrt(Table_RNA_Enterococcus_faecalis$Abunda_bino))
Table_DNA_Enterococcus_faecalis<-Table_DNA[which(Table_DNA$Species=="Enterococcus_faecalis"),]
Table_DNA_Enterococcus_faecalis$Abunda_asin<-asin(sqrt(Table_DNA_Enterococcus_faecalis$Abunda_bino))

Table_RNA_DNA_Enterococcus_faecalis<-rbind(Table_RNA_Enterococcus_faecalis,Table_DNA_Enterococcus_faecalis)

Table_RNA_DNA_Enterococcus_faecalis$Patient<-as.factor(Table_RNA_DNA_Enterococcus_faecalis$Patient)
Table_RNA_DNA_Enterococcus_faecalis$Type<-as.factor(Table_RNA_DNA_Enterococcus_faecalis$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Enterococcus_faecalis, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Enterococcus_faecalis, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


################################
#######Escherichia_coli#########

Table_RNA_Escherichia_coli<-Table_RNA[which(Table_RNA$Species=="Escherichia_coli"),]
Table_RNA_Escherichia_coli$Abunda_asin<-asin(sqrt(Table_RNA_Escherichia_coli$Abunda_bino))
Table_DNA_Escherichia_coli<-Table_DNA[which(Table_DNA$Species=="Escherichia_coli"),]
Table_DNA_Escherichia_coli$Abunda_asin<-asin(sqrt(Table_DNA_Escherichia_coli$Abunda_bino))

Table_RNA_DNA_Escherichia_coli<-rbind(Table_RNA_Escherichia_coli,Table_DNA_Escherichia_coli)

Table_RNA_DNA_Escherichia_coli$Patient<-as.factor(Table_RNA_DNA_Escherichia_coli$Patient)
Table_RNA_DNA_Escherichia_coli$Type<-as.factor(Table_RNA_DNA_Escherichia_coli$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Escherichia_coli, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Escherichia_coli, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")




################################
#######Finegoldia_magna#########

Table_RNA_Finegoldia_magna<-Table_RNA[which(Table_RNA$Species=="Finegoldia_magna"),]
Table_RNA_Finegoldia_magna$Abunda_asin<-asin(sqrt(Table_RNA_Finegoldia_magna$Abunda_bino))
Table_DNA_Finegoldia_magna<-Table_DNA[which(Table_DNA$Species=="Finegoldia_magna"),]
Table_DNA_Finegoldia_magna$Abunda_asin<-asin(sqrt(Table_DNA_Finegoldia_magna$Abunda_bino))

Table_RNA_DNA_Finegoldia_magna<-rbind(Table_RNA_Finegoldia_magna,Table_DNA_Finegoldia_magna)

Table_RNA_DNA_Finegoldia_magna$Patient<-as.factor(Table_RNA_DNA_Finegoldia_magna$Patient)
Table_RNA_DNA_Finegoldia_magna$Type<-as.factor(Table_RNA_DNA_Finegoldia_magna$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Finegoldia_magna, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Finegoldia_magna, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


########################################
#######Klebsiella_michiganensis#########

Table_RNA_Klebsiella_michiganensis<-Table_RNA[which(Table_RNA$Species=="Klebsiella_michiganensis"),]
Table_RNA_Klebsiella_michiganensis$Abunda_asin<-asin(sqrt(Table_RNA_Klebsiella_michiganensis$Abunda_bino))
Table_DNA_Klebsiella_michiganensis<-Table_DNA[which(Table_DNA$Species=="Klebsiella_michiganensis"),]
Table_DNA_Klebsiella_michiganensis$Abunda_asin<-asin(sqrt(Table_DNA_Klebsiella_michiganensis$Abunda_bino))

Table_RNA_DNA_Klebsiella_michiganensis<-rbind(Table_RNA_Klebsiella_michiganensis,Table_DNA_Klebsiella_michiganensis)

Table_RNA_DNA_Klebsiella_michiganensis$Patient<-as.factor(Table_RNA_DNA_Klebsiella_michiganensis$Patient)
Table_RNA_DNA_Klebsiella_michiganensis$Type<-as.factor(Table_RNA_DNA_Klebsiella_michiganensis$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Klebsiella_michiganensis, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Klebsiella_michiganensis, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


##################################
#######Klebsiella_oxytoca#########

Table_RNA_Klebsiella_oxytoca<-Table_RNA[which(Table_RNA$Species=="Klebsiella_oxytoca"),]
Table_RNA_Klebsiella_oxytoca$Abunda_asin<-asin(sqrt(Table_RNA_Klebsiella_oxytoca$Abunda_bino))
Table_DNA_Klebsiella_oxytoca<-Table_DNA[which(Table_DNA$Species=="Klebsiella_oxytoca"),]
Table_DNA_Klebsiella_oxytoca$Abunda_asin<-asin(sqrt(Table_DNA_Klebsiella_oxytoca$Abunda_bino))

Table_RNA_DNA_Klebsiella_oxytoca<-rbind(Table_RNA_Klebsiella_oxytoca,Table_DNA_Klebsiella_oxytoca)

Table_RNA_DNA_Klebsiella_oxytoca$Patient<-as.factor(Table_RNA_DNA_Klebsiella_oxytoca$Patient)
Table_RNA_DNA_Klebsiella_oxytoca$Type<-as.factor(Table_RNA_DNA_Klebsiella_oxytoca$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Klebsiella_oxytoca, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Klebsiella_oxytoca, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


#####################################
#######Klebsiella_pneumoniae#########

Table_RNA_Klebsiella_pneumoniae<-Table_RNA[which(Table_RNA$Species=="Klebsiella_pneumoniae"),]
Table_RNA_Klebsiella_pneumoniae$Abunda_asin<-asin(sqrt(Table_RNA_Klebsiella_pneumoniae$Abunda_bino))
Table_DNA_Klebsiella_pneumoniae<-Table_DNA[which(Table_DNA$Species=="Klebsiella_pneumoniae"),]
Table_DNA_Klebsiella_pneumoniae$Abunda_asin<-asin(sqrt(Table_DNA_Klebsiella_pneumoniae$Abunda_bino))

Table_RNA_DNA_Klebsiella_pneumoniae<-rbind(Table_RNA_Klebsiella_pneumoniae,Table_DNA_Klebsiella_pneumoniae)

Table_RNA_DNA_Klebsiella_pneumoniae$Patient<-as.factor(Table_RNA_DNA_Klebsiella_pneumoniae$Patient)
Table_RNA_DNA_Klebsiella_pneumoniae$Type<-as.factor(Table_RNA_DNA_Klebsiella_pneumoniae$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Klebsiella_pneumoniae, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Klebsiella_pneumoniae, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")

##########################################
#######Klebsiella_quasipneumoniae#########

Table_RNA_Klebsiella_quasipneumoniae<-Table_RNA[which(Table_RNA$Species=="Klebsiella_quasipneumoniae"),]
Table_RNA_Klebsiella_quasipneumoniae$Abunda_asin<-asin(sqrt(Table_RNA_Klebsiella_quasipneumoniae$Abunda_bino))
Table_DNA_Klebsiella_quasipneumoniae<-Table_DNA[which(Table_DNA$Species=="Klebsiella_quasipneumoniae"),]
Table_DNA_Klebsiella_quasipneumoniae$Abunda_asin<-asin(sqrt(Table_DNA_Klebsiella_quasipneumoniae$Abunda_bino))

Table_RNA_DNA_Klebsiella_quasipneumoniae<-rbind(Table_RNA_Klebsiella_quasipneumoniae,Table_DNA_Klebsiella_quasipneumoniae)

Table_RNA_DNA_Klebsiella_quasipneumoniae$Patient<-as.factor(Table_RNA_DNA_Klebsiella_quasipneumoniae$Patient)
Table_RNA_DNA_Klebsiella_quasipneumoniae$Type<-as.factor(Table_RNA_DNA_Klebsiella_quasipneumoniae$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Klebsiella_quasipneumoniae, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Klebsiella_quasipneumoniae, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


####################################
#######Klebsiella_variicola#########

Table_RNA_Klebsiella_variicola<-Table_RNA[which(Table_RNA$Species=="Klebsiella_variicola"),]
Table_RNA_Klebsiella_variicola$Abunda_asin<-asin(sqrt(Table_RNA_Klebsiella_variicola$Abunda_bino))
Table_DNA_Klebsiella_variicola<-Table_DNA[which(Table_DNA$Species=="Klebsiella_variicola"),]
Table_DNA_Klebsiella_variicola$Abunda_asin<-asin(sqrt(Table_DNA_Klebsiella_variicola$Abunda_bino))

Table_RNA_DNA_Klebsiella_variicola<-rbind(Table_RNA_Klebsiella_variicola,Table_DNA_Klebsiella_variicola)

Table_RNA_DNA_Klebsiella_variicola$Patient<-as.factor(Table_RNA_DNA_Klebsiella_variicola$Patient)
Table_RNA_DNA_Klebsiella_variicola$Type<-as.factor(Table_RNA_DNA_Klebsiella_variicola$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Klebsiella_variicola, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Klebsiella_variicola, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")



####################################
#######Klebsiella_variicola#########

Table_RNA_Klebsiella_variicola<-Table_RNA[which(Table_RNA$Species=="Klebsiella_variicola"),]
Table_RNA_Klebsiella_variicola$Abunda_asin<-asin(sqrt(Table_RNA_Klebsiella_variicola$Abunda_bino))
Table_DNA_Klebsiella_variicola<-Table_DNA[which(Table_DNA$Species=="Klebsiella_variicola"),]
Table_DNA_Klebsiella_variicola$Abunda_asin<-asin(sqrt(Table_DNA_Klebsiella_variicola$Abunda_bino))

Table_RNA_DNA_Klebsiella_variicola<-rbind(Table_RNA_Klebsiella_variicola,Table_DNA_Klebsiella_variicola)

Table_RNA_DNA_Klebsiella_variicola$Patient<-as.factor(Table_RNA_DNA_Klebsiella_variicola$Patient)
Table_RNA_DNA_Klebsiella_variicola$Type<-as.factor(Table_RNA_DNA_Klebsiella_variicola$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Klebsiella_variicola, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Klebsiella_variicola, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


#####################################
#######Staphylococcus_aureus#########

Table_RNA_Staphylococcus_aureus<-Table_RNA[which(Table_RNA$Species=="Staphylococcus_aureus"),]
Table_RNA_Staphylococcus_aureus$Abunda_asin<-asin(sqrt(Table_RNA_Staphylococcus_aureus$Abunda_bino))
Table_DNA_Staphylococcus_aureus<-Table_DNA[which(Table_DNA$Species=="Staphylococcus_aureus"),]
Table_DNA_Staphylococcus_aureus$Abunda_asin<-asin(sqrt(Table_DNA_Staphylococcus_aureus$Abunda_bino))

Table_RNA_DNA_Staphylococcus_aureus<-rbind(Table_RNA_Staphylococcus_aureus,Table_DNA_Klebsiella_variicola)

Table_RNA_DNA_Staphylococcus_aureus$Patient<-as.factor(Table_RNA_DNA_Staphylococcus_aureus$Patient)
Table_RNA_DNA_Staphylococcus_aureus$Type<-as.factor(Table_RNA_DNA_Staphylococcus_aureus$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Staphylococcus_aureus, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Staphylococcus_aureus, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


######################################
#######Staphylococcus_hominis#########

Table_RNA_Staphylococcus_hominis<-Table_RNA[which(Table_RNA$Species=="Staphylococcus_hominis"),]
Table_RNA_Staphylococcus_hominis$Abunda_asin<-asin(sqrt(Table_RNA_Staphylococcus_hominis$Abunda_bino))
Table_DNA_Staphylococcus_hominis<-Table_DNA[which(Table_DNA$Species=="Staphylococcus_hominis"),]
Table_DNA_Staphylococcus_hominis$Abunda_asin<-asin(sqrt(Table_DNA_Staphylococcus_hominis$Abunda_bino))

Table_RNA_DNA_Staphylococcus_hominis<-rbind(Table_RNA_Staphylococcus_hominis,Table_DNA_Staphylococcus_hominis)

Table_RNA_DNA_Staphylococcus_hominis$Patient<-as.factor(Table_RNA_DNA_Staphylococcus_hominis$Patient)
Table_RNA_DNA_Staphylococcus_hominis$Type<-as.factor(Table_RNA_DNA_Staphylococcus_hominis$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Staphylococcus_hominis, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Staphylococcus_hominis, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


######################################
#######Staphylococcus_warneri#########

Table_RNA_Staphylococcus_warneri<-Table_RNA[which(Table_RNA$Species=="Staphylococcus_warneri"),]
Table_RNA_Staphylococcus_warneri$Abunda_asin<-asin(sqrt(Table_RNA_Staphylococcus_warneri$Abunda_bino))
Table_DNA_Staphylococcus_warneri<-Table_DNA[which(Table_DNA$Species=="Staphylococcus_warneri"),]
Table_DNA_Staphylococcus_warneri$Abunda_asin<-asin(sqrt(Table_DNA_Staphylococcus_warneri$Abunda_bino))

Table_RNA_DNA_Staphylococcus_warneri<-rbind(Table_RNA_Staphylococcus_warneri,Table_DNA_Staphylococcus_warneri)

Table_RNA_DNA_Staphylococcus_warneri$Patient<-as.factor(Table_RNA_DNA_Staphylococcus_warneri$Patient)
Table_RNA_DNA_Staphylococcus_warneri$Type<-as.factor(Table_RNA_DNA_Staphylococcus_warneri$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Staphylococcus_warneri, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Staphylococcus_warneri, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")




###################################
#######Streptococcus_mitis#########

Table_RNA_Streptococcus_mitis<-Table_RNA[which(Table_RNA$Species=="Streptococcus_mitis"),]
Table_RNA_Streptococcus_mitis$Abunda_asin<-asin(sqrt(Table_RNA_Streptococcus_mitis$Abunda_bino))
Table_DNA_Streptococcus_mitis<-Table_DNA[which(Table_DNA$Species=="Streptococcus_mitis"),]
Table_DNA_Streptococcus_mitis$Abunda_asin<-asin(sqrt(Table_DNA_Streptococcus_mitis$Abunda_bino))

Table_RNA_DNA_Streptococcus_mitis<-rbind(Table_RNA_Streptococcus_mitis,Table_DNA_Streptococcus_mitis)

Table_RNA_DNA_Streptococcus_mitis$Patient<-as.factor(Table_RNA_DNA_Streptococcus_mitis$Patient)
Table_RNA_DNA_Streptococcus_mitis$Type<-as.factor(Table_RNA_DNA_Streptococcus_mitis$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Streptococcus_mitis, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Streptococcus_mitis, aes(x=DOL, y=log(Abundance+0.0000001),col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


##########################################
#######Streptococcus_thermophilus#########

Table_RNA_Streptococcus_thermophilus<-Table_RNA[which(Table_RNA$Species=="Streptococcus_thermophilus"),]
Table_RNA_Streptococcus_thermophilus$Abunda_asin<-asin(sqrt(Table_RNA_Streptococcus_thermophilus$Abunda_bino))
Table_DNA_Streptococcus_thermophilus<-Table_DNA[which(Table_DNA$Species=="Streptococcus_thermophilus"),]
Table_DNA_Streptococcus_thermophilus$Abunda_asin<-asin(sqrt(Table_DNA_Streptococcus_thermophilus$Abunda_bino))

Table_RNA_DNA_Streptococcus_thermophilus<-rbind(Table_RNA_Streptococcus_thermophilus,Table_DNA_Streptococcus_thermophilus)

Table_RNA_DNA_Streptococcus_thermophilus$Patient<-as.factor(Table_RNA_DNA_Streptococcus_thermophilus$Patient)
Table_RNA_DNA_Streptococcus_thermophilus$Type<-as.factor(Table_RNA_DNA_Streptococcus_thermophilus$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Streptococcus_thermophilus, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Streptococcus_thermophilus, aes(x=DOL, y=Abunda_asin,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")



##################################
#######Veillonella_dispar#########

Table_RNA_Veillonella_dispar<-Table_RNA[which(Table_RNA$Species=="Veillonella_dispar"),]
Table_RNA_Veillonella_dispar$Abunda_asin<-asin(sqrt(Table_RNA_Veillonella_dispar$Abunda_bino))
Table_DNA_Veillonella_dispar<-Table_DNA[which(Table_DNA$Species=="Veillonella_dispar"),]
Table_DNA_Veillonella_dispar$Abunda_asin<-asin(sqrt(Table_DNA_Veillonella_dispar$Abunda_bino))

Table_RNA_DNA_Veillonella_dispar<-rbind(Table_RNA_Veillonella_dispar,Table_DNA_Veillonella_dispar)

Table_RNA_DNA_Veillonella_dispar$Patient<-as.factor(Table_RNA_DNA_Veillonella_dispar$Patient)
Table_RNA_DNA_Veillonella_dispar$Type<-as.factor(Table_RNA_DNA_Veillonella_dispar$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Veillonella_dispar, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Veillonella_dispar, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


#####################################
#######Veillonella_infantium#########

Table_RNA_Veillonella_infantium<-Table_RNA[which(Table_RNA$Species=="Veillonella_infantium"),]
Table_RNA_Veillonella_infantium$Abunda_asin<-asin(sqrt(Table_RNA_Veillonella_infantium$Abunda_bino))
Table_DNA_Veillonella_infantium<-Table_DNA[which(Table_DNA$Species=="Veillonella_infantium"),]
Table_DNA_Veillonella_infantium$Abunda_asin<-asin(sqrt(Table_DNA_Veillonella_infantium$Abunda_bino))

Table_RNA_DNA_Veillonella_infantium<-rbind(Table_RNA_Veillonella_infantium,Table_DNA_Veillonella_infantium)

Table_RNA_DNA_Veillonella_infantium$Patient<-as.factor(Table_RNA_DNA_Veillonella_infantium$Patient)
Table_RNA_DNA_Veillonella_infantium$Type<-as.factor(Table_RNA_DNA_Veillonella_infantium$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Veillonella_infantium, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Veillonella_infantium, aes(x=DOL, y=Abunda_asin,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


###################################
#######Veillonella_parvula#########

Table_RNA_Veillonella_parvula<-Table_RNA[which(Table_RNA$Species=="Veillonella_parvula"),]
Table_RNA_Veillonella_parvula$Abunda_asin<-asin(sqrt(Table_RNA_Veillonella_parvula$Abunda_bino))
Table_DNA_Veillonella_parvula<-Table_DNA[which(Table_DNA$Species=="Veillonella_parvula"),]
Table_DNA_Veillonella_parvula$Abunda_asin<-asin(sqrt(Table_DNA_Veillonella_parvula$Abunda_bino))

Table_RNA_DNA_Veillonella_parvula<-rbind(Table_RNA_Veillonella_parvula,Table_DNA_Veillonella_parvula)

Table_RNA_DNA_Veillonella_parvula$Patient<-as.factor(Table_RNA_DNA_Veillonella_parvula$Patient)
Table_RNA_DNA_Veillonella_parvula$Type<-as.factor(Table_RNA_DNA_Veillonella_parvula$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Veillonella_parvula, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Veillonella_parvula, aes(x=DOL, y=Abundance,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")


#######################################
#######Veillonella_sp_T11011_6#########

Table_RNA_Veillonella_sp_T11011_6<-Table_RNA[which(Table_RNA$Species=="Veillonella_sp_T11011_6"),]
Table_RNA_Veillonella_sp_T11011_6$Abunda_asin<-asin(sqrt(Table_RNA_Veillonella_sp_T11011_6$Abunda_bino))
Table_DNA_Veillonella_sp_T11011_6<-Table_DNA[which(Table_DNA$Species=="Veillonella_sp_T11011_6"),]
Table_DNA_Veillonella_sp_T11011_6$Abunda_asin<-asin(sqrt(Table_DNA_Veillonella_sp_T11011_6$Abunda_bino))

Table_RNA_DNA_Veillonella_sp_T11011_6<-rbind(Table_RNA_Veillonella_sp_T11011_6,Table_DNA_Veillonella_sp_T11011_6)

Table_RNA_DNA_Veillonella_sp_T11011_6$Patient<-as.factor(Table_RNA_DNA_Veillonella_sp_T11011_6$Patient)
Table_RNA_DNA_Veillonella_sp_T11011_6$Type<-as.factor(Table_RNA_DNA_Veillonella_sp_T11011_6$Type)



fit_zinb <- glmmadmb(Abunda_asin ~ Type+DOL_log+(0+DOL_log|Patient/Type),
                     data=Table_RNA_DNA_Veillonella_sp_T11011_6, 
                     zeroInflation=TRUE, 
                     family="gaussian")

summary(fit_zinb)



#fit_zinb <- glmmadmb(Abunda_asin ~ Type*DOL_log+ ABX_recent+Recent_penicillin+Recent_amingoglycoside+Recent_glycopeptide+Recent_lincosamide+
#                     Recent_carbapenem+Recent_macrolide+Recent_3rd_gen_cephalosporin+Recent_4th_gen_cephalosporin+
#                     Abx_cumulative+milk_current+milk_cumulative+formula_cumulative+feeding_current+Mat_preeclampsia+(0+DOL_log|Patient/Type),
#                        data=Table_RNA_DNA_Staph_epi, 
#                        zeroInflation=TRUE, 
#                        family="nbinom")

#summary(fit_zinb)


ggplot(Table_RNA_DNA_Veillonella_parvula, aes(x=DOL, y=Abunda_asin,col=Type)) +
  geom_point() +
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold")) +
  geom_smooth( method = "loess")
