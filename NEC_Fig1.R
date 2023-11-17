####################################
############ Maaslin2 ##############
####################################

library("Maaslin2")

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/211117_Maaslin_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/function_PERMANOVA/211128_control_human_community.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  output = "./Desktop/Projects/NEC/final_analysis/control_analysis/Maaslin2_functions_0_1", 
  fixed_effects = c("DOL_log","Recent_penicillin","Recent_amingoglycoside","Recent_lincosamide","Recent_1st_gen_cephalosporin",
                    "Recent_carbapenem","Recent_3rd_gen_cephalosporin","Recent_4th_gen_cephalosporin","Other_ANTIM_recent",
                    "feeding_current","Mat_preeclampsia","Abx_cumulative","ANTIM_cumulative",
                    "milk_cumulative","formula_cumulative","mat_steroids","gravida"),
  random_effects = c("Patient_ID"))


####################################
############### PCA ################
####################################
library("compositions")
library("viridis")
library("vegan")
library("wesanderson")
library("metR")
library(RColorBrewer)
library(scico)
library(grid)

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/211117_Maaslin_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_PERMANOVA/211128_control_human_community.txt", header = TRUE,row.names = 1)

#Calculate SD
DNA_F_SD <- as.data.frame(diversity(SM_DF, index="shannon"))
colnames(DNA_F_SD)[1] <- "ShannonDiversity"


DNA_F_SD_meta<-merge(DNA_F_SD, SM_meta,by.x = "row.names", by.y = "row.names")

#CLR transform functional data for PCA
DNA_F_clr<-clr(SM_DF)

DNA_F_clr.rda <- prcomp(DNA_F_clr)
summary(DNA_F_clr.rda)

DNA_F_clr.rda_dims<-DNA_F_clr.rda$x

ordi <- ordisurf(as.matrix(DNA_F_clr.rda_dims[,1:2]), choices = c(1, 2), DNA_F_SD$ShannonDiversity,
                 col = "forestgreen",knots = 2,isotropic = FALSE,method = 'P-REML',
                 family = "gaussian",nlevels = 20) #created the ordisurf object

ordi.grid <- ordi$grid #extracts the ordisurf object
ordi.bug <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.bug$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.bug.na <- data.frame(na.omit(ordi.bug)) #gets rid of the nas

breaks_grad = seq(min(ordi.bug.na$z), max(ordi.bug.na$z), by = (max(ordi.bug.na$z)-min(ordi.bug.na$z))/10)
breaks_pal = seq(min(ordi.bug.na$z), max(ordi.bug.na$z), by = (max(ordi.bug.na$z)-min(ordi.bug.na$z))/10)
pal <- brewer.pal(n = 8, name = "Greys")

#VERSION 1
#gplotz <- ggplot(data = ordi.bug.na,
#                 aes(x = x,
#                     y = y,
#                     z = z)) +
#  theme_bw()+
#  theme(legend.position="none",
#        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=0.5), 
#        axis.text.x= element_blank(),axis.title.x = element_blank(), 
#        axis.title.y = element_blank(),
#        axis.text.y= element_blank())+
#  geom_contour_filled(breaks = breaks_grad,
#                    na.fill = TRUE,color="black")+
#  xlim(-15,15)+
#  ylim(-25,10)

#gplotz
#ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/Functions_control_PCA_ALT.png", gplotz, width=6, height=6, bg = "transparent")

#VERSION 2
gplotz2 <- ggplot(data = ordi.bug.na,
                 aes(x = x,
                     y = y,
                     z = z)) +
  theme_bw()+
  theme(legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=0.5), 
        axis.text.x= element_blank(),axis.title.x = element_blank(), axis.ticks=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  xlim(-15,15)+
  ylim(-25,10)+
  geom_contour_fill(breaks = breaks_grad,
                    na.fill = TRUE)+ 
  scale_fill_gradientn(colors=pal,name='Shannon Div', breaks=breaks_pal)
gplotz2
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/Functions_control_PCA.png", gplotz2, width=6, height=6, bg = "transparent")


DNA_F_clr.rda_dims_DF<-as.data.frame(DNA_F_clr.rda_dims)
DNA_F_clr.rda_dims_DF_meta<-merge(DNA_F_clr.rda_dims_DF, SM_meta,by.x = "row.names", by.y = "row.names")


gplot_points<-ggplot(DNA_F_clr.rda_dims_DF_meta, aes(x=PC1, y=PC2,fill=DOL_log)) +
  geom_point(size=2.5,pch=21,stroke=1)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),panel.background = element_rect(fill = "transparent"),
        axis.text.x= element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  ylab('PCO2')+
  xlab('PCO1')+
  xlim(-15,15)+
  ylim(-25,10)+ scale_fill_viridis_c(option = "magma")
gplot_points
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/Functions_control_PCA_points.png", gplot_points, width=6, height=6, bg = "transparent")



####################################
############ Taxonomy ##############
####################################

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/211117_Maaslin_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF_species <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/211117_metaphlan_table_controls.txt", header = TRUE,row.names = 1)


#Calculate SD
DNA_S_SD <- as.data.frame(diversity(SM_DF_species, index="shannon"))
colnames(DNA_S_SD)[1] <- "ShannonDiversity"


DNA_S_SD_meta<-merge(DNA_S_SD, SM_meta,by.x = "row.names", by.y = "row.names")

#CLR transform functional data for PCA
DNA_S_clr<-clr(SM_DF_species)

DNA_S_clr.rda <- prcomp(DNA_S_clr)
summary(DNA_S_clr.rda)

DNA_S_clr.rda_dims<-DNA_S_clr.rda$x

ordi_S <- ordisurf(as.matrix(DNA_S_clr.rda_dims[,1:2]), choices = c(1, 2), DNA_S_SD$ShannonDiversity,
                 col = "forestgreen",knots = 2,isotropic = FALSE,method = 'P-REML',
                 family = "gaussian",nlevels = 20) #created the ordisurf object

ordi_S.grid <- ordi_S$grid #extracts the ordisurf object
ordi_S.bug <- expand.grid(x = ordi_S.grid$x, y = ordi_S.grid$y) #get x and ys
ordi_S.bug$z <- as.vector(ordi_S.grid$z) #unravel the matrix for the z scores
ordi_S.bug.na <- data.frame(na.omit(ordi_S.bug)) #gets rid of the nas

breaks_grad_S = seq(min(ordi_S.bug.na$z), max(ordi_S.bug.na$z), by = (max(ordi_S.bug.na$z)-min(ordi_S.bug.na$z))/10)
breaks_pal_S = seq(min(ordi_S.bug.na$z), max(ordi_S.bug.na$z), by = (max(ordi_S.bug.na$z)-min(ordi_S.bug.na$z))/10)
pal_S <- brewer.pal(n = 8, name = "Greys")

#VERSION 2
gplotz2 <- ggplot(data = ordi_S.bug.na,
                  aes(x = x,
                      y = y,
                      z = z)) +
  theme_bw()+
  theme(legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=0.5), 
        axis.text.x= element_blank(),axis.title.x = element_blank(), axis.ticks=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  geom_contour_fill(breaks = breaks_grad_S,
                    na.fill = TRUE)+ 
  scale_fill_gradientn(colors=pal_S,name='Shannon Div', breaks=breaks_pal_S)+
  xlim(-8,8)+
  ylim(-6,6)
gplotz2
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/Species_control_PCA.png", gplotz2, width=6, height=6, bg = "transparent")


DNA_S_clr.rda_dims_DF<-as.data.frame(DNA_S_clr.rda_dims)
DNA_S_clr.rda_dims_DF_meta<-merge(DNA_S_clr.rda_dims_DF, SM_meta,by.x = "row.names", by.y = "row.names")


gplot_points<-ggplot(DNA_S_clr.rda_dims_DF_meta, aes(x=PC1, y=PC2,fill=DOL_log)) +
  geom_point(size=3,pch=21,stroke=1)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),panel.background = element_rect(fill = "transparent"),
        axis.text.x= element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  ylab('PCO2')+
  xlab('PCO1')+
  xlim(-8,8)+
  ylim(-6,6)+ 
  scale_fill_viridis_c(option = "magma")
gplot_points
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/Species_control_PCA_points.png", gplot_points, width=6, height=6, bg = "transparent")

#############################
######## Procrustes #########
#############################

pro <- procrustes(X = DNA_S_clr.rda, Y = DNA_F_clr.rda,symmetric=TRUE)
pro


protest(X = DNA_S_clr.rda, Y = DNA_F_clr.rda, scores = "sites", permutations = 999)

ctest <- data.frame(rda1=pro$Yrot[,1],
                    rda2=pro$Yrot[,2],xrda1=pro$X[,1],
                    xrda2=pro$X[,2])

ctest_DF<-as.data.frame(ctest)
ctest_DF_meta<-merge(ctest_DF, SM_meta,by.x = "row.names", by.y = "row.names")

procrustes<-ggplot(ctest_DF_meta) +
  geom_point(aes(x=rda1, y=rda2),pch=21,size=3,fill="#ecc19c",) +
  geom_point(aes(x=xrda1, y=xrda2),pch=21,size=3,fill="black",alpha=0.3) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2),lwd=0.05)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black",size=0.5, fill = NA),
        axis.text.x= element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())
procrustes
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/FigS2A.png", procrustes, width=8, height=8, bg = "transparent")

##################################
######## DOL vs GA_DOL ###########
##################################
library(ggplot2)
library(reshape)
library(vegan)

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/metadata/DOL_GADOL.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF_species <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/211117_metaphlan_table_controls.txt", header = TRUE,row.names = 1)
species_DF_BC <-vegdist(SM_DF_species, method = 'bray')

species_BC_matrix <- as.matrix(species_DF_BC)
#Consecutive samples
species_BC_longform <- melt(species_BC_matrix, varnames = c("ID1", "ID2"))
SM_DF_matrix_meta1<-merge(species_BC_longform, SM_meta,by.x = "ID1", by.y = "row.names")
colnames(SM_DF_matrix_meta1)[4:ncol(SM_DF_matrix_meta1)] <- paste("1", colnames(SM_meta)[1:ncol(SM_meta)], sep = "_")
SM_DF_matrix_meta2<-merge(SM_DF_matrix_meta1, SM_meta,by.x = "ID2", by.y = "row.names")
colnames(SM_DF_matrix_meta2)[10:ncol(SM_DF_matrix_meta2)] <- paste("2", colnames(SM_meta)[1:ncol(SM_meta)], sep = "_")

#subset to consecutive samples
SM_BC_C<-subset(SM_DF_matrix_meta2,SM_DF_matrix_meta2$`1_Patient`!=SM_DF_matrix_meta2$`2_Patient`)
SM_BC_C$DOLdiff<-abs(SM_BC_C$`1_DOL`-SM_BC_C$`2_DOL`)
SM_BC_C$GADOLdiff<-abs(SM_BC_C$`1_GA_DOL`-SM_BC_C$`2_GA_DOL`)

SM_BC_C_DOL_1_5<-subset(SM_BC_C,SM_BC_C$`1_DOL`<=5&SM_BC_C$`2_DOL`<=5)
SM_BC_C_DOL_3_7<-subset(SM_BC_C,SM_BC_C$`1_DOL`<=7&SM_BC_C$`2_DOL`<=7&SM_BC_C$`1_DOL`>=3&SM_BC_C$`2_DOL`>=3)

SM_BC_C_DOL_7_above<-subset(SM_BC_C,SM_BC_C$`1_DOL`>7&SM_BC_C$`2_DOL`>7)
SM_BC_C_DOL_14<-subset(SM_BC_C_DOL_7_above,SM_BC_C_DOL_7_above$`1_DOL`<=14&SM_BC_C_DOL_7_above$`2_DOL`<=14)
SM_BC_C_DOL_14_above<-subset(SM_BC_C_DOL_7_above,SM_BC_C_DOL_7_above$`1_DOL`>14&SM_BC_C_DOL_7_above$`2_DOL`>14)
SM_BC_C_DOL_21<-subset(SM_BC_C_DOL_14_above,SM_BC_C_DOL_14_above$`1_DOL`<=21&SM_BC_C_DOL_14_above$`2_DOL`<=21)
SM_BC_C_DOL_21_above<-subset(SM_BC_C_DOL_14_above,SM_BC_C_DOL_14_above$`1_DOL`>21&SM_BC_C_DOL_14_above$`2_DOL`>21)
SM_BC_C_DOL_28<-subset(SM_BC_C_DOL_21_above,SM_BC_C_DOL_21_above$`1_DOL`<=28&SM_BC_C_DOL_21_above$`2_DOL`<=28)
SM_BC_C_DOL_28_above<-subset(SM_BC_C_DOL_21_above,SM_BC_C_DOL_21_above$`1_DOL`>28&SM_BC_C_DOL_21_above$`2_DOL`>28)
SM_BC_C_DOL_35<-subset(SM_BC_C_DOL_28_above,SM_BC_C_DOL_28_above$`1_DOL`<=35&SM_BC_C_DOL_28_above$`2_DOL`<=35)
SM_BC_C_DOL_35_above<-subset(SM_BC_C_DOL_28_above,SM_BC_C_DOL_28_above$`1_DOL`>35&SM_BC_C_DOL_28_above$`2_DOL`>35)
SM_BC_C_DOL_42<-subset(SM_BC_C_DOL_35_above,SM_BC_C_DOL_35_above$`1_DOL`<=42&SM_BC_C_DOL_35_above$`2_DOL`<=42)
SM_BC_C_DOL_42_above<-subset(SM_BC_C_DOL_35_above,SM_BC_C_DOL_35_above$`1_DOL`>42&SM_BC_C_DOL_35_above$`2_DOL`>42)
SM_BC_C_DOL_49<-subset(SM_BC_C_DOL_42_above,SM_BC_C_DOL_42_above$`1_DOL`<=49&SM_BC_C_DOL_42_above$`2_DOL`<=49)
SM_BC_C_DOL_49_above<-subset(SM_BC_C_DOL_42_above,SM_BC_C_DOL_42_above$`1_DOL`>49&SM_BC_C_DOL_42_above$`2_DOL`>49)
SM_BC_C_DOL_56<-subset(SM_BC_C_DOL_49_above,SM_BC_C_DOL_49_above$`1_DOL`<=56&SM_BC_C_DOL_49_above$`2_DOL`<=56)



SM_BC_C<-SM_BC_C[sample(nrow(SM_BC_C), 35000),]

BC<-ggplot(SM_BC_C_DOL_3_7, aes(x=GADOLdiff, y=value)) +
  geom_point(pch='.',alpha=1)+geom_smooth(method="lm",col='#253494')+theme_bw()+xlim(0,80)+
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.text.x= element_blank(),axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits = c(0,1.1))+
  ylab("BC distance")+
  xlab("∆GA adjusted DOL")
BC  
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/GA/DOL_56_BC.png", BC, width=6, height=5, bg = "transparent")


#DOL
BC<-ggplot(SM_BC_C, aes(x=DOLdiff, y=value)) +
  geom_point(pch='.',alpha=0)+geom_smooth(method="loess",col='black')+theme_bw()+xlim(0,80)+
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.text.x= element_blank(),axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits = c(0,1.1))+
  ylab("BC distance")+
  xlab("∆DOL")
BC  
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/GA/DOL_BC.png", BC, width=6, height=5, bg = "transparent")






SM_BC_C_DOL_20_within<-subset(SM_BC_C_DOL_20,SM_BC_C_DOL_20$`1_Hospital`==SM_BC_C_DOL_20$`2_Hospital`)
SM_BC_C_DOL_20_within<-SM_BC_C_DOL_20_within[sample(nrow(SM_BC_C_DOL_20_within), 35000),]
SM_BC_C_DOL_20_between<-subset(SM_BC_C_DOL_20,SM_BC_C_DOL_20$`1_Hospital`!=SM_BC_C_DOL_20$`2_Hospital`)
SM_BC_C_DOL_20_between<-SM_BC_C_DOL_20_between[sample(nrow(SM_BC_C_DOL_20_between), 35000),]

BC<-ggplot(SM_BC_C_DOL_20_between, aes(x=GADOLdiff, y=value)) +
  geom_point(pch='.',alpha=0)+geom_smooth(method="loess",col='red')+ylim(0,1)+theme_bw()+xlim(0,80)+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),panel.background = element_rect(fill = "transparent"),
        axis.text.x= element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  ylab("BC distance")+
  xlab("∆GA adjusted DOL")
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/GA/DOL_20_BC_between.png", BC, width=6, height=5, bg = "transparent")


SM_BC_C_DOL_40_within<-subset(SM_BC_C_DOL_40,SM_BC_C_DOL_40$`1_Hospital`==SM_BC_C_DOL_40$`2_Hospital`)
#SM_BC_C_DOL_40_within<-SM_BC_C_DOL_40_within[sample(nrow(SM_BC_C_DOL_40_within), 35000),]
SM_BC_C_DOL_40_between<-subset(SM_BC_C_DOL_40,SM_BC_C_DOL_40$`1_Hospital`!=SM_BC_C_DOL_40$`2_Hospital`)
#SM_BC_C_DOL_40_between<-SM_BC_C_DOL_40_between[sample(nrow(SM_BC_C_DOL_40_between), 35000),]

BC<-ggplot(SM_BC_C_DOL_40_within, aes(x=GADOLdiff, y=value)) +
  geom_point(pch='.',alpha=0)+geom_smooth(method="loess",col='purple')+ylim(0,1)+theme_bw()+xlim(0,80)+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),panel.background = element_rect(fill = "transparent"),
        axis.text.x= element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  ylab("BC distance")+
  xlab("∆GA adjusted DOL")
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/GA/DOL_40_BC_within.png", BC, width=6, height=5, bg = "transparent")


##################################
####### Weight/LengthGain ########
##################################

WG_week <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_weight_gain/old/weight_gain_by_week_ALL_WEEKS.txt", header = TRUE)


w_WOL<-ggplot(WG_week, aes(x = Week, y = WG_per)) + 
  geom_line(aes(color = Patient),alpha=0.2)+
  geom_smooth(method="loess",col="black")+
  geom_hline(yintercept=0)+
  theme_classic()+
  ylab("∆weight to previous week")+
  xlab("Week of life")+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black",size=0.5, fill = NA))
w_WOL
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/∆weight_WOL.png", w_WOL, width=6, height=4, bg = "transparent")



w_hist<-ggplot(WG_week, aes(x=WG_per))+
  geom_histogram(color="black", fill="white",binwidth=2)+
  facet_grid(cols = vars(Week), scales = "free_x")+
  xlim(-40,40)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black",size=0.5, fill = NA))
w_hist
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/weightHist_WOL.png", w_hist, width=10, height=1.5, bg = "transparent")



by(WG_week, WG_week$Week, summary)

Length_gain <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/length_gain/length_gain.txt", header = TRUE)
L_WOL<-ggplot(Length_gain, aes(x = week, y = lngth_gain)) + 
  geom_line(aes(color = study_id),alpha=0.2)+
  geom_hline(yintercept=0)+
  geom_smooth(method="loess",col="black")+
  theme_classic()+
  ylab("∆length to previous week")+
  xlab("Week of life")+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black",size=0.5, fill = NA))
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/∆length_WOL.png", L_WOL, width=6, height=4, bg = "transparent")


l_hist<-ggplot(Length_gain, aes(x=lngth_gain))+
  geom_histogram(color="black", fill="white",binwidth=2)+
  facet_grid(cols = vars(week), scales = "free_x")+
  xlim(-20,20)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black",size=0.5, fill = NA))
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/lengthHist_WOL.png", l_hist, width=10, height=1.5, bg = "transparent")


WG_contin_week <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/weight_gain/subject_continuous.txt", header = TRUE, na.strings = "NA")


#GLMMM
names(diversity_LA_plot)[2]<- 'SD'
diversity_LA_plot$Cohort<-as.factor(diversity_LA_plot$Cohort)

model.b = lme(SD ~ Cohort*DOL, 
              random = ~1|Patient, data=diversity_LA_plot)
summary(model.b)

model.null = lme(SD ~ DOL, 
                 random = ~1|Patient, data=diversity_LA_plot)

anova(update(model.b, . ~ ., method = 'ML'),
      update(model.null, . ~ ., method = 'ML'))

library (multcomp)
summary(glht(model.b, mcp(Cohort='Tukey')))




IDs<-rev(unique(WG_contin_week$Week))
results <- list()
for(i in 1:10){
  temp<-WG_contin_week[WG_contin_week$Week==IDs[i],]
  temp$Quartile <- factor(temp$Quartile, levels=c("1","2","3","4"))
  temp_formulae <- lapply(colnames(temp)[4:ncol(temp)], function(x) as.formula(paste0(x, " ~ Quartile")))
  temp_res <- lapply(temp_formulae, function(x) summary(aov(x, data = temp)))
  names(temp_res) <- format(temp_formulae)
  results[[i]] <- temp_res
}

new_results <- plyr::adply(results,1,unlist,.id = NULL)
write.csv(new_results,"./Desktop/Projects/NEC/final_analysis/control_analysis/weight_gain/subject_continuous_results.csv")

R<-ggplot(WG_contin_week_3, aes(x=WG_contin_week_3$Quartile, y=WG_contin_week_3$Birthweight)) +
  theme_bw()+
  geom_boxplot()
R

R<-ggplot(WG_contin_week_3, aes(x=WG_contin_week_3$Quartile, y=WG_contin_week_3$GA_birth)) +
  theme_bw()+
  geom_boxplot()
R

R<-ggplot(WG_contin_week_5, aes(x=WG_contin_week_5$Quartile, y=WG_contin_week_5$GA_birth)) +
  theme_bw()+
  geom_boxplot()
R

R<-ggplot(WG_contin_week_6, aes(x=WG_contin_week_6$Quartile, y=WG_contin_week_6$mat_age)) +
  theme_bw()+
  geom_boxplot()
R



stool_col_type_s_GO_ORA <- read.delim('Desktop/Projects/MDRO_UTI_COHORT/Paper_split/Paper2/Revision CHM/MGE_blast_results/urine_species.txt', check.names = FALSE, header = TRUE, row.names = 1)
s_col_s<-as.data.frame(apply(stool_col_type_s_GO_ORA, 1, 
                             function(x) {
                               tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                               fisher.test(tbl, alternative="less")$p.value
                             }))
s_col_s$p.adjust<-sapply(s_col_s,p.adjust,method="fdr")
write.table(s_col_s, "Desktop/Projects/MDRO_UTI_COHORT/Paper_split/Paper2/Revision CHM/MGE_blast_results/urine_underrepresented_taxa_stats.txt", sep="\t")




##################################
####### DNA/RNA diversity ########
##################################

#####Species DNA
SM_DF_DNA <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Metaphlan/211117_DNA_metaphlan_control_matrix.txt", header = TRUE,row.names = 1)

#Calculate SD
DNA_Tax_SD <- as.data.frame(diversity(SM_DF_DNA, index="shannon"))
colnames(DNA_Tax_SD)[1] <- "ShannonDiversity_DNA"


#####Species RNA
SM_DF_RNA <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Metaphlan/211108_RNA_metaphlan_control_matrix.txt", header = TRUE,row.names = 1)

#Calculate SD
RNA_Tax_SD <- as.data.frame(diversity(SM_DF_RNA, index="shannon"))
colnames(RNA_Tax_SD)[1] <- "ShannonDiversity_RNA"


SD_RNA_DNA<-merge(DNA_Tax_SD, RNA_Tax_SD,by.x = "row.names", by.y = "row.names")



RNA_DNA_SD_plot<-ggplot(SD_RNA_DNA, aes(x=ShannonDiversity_DNA, y=ShannonDiversity_RNA)) +
  geom_point(pch='.')+geom_smooth(method="lm",col='black')+theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))+
  ylab("RNA SD")+
  xlab("DNA SD")
RNA_DNA_SD_plot

ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/GA/DOL_40_BC_within.png", BC, width=6, height=5, bg = "transparent")


##########################
#### Species dyanmics ####
##########################
library(scales)
Taxa_dynamics <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_taxa_dynamics/FigS1_C_F.txt", header = TRUE)


p <-ggplot(Taxa_dynamics, aes(x=DOL, y=Prevalence, col=Genus, linetype=as.factor(line_type))) +
  geom_point(pch='.',alpha=0)+geom_smooth(method="loess")+theme_bw()+
  ylab("")+
  xlab("Day of life")+
  theme_classic()+
  theme(legend.text=element_text(size=12, face='bold'),legend.title=element_blank(),legend.position = "top",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x= element_text(size=14, face='bold',colour='black'),
        axis.title.x = element_text(size=18, vjust = 1,face='bold'),
        axis.title.y = element_text(size=18, hjust = 1,face='bold'),
        axis.text.y= element_text(size=14,face='bold',colour='black'),
        strip.text = element_text(size=14,face="bold"))+
  scale_y_continuous(limit=c(0,1),breaks=c(0,0.5,1),oob=squish)+ facet_wrap(~Genus, ncol = 3)+scale_x_continuous(breaks=c(0,20,40,60,80),limits=c(0,80))+
  scale_color_manual(values=c("#ffc13b","#ecc19c","#aed6dc","#81b7d2","#4a536b","#77c593","#B85042","#ff6e40","#d13ca4"))
p  
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_taxa_dynamics/Suppl_figure_species.png", p, width=10, height=8, bg = "transparent")



