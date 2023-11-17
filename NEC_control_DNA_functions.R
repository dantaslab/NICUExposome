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

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/211117_Maaslin_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/function_PERMANOVA/211128_control_human_community.txt", header = TRUE,row.names = 1)

#Calculate SD
DNA_F_SD <- as.data.frame(diversity(SM_DF, index="shannon"))
colnames(DNA_F_SD)[1] <- "ShannonDiversity"


DNA_F_SD_meta<-merge(DNA_F_SD, SM_meta,by.x = "row.names", by.y = "row.names")

#CLR transform functional data for PCA
DNA_F_clr<-clr(SM_DF)

DNA_F_clr.rda <- prcomp(DNA_F_clr)
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
pal <- brewer.pal(n = 8, name = "Blues")

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
        axis.text.x= element_blank(),axis.title.x = element_blank(), 
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
        axis.text.x= element_blank(),axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  geom_contour_fill(breaks = breaks_grad_S,
                    na.fill = TRUE)+ 
  scale_fill_gradientn(colors=pal_S,name='Shannon Div', breaks=breaks_pal_S)
gplotz2
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/Species_control_PCA.png", gplotz2, width=6, height=6, bg = "transparent")


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
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/Species_control_PCA_points.png", gplot_points, width=6, height=6, bg = "transparent")



