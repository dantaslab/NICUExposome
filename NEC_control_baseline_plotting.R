library("compositions")
library("viridis")
library("vegan")
library("wesanderson")
library("metR")
library(RColorBrewer)
library(scico)
library(grid)

SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/metaphlan_community_STL.txt", header = TRUE,check.names=FALSE,row.names = 1)

metadata <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/metaphlan_metadata_STL.txt", header = TRUE,check.names=FALSE,row.names = 1)

SM_DF_taxa_clr<-clr(SM_DF)



####################################
####### Composition over DOL #######
####################################
library("compositions")
library("viridis")
library("vegan")
library("wesanderson")
library("metR")
library(RColorBrewer)
library(scico)
library(grid)

#Calculate SD
Twins_taxa_SD <- as.data.frame(diversity(SM_DF, index="shannon"))
colnames(Twins_taxa_SD)[1] <- "ShannonDiversity"

Twins_taxa_SD_meta<-merge(Twins_taxa_SD, metadata,by.x = "row.names", by.y = "row.names")

#CLR transform functional data for PCA
Twins_taxa_clr<-clr(SM_DF)

Twins_taxa_clr.rda <- prcomp(Twins_taxa_clr)
summary(Twins_taxa_clr.rda)

Twins_taxa_clr.rda_dims<-Twins_taxa_clr.rda$x

ordi <- ordisurf(as.matrix(Twins_taxa_clr.rda_dims[,1:2]), choices = c(1, 2), Twins_taxa_SD$ShannonDiversity,
                 col = "forestgreen",knots = 2,isotropic = FALSE,method = 'P-REML',
                 family = "gaussian",nlevels = 20) #created the ordisurf object

ordi.grid <- ordi$grid #extracts the ordisurf object
ordi.bug <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.bug$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.bug.na <- data.frame(na.omit(ordi.bug)) #gets rid of the nas

breaks_grad = c(-0.2,0,0.2,0.4,0.6,0.8,1,1.2)
  seq(0, max(ordi.bug.na$z), by = (max(ordi.bug.na$z)-0)/6)
breaks_pal = c(-0.2,0,0.2,0.4,0.6,0.8,1,1.2)
  seq(0, max(ordi.bug.na$z), by = (max(ordi.bug.na$z)-0)/6)
pal <- brewer.pal(n = 8, name = "Greens")
colfunc <- colorRampPalette(c("#a9c0a6","#1d3c45"))
pal<-colfunc(6)

#VERSION 2
gplotz2 <- ggplot(data = ordi.bug.na,
                  aes(x = x,
                      y = y,
                      z = z)) +
  theme_bw()+
  theme(legend.position="none",legend.key.width = unit(2, 'cm'),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=0.5), 
        axis.text.x= element_blank(),axis.title.x = element_blank(), axis.ticks=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  geom_contour_fill(breaks = breaks_grad,
                    na.fill = TRUE)+ 
  xlim(-7.5,7.5)+
  ylim(-6,6)+
  scale_fill_gradientn(colors=pal,name='Shannon Div', breaks=breaks_pal)
gplotz2
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/PCoA_taxa_background.png", gplotz2, width=4.2, height=4.2, bg = "transparent")


Twins_taxa_clr.rda_dims_DF<-as.data.frame(Twins_taxa_clr.rda_dims)
Twins_taxa_clr.rda_dims_DF_meta<-merge(Twins_taxa_clr.rda_dims_DF, metadata,by.x = "row.names", by.y = "row.names")
Twins_taxa_clr.rda_dims_DF_meta

colfunc <- colorRampPalette(c("#fff1e1","#d2601a"))
colfunc(6)

gplot_points_infants<-ggplot(Twins_taxa_clr.rda_dims_DF_meta, aes(x=PC1, y=PC2,fill=Group)) +
  geom_point(size=2, pch=21,alpha=0.2)+
  geom_point(data=Twins_taxa_clr.rda_dims_DF_meta %>% 
               group_by(Group) %>% 
               summarise_at(vars(matches("PC")), mean),
             size=5, shape=21)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = "transparent"),
        axis.text.x= element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  xlim(-7.5,7.5)+
  ylim(-6,6)+
  ylab('PCO2')+
  xlab('PCO1')+
  scale_fill_manual(values = c("#FFF1E1","#F6D3B9","#EDB791","#E39969","#DB7D41","#D2601A"))+
  scale_shape_manual(values=c(21))
gplot_points_infants
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/PCoA_taxa_points.png", gplot_points_infants, width=4, height=4, bg = "transparent")



##################################
####### Genus introduction #######
##################################

genus_abundance <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_NEC_control_PA_community_genus_selection.txt", header = TRUE)
metadata <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_NEC_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")



genus_abundance_melt <- melt(genus_abundance, varnames = c("ID1", "ID2"))
genus_abundance_meta<-merge(genus_abundance_melt, metadata,by.x = "X", by.y = "row.names")

p <-ggplot(genus_abundance_meta, aes(x=DOL, y=value, col=variable)) +
  geom_point(alpha=0)+geom_smooth(method="loess",alpha=0.25)+theme_bw()+
  ylab("Genus present [%]")+
  xlab("DOL")+
  theme_classic()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line =  element_line(size=0.25),
        panel.background = element_blank())+
  scale_color_manual(values = c("#ffc13b","#ecc19c","#aed6dc","#81b7d2","#4a536b","#77c593","#B85042","#ff6e40","#d13ca4"))+
  scale_x_continuous(breaks=c(0,20,40,60),limits=c(0,60))+
  scale_y_continuous(limit=c(-0.2,1.2))+
  coord_cartesian(ylim = c(0,1))
p  
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/Genus_progression.png", p, width=6, height=3, bg = "transparent")


##################################
######## Genus abundance #########
##################################

genus_abundance <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/Genus_abundance.txt", header = TRUE)
metadata <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_NEC_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

genus_abundance_melt <- melt(genus_abundance, varnames = c("ID1", "ID2"))
genus_abundance_meta<-merge(genus_abundance_melt, metadata,by.x = "SampleID", by.y = "row.names")

write.csv(genus_abundance_meta,"./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/DOL_abundance.csv")

genus_abundance <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/DOL_abundance_average.txt", header = TRUE)
genus_abundance$variable <- factor(genus_abundance$variable, levels=c("Bifidobacterium","Clostridium","Enterobacter","Enterococcus",
                                                                      "Escherichia","Klebsiella","Staphylococcus","Streptococcus","Veillonella","Others"))


barplot<-ggplot(genus_abundance, aes(fill=variable, y=value, x=DOL)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("#d9a5b3","#ffc13b","#aed6dc","#81b7d2","#4a536b","#77c593","#B85042","#ff6e40","#d13ca4","#818181"))+
  theme_classic()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),axis.line.y =  element_line(size=0.5),axis.line.x =  element_line(size=0.0),
        panel.background = element_blank())+
  coord_cartesian(ylim = c(0,100))+
  scale_y_continuous(expand = c(0.01,0.01))+
  scale_x_continuous(expand = c(0.01,0.01))
barplot  
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/Relative_abundance_DOL.png", barplot, width=6, height=3, bg = "transparent")





  
  geom_point(alpha=0)+geom_smooth(method="loess",alpha=0.25)+theme_bw()+
  ylab("Genus present [%]")+
  xlab("DOL")+
  theme_classic()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line =  element_line(size=0.25),
        panel.background = element_blank())+
  scale_color_manual(values = c("#ffc13b","#ecc19c","#aed6dc","#81b7d2","#4a536b","#77c593","#B85042","#ff6e40","#d13ca4"))+
  scale_x_continuous(breaks=c(0,20,40,60),limits=c(0,60))+
  scale_y_continuous(limit=c(-0.2,1.2))+
  coord_cartesian(ylim = c(0,1))
p  




#Calculate SD
SD_time <- as.data.frame(diversity(SM_DF, index="shannon"))
colnames(SD_time)[1] <- "ShannonDiversity"

SD_time_meta<-merge(SD_time, metadata,by.x = "row.names", by.y = "row.names")

BC<-ggplot(SD_time_meta, aes(x=DOL, y=ShannonDiversity)) +
  geom_point(size=2, pch=21,alpha=0.2)+geom_smooth(method="loess",col="red")+theme_bw()+xlim(0,80)+
  theme_classic()+
  ylim(0,2.5)+
  theme(legend.text=element_text(size=14, face='bold'),legend.title=element_text(size=14, face='bold'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x= element_text(size=14,angle=45,hjust=1,vjust=1,face='bold',colour='black'),
        axis.title.x = element_text(size=18, face='bold'),
        axis.title.y = element_text(size=18, face='bold'),
        axis.text.y= element_text(size=14,face='bold',colour='black'))
BC


SM_DF_richness<-as.data.frame(rowSums(SM_DF[,1:ncol(SM_DF)]>0.0))
colnames(SM_DF_richness)[1] <-"Richness"
SM_DF_richness_meta<-merge(SM_DF_richness, metadata,by.x = 'row.names', by.y = "row.names")

rare_R_SB<-ggplot(SM_DF_richness_meta, aes(x=DOL, y=Richness))+
  geom_point(size=2, pch=21,alpha=0)+
  theme_classic()+
  theme(legend.text=element_text(size=14, face='bold'),legend.title=element_text(size=14, face='bold'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  ylab('Richness')+
  xlab('DOL')+geom_smooth(method="loess",col="red")+
  scale_y_continuous(limits=c(0,30),expand = c(0.01,0.01))+
  scale_x_continuous(limits=c(0,80),expand = c(0.01,0.01))
rare_R_SB
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/Richness.png", rare_R_SB, width=6, height=2, bg = "transparent")


library(tidyverse)
SM_DF_BC <-vegdist(SM_DF_taxa_clr, method = 'euclidean')
SM_DF_BC_pco.bray<-pco(SM_DF_BC,k=4)
SM_DF_BC_pco.bray_DF<-as.data.frame(SM_DF_BC_pco.bray$points)
SM_DF_BC_meta<-merge(SM_DF_BC_pco.bray_DF, metadata,by.x = 'row.names', by.y = 'row.names')

nt<-ggplot(SM_DF_BC_meta, aes(x=V1, y=V2, fill=Group)) +
  geom_point(size=2, pch=21,alpha=0.2)+
  geom_point(data=SM_DF_BC_meta %>% 
               group_by(Group) %>% 
               summarise_at(vars(matches("V")), mean),
             size=5, shape=21)+
  theme_classic()+
  theme(legend.position ="right",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black",size=0.5, fill = NA),
        axis.text.x= element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  ylab('PCO2')+
  xlab('PCO1')
nt

SM_DF_BC_pco.bray$eig[4]/sum(SM_DF_BC_pco.bray$eig)

#####Species_introduction
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_control_PA_community.txt", header = TRUE,row.names = 1)
SM_meta<-read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_metadata.txt", header = TRUE,row.names = 1)
SM_DF_meta<-merge(SM_DF, SM_meta,by.x = "row.names", by.y = "row.names")

write.csv(SM_DF_meta,"./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_PA/plotting.csv")
plotting_Cdiff <- read.csv("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_PA/plotting_Cdiff.csv", header = TRUE,check.names=FALSE)


gg<-ggplot(plotting_Cdiff, aes(x=as.numeric(DOL), y=as.numeric(Clostridioides_difficile), linetype=Antibiotic)) +
  #scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,80)) +
  #scale_x_reverse(limits=c(30,0)) +
  labs(y = "", x ="") +
  geom_smooth(method = "loess",color="#A7BEAE",size=1) +
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line = element_line(size=0.25),
        panel.background = element_blank())+
  scale_linetype_manual(values=c("dashed","dotted","solid","dotdash"))+
  coord_cartesian(ylim = c(0,1))
gg
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_PA/Cdiff_plot.png", gg, width=4, height=1.8, bg = "transparent")

plotting_Sepi <- read.csv("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_PA/plotting_Sepi.csv", header = TRUE,check.names=FALSE)


gg<-ggplot(plotting_Sepi, aes(x=as.numeric(DOL), y=as.numeric(Staphylococcus_epidermidis), linetype=Antibiotic)) +
  #scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,80)) +
  #scale_x_reverse(limits=c(30,0)) +
  labs(y = "", x ="") +
  geom_smooth(method = "loess",color="#B85042",size=1) +
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line.y =  element_line(size=0.25),axis.line.x = element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank())+
  scale_linetype_manual(values=c("twodash","solid","dotdash","longdash"))+
  coord_cartesian(ylim = c(0,1))
gg
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/ALL_PA/Sepi_plot.png", gg, width=4, height=1.8, bg = "transparent")





####Shift_persistence
SM_DF_BC <- vegdist(SM_DF, method = "bray")

species_BC_matrix <- as.matrix(SM_DF_BC)
species_BC_matrix[lower.tri(species_BC_matrix)] <- NA
species_BC_longform <- melt(species_BC_matrix, varnames = c("ID1", "ID2"))
SM_DF_matrix_meta1<-merge(species_BC_longform, metadata,by.x = "ID1", by.y = "row.names")
colnames(SM_DF_matrix_meta1)[4:ncol(SM_DF_matrix_meta1)] <- paste("1", colnames(metadata)[1:ncol(metadata)], sep = "_")
SM_DF_matrix_meta2<-merge(SM_DF_matrix_meta1, metadata,by.x = "ID2", by.y = "row.names")
colnames(SM_DF_matrix_meta2)[7:ncol(SM_DF_matrix_meta2)] <- paste("2", colnames(metadata)[1:ncol(metadata)], sep = "_")

SM_DF_matrix_meta2<-subset(SM_DF_matrix_meta2,SM_DF_matrix_meta2$'1_PatientID'==SM_DF_matrix_meta2$'2_PatientID')
write.csv(SM_DF_matrix_meta2,"./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_tax_shift_analysis/BC_shift.csv")




fisher <- matrix(c(15, 6, 12,9,2,7,25,40,35,37,38,36), nrow = 6)
fisher
fisher.test(fisher,alternative="two.sided")
