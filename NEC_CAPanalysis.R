#####################################
######### NEC CAP analysis ##########
#####################################

library(ggplot2)
library(labdsv)
library(vegan)
library(ggpubr)
library(reshape2)
library(reshape)
library(permute)
library(dplyr)
library(rsample)
library(purrr)
library(tidyr)
library(rowr)
library(rcompanion)


NEC_metaphlan <- read.delim("Downloads/210201_metaphlan_table (3).txt", header = TRUE, row.names = 1, check.names = FALSE)
NEC_metadata <- read.delim("Downloads/210429_nec_metadata (1).txt", header = TRUE)
row.names(NEC_metadata)<-NEC_metadata$DNA_Sample_ID

vare.cap <- capscale(NEC_metaphlan ~ sqrtDOL, NEC_metadata, dist="bray")
#vare.cap
#anova(vare.cap)
#plot(vare.cap)
smry_cap <- summary(vare.cap)
CAP_for_plot  <- data.frame(smry_cap$sites[,1:4])
CAP_for_plot$Sample_ID<-row.names(CAP_for_plot)
Merged_CAP<-merge(CAP_for_plot, NEC_metadata,by.x = "Sample_ID", by.y = "DNA_Sample_ID")

ggplot(Merged_CAP, aes(x=CAP1, y=MDS1, fill=sqrtDOL)) +geom_point(size=5,pch=21)+
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_rect(colour = "black"), 
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank(), axis.ticks = element_blank())  