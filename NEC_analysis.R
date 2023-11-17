###### Packages ######
library(ggplot2)
library(reshape)
library(labdsv)
library(vegan)
library(ggpubr)
library(permute)
library(dplyr)
library(rsample)
library(purrr)
library(tidyr)
library(rowr)
library(rcompanion)
library(ggpubr)
library("BiodiversityR")
library("Maaslin2")

###### Load DFs #######
MetaPhlan_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/210201_metaphlan_table.txt", header = TRUE, check.names = FALSE, row.names=1)
MetaPhlan_meta<- read.delim("Desktop/Projects/NEC/final_analysis/taxonomy/210429_nec_metadata.txt", header = TRUE, check.names = FALSE, row.names=1)

###### Diversity #######
Metaphlan_SD <-as.data.frame(diversity(MetaPhlan_DF, index="shannon"))

colnames(Metaphlan_SD)[1] <- "ShannonDiversity"
Metaphlan_SD_meta <- merge(Metaphlan_SD, MetaPhlan_meta,by.x = 'row.names', by.y = "row.names")
row.names(Metaphlan_SD_meta)<-Metaphlan_SD_meta[,1]
Metaphlan_SD_meta<-Metaphlan_SD_meta[,-1]

SD<-ggplot(Metaphlan_SD_meta, aes(x=DOL, y=ShannonDiversity, fill=as.factor(NEC_status))) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=0.8)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=1), 
        axis.text.x= element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  ylab("Shannon Diversity")+
  xlab("DOL")+
  ylim(0,2.5)+
  scale_fill_manual(values=c('#0d1137','#e52165'))+
  geom_smooth(aes(group=as.factor(NEC_status)),method="loess",color="white")
SD
ggsave('./Desktop/Projects/NEC/final_analysis/taxonomy/plots/SD_DOL_limited.png', SD, width=12, height=5, bg = 'transparent')

onset_DF<- read.delim("Desktop/Projects/NEC/final_analysis/taxonomy/onset.txt", header = TRUE, check.names = FALSE)

onset<-ggplot(onset_DF, aes(x=as.factor(Dummy), y=NEC_onset,fill=as.factor(NEC_status))) +
  theme_bw()+
  geom_jitter(position=position_jitter(0.3), pch=21, size=3)+
  geom_boxplot(color="black",fill="white",alpha=0.2,size=1, outlier.shape=NA)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=1), 
        axis.text.x= element_blank(),axis.title.x = element_blank(),axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  ylab("Shannon Diversity")+
  xlab("DOL")+
  ylim(1,81)+
  scale_fill_manual(values=c('#e52165'))+coord_flip()
onset
ggsave('./Desktop/Projects/NEC/final_analysis/taxonomy/plots/onset.png', onset, width=12, height=1.3, bg = 'transparent')

####SUBSET EARLY/LATE####
Metaphlan_SD_meta_early<-Metaphlan_SD_meta[which(Metaphlan_SD_meta$NEC_onset<40),]
Metaphlan_SD_meta_late<-Metaphlan_SD_meta[which(Metaphlan_SD_meta$NEC_onset>40),]

SD<-ggplot(Metaphlan_SD_meta, aes(x=DOL, y=ShannonDiversity, fill=as.factor(NEC_status))) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=0.8)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=1), 
        axis.text.x= element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  ylab("Shannon Diversity")+
  xlab("DOL")+
  ylim(0,2.5)+
  scale_fill_manual(values=c('#0d1137','#e52165'))+
  geom_smooth(aes(group=as.factor(NEC_status)),method="loess",color="white")
SD
