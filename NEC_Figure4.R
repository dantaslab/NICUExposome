library(scales)
library(vegan)
library("labdsv")
library(ggplot2)

DNA_metaphlan<-read.delim("Desktop/Projects/NEC/final_analysis/datasets/DNA_metaphlan_relAbundance.txt", header = TRUE, sep='\t',check.names = FALSE,row.names = 1)
NEC_metadata<-read.delim("Desktop/Projects/NEC/final_analysis/metadata/210429_nec_metadata.txt", header = TRUE, sep='\t',check.names = FALSE,row.names = 1)

RNA_metaphlan<-read.delim("Desktop/Projects/NEC/final_analysis/datasets/RNA_metaphlan_relAbundance.txt", header = TRUE, sep='\t',check.names = FALSE,row.names = 1)


##############################
###### NEC_composition #######
##############################

##############
## Taxonomy ##
##############

#### DNA #####

DNA_metaphlan_BC <-vegdist(DNA_metaphlan, method = "bray")

DNA_metaphlan_pco.bray<-pco(DNA_metaphlan_BC,k=4)
DNA_metaphlan_pco.bray_DF<-as.data.frame(DNA_metaphlan_pco.bray$points)
NEC_metaphlan_BC_meta <- merge(DNA_metaphlan_pco.bray_DF, NEC_metadata,by.x = 'row.names', by.y = "DNA_Sample_ID")

#########PCoA############
g<-ggplot(NEC_metaphlan_BC_meta, aes(x=V1, y=V2, fill=as.factor(NEC_status))) +geom_point(size=5, pch=21)+
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=1), 
        axis.text.x= element_blank(),axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  ylab("PCO2")+
  xlab("PCO1")+
  scale_fill_manual(values=c("#3b4d61","#ef9d10"))+
  coord_fixed(1/1.25)
g
ggsave("./Desktop/Projects/NEC/final_analysis/figures/PCoA_metaphlan.png", g, width=10, height=5, bg = "transparent")

##### RNA #####

RNA_metaphlan_BC <-vegdist(RNA_metaphlan, method = "bray")

RNA_metaphlan_pco.bray<-pco(RNA_metaphlan_BC,k=4)
RNA_metaphlan_pco.bray_DF<-as.data.frame(RNA_metaphlan_pco.bray$points)
NEC_RNA_metaphlan_BC_meta <- merge(RNA_metaphlan_pco.bray_DF, NEC_metadata,by.x = 'row.names', by.y = "DNA_Sample_ID")

#########PCoA############
g<-ggplot(NEC_RNA_metaphlan_BC_meta, aes(x=V1, y=V2, fill=as.factor(NEC_status))) +geom_point(size=5, pch=21)+
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=1), 
        axis.text.x= element_blank(),axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  ylab("PCO2")+
  xlab("PCO1")+
  scale_fill_manual(values=c("#3b4d61","#ef9d10"))+
  coord_fixed(1/1.25)
g
ggsave("./Desktop/Projects/NEC/final_analysis/figures/PCoA_RNA_metaphlan.png", g, width=10, height=5, bg = "transparent")


#################
### Functions ###
#################

#### DNA #####
DNA_humann<-read.delim("Desktop/Projects/NEC/final_analysis/datasets/DNA_humann_relAbundance.txt", header = TRUE, sep='\t',check.names = FALSE,row.names = 1)

DNA_humann_BC <-vegdist(DNA_humann, method = "bray")

DNA_humann_pco.bray<-pco(DNA_humann_BC,k=4)
DNA_humann_pco.bray_DF<-as.data.frame(DNA_humann_pco.bray$points)
DNA_humann_BC_meta <- merge(DNA_humann_pco.bray_DF, NEC_metadata,by.x = 'row.names', by.y = "DNA_Sample_ID")

#########PCoA############
g<-ggplot(DNA_humann_BC_meta, aes(x=V1, y=V2, fill=as.factor(NEC_status))) +geom_point(size=5, pch=21)+
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=1), 
        axis.text.x= element_blank(),axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  ylab("PCO2")+
  xlab("PCO1")+
  scale_fill_manual(values=c("#3b4d61","#ef9d10"))+
  coord_fixed(1/1.1)
g
ggsave("./Desktop/Projects/NEC/final_analysis/figures/PCoA_DNA_humann.png", g, width=10, height=5, bg = "transparent")

##### RNA #####
RNA_humann<-read.delim("Desktop/Projects/NEC/final_analysis/datasets/RNA_humann_relAbundance.txt", header = TRUE, sep='\t',check.names = FALSE,row.names = 1)

RNA_humann_BC <-vegdist(RNA_humann, method = "bray")

RNA_humann_pco.bray<-pco(RNA_humann_BC,k=4)
RNA_humann_pco.bray_DF<-as.data.frame(RNA_humann_pco.bray$points)
NEC_RNA_humann_BC_meta <- merge(RNA_humann_pco.bray_DF, NEC_metadata,by.x = 'row.names', by.y = "DNA_Sample_ID")

#########PCoA############
g<-ggplot(NEC_RNA_humann_BC_meta, aes(x=V1, y=V2, fill=as.factor(NEC_status))) +geom_point(size=5, pch=21)+
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black",size=1), 
        axis.text.x= element_blank(),axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())  +
  ylab("PCO2")+
  xlab("PCO1")+
  scale_fill_manual(values=c("#3b4d61","#ef9d10"))+
  coord_fixed(1/1.25)
g
ggsave("./Desktop/Projects/NEC/final_analysis/figures/PCoA_RNA_humann.png", g, width=10, height=5, bg = "transparent")


############################
###### NEC_diversity #######
############################

DNA_metaphlan_SD<- as.data.frame(diversity(DNA_metaphlan, index="shannon"))
colnames(DNA_metaphlan_SD)[1] <- "ShannonDiversity"
NEC_metaphlan_SD_meta <- merge(DNA_metaphlan_SD, NEC_metadata,by.x = 'row.names', by.y = "DNA_Sample_ID")

NEC_div_plot <-ggplot(NEC_metaphlan_SD_meta, aes(x=DOL, y=ShannonDiversity, col=as.factor(NEC_status))) +
  geom_point(pch='.',alpha=0.6)+
  geom_smooth(method="loess",size=1.5)+theme_bw()+
  ylab("")+
  xlab("Day of life")+
  theme_classic()+
  theme(legend.text=element_text(size=12, face='bold'),legend.title=element_blank(),legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.line=element_line(size=1),
        axis.text.x= element_blank(),axis.ticks.y=element_line(size=1),axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  scale_y_continuous(limit=c(0,2.5),breaks=c(0,1,2),oob=squish)+
  scale_color_manual(values=c("#3b4d61","#ef9d10"))
NEC_div_plot 
ggsave("./Desktop/Projects/NEC/final_analysis/figures/Shannon_diversity_DOL.png", NEC_div_plot, width=10, height=5, bg = "transparent")


NEC_onset<-NEC_metaphlan_SD_meta[which(NEC_metaphlan_SD_meta$NEC_status=="1"),]
NEC_onset<-NEC_onset[which(NEC_onset$NEC_onset<=80),]
NEC_onset_unique<-unique(NEC_onset[,c('Patient','NEC_onset','NEC_status')])

NEC_onset_plot <-ggplot(NEC_onset_unique, aes(y=NEC_onset,x=as.factor(NEC_status))) +
  geom_violin(size=1)+geom_jitter(pch=21,size=4,fill="#ef9d10")+
  geom_boxplot(width=0.2, color="black",outlier.shape = NA,fill=NA,size=1)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme_classic()+
  theme(legend.text=element_text(size=12, face='bold'),legend.title=element_blank(),legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.line=element_line(size=1),
        axis.text.x= element_blank(),axis.ticks.x=element_line(size=1),axis.ticks.y=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+coord_flip()+
  scale_y_continuous(limit=c(0,80),breaks=c(0,20,40,60,80))
NEC_onset_plot 
ggsave("./Desktop/Projects/NEC/final_analysis/figures/Shannon_onset_DOL.png", NEC_onset_plot, width=10, height=1.5, bg = "transparent")

###########################
####### Enterotype ########
###########################

Late_onset <-
  matrix(c(0,0,0,1,2,4,3,0,2,1,3,0,0,5,3,6),
         nrow = 2)
fisher.test(Late_onset, alternative = "two.sided")

Early_onset <-
  matrix(c(1,1,1,6,5,3,5,11,2,1,1,13,8,10,11,22),
         nrow = 2)
fisher.test(Early_onset, alternative = "two.sided")

############################
##### Species_plotting #####
############################

species_plotting <- read.delim("taxonomy/Maaslin2/species_plotting.txt", header = TRUE,row.names = 1)

Klebsiella<-ggplot(species_plotting, aes(x=Days_pre_onset, y=Klebsiella,col=Case_control,line))+
  geom_point(size=1, pch=20,alpha=0)+xlab('Days pre onset')+geom_smooth(aes(linetype = Time),method="loess")+
  scale_color_manual(values = c("#ef9d10","#3b4d61"))+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.border = element_blank(),panel.grid.major.y = element_line(size=0.25,color="#818181"),panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks =  element_line(size=0.75),
        axis.title.y = element_blank(),axis.line =  element_line(size=0.75),
        panel.background = element_blank())+
  ylab('Klebsiella')+
  scale_x_continuous(breaks=c(0,20,40,60),limits=c(0,60))+
  scale_y_continuous(breaks=c(0,50,100),limits=c(-5,100),oob=squish)
Klebsiella
ggsave("figures/Klebsiella_abundance.png", Klebsiella, width=5, height=3.5, bg = "transparent")



V_parvula<-ggplot(species_plotting, aes(x=Days_pre_onset, y=Veillonella_parvula,col=Case_control,line))+
  geom_point(size=1, pch=20,alpha=0)+xlab('Days pre onset')+geom_smooth(aes(linetype = Time),method="loess")+
  scale_color_manual(values = c("#ef9d10","#3b4d61"))+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.border = element_blank(),panel.grid.major.y = element_line(size=0.25,color="#818181"),panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks.y =  element_line(size=0.75),axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),axis.line.y =  element_line(size=0.75),
        panel.background = element_blank())+
  ylab('Klebsiella')+
  scale_x_continuous(breaks=c(0,20,40,60),limits=c(0,60))+
  scale_y_continuous(breaks=c(0,1,10,100),limits=c(-5,100),oob=squish,trans = scales::pseudo_log_trans())
V_parvula

ggsave("figures/V_parvula_abundance.png", V_parvula, width=5, height=4.7, bg = "transparent")


####################
###### Shifts ######
####################

shift_plot <- read.delim("shifts/shift_plot.txt", header = TRUE,row.names = 1)


hist(shift_plot$DOL,col="#B85042",breaks=40,xlim=c(1,80))


shift<-ggplot(shift_plot, aes(x = relative_onset, y = Group, fill=Group)) + 
  geom_density_ridges()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.border = element_blank(),panel.grid.major.y = element_line(size=0.25,color="#818181"),panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks.x =  element_line(size=0.5),axis.ticks.y =  element_blank(),
        axis.title.y = element_blank(),axis.line =  element_line(size=0.5),
        panel.background = element_blank())+xlim(0,80)+
  scale_fill_manual(values=c("#ef9d10","#3b4d61"))
shift

ggsave("figures/shifts.png", shift, width=5, height=3, bg = "transparent")

wilcox.test(relative_onset ~ Group, data=shift_plot)
  
  
  