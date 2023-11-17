#####################################
######## NEC summary figure #########
#####################################

library('ggplot2')
library('ggridges')

NEC_summary <- read.delim("./Desktop/Projects/NEC/final_analysis/summary/Fig1_DNA_samples.txt", header = TRUE,row.names = 1, na.strings = "NA")

NEC_summary$Type  = factor(NEC_summary$Type, levels=c("Cytokine profiles", "Metatranscriptomics", "Metagenomics"))

panelA<-ggplot(NEC_summary, aes(x=Days_pre_onset, y=Type, fill=NEC_status)) +
  theme_bw()+
  geom_density_ridges2(aes(height =..ndensity..),scale = 0.9,bandwidth = 1,alpha=0.6,cex=0.5)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.border = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(),
        axis.text.x= element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y=element_blank(),
        axis.text.y= element_blank())+scale_x_reverse(limits = c(120, 0))+
  scale_fill_manual(values=c("#e52165","#0d1137"))
ggsave('./Desktop/Projects/NEC/final_analysis/summary/summary.png', panelA, width=8, height=4, bg = 'transparent')
