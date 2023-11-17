#####################################
######### NEC enterotypes ###########
#####################################

library('ggplot2')
library('ggridges')

NEC_enterotype <- read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/DMM/dmn_SampleEnterotype.txt", header = TRUE,row.names = 1, na.strings = "NA")

NEC_enterotype$Enterotype_relabel  = factor(NEC_enterotype$Enterotype_relabel, levels=c("E8_new", "E7_new", "E6_new","E5_new", "E4_new", "E3_new","E2_new","E1_new"))

panelA<-ggplot(NEC_enterotype, aes(x=Days_pre_onset, y=Enterotype, fill=NEC_status)) +
  theme_bw()+
  geom_density_ridges2(aes(height =..ndensity..),scale = 0.9,bandwidth = 1,alpha=0.6,cex=0.5)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.border = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(),
        axis.text.x= element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y=element_blank(),
        axis.text.y= element_blank())+scale_x_reverse(limits = c(120, 0))+
  scale_fill_manual(values=c("#e52165","#0d1137"))
panelA
ggsave('./Desktop/Projects/NEC/final_analysis/taxonomy/DMM/density_preonset.png', panelA, width=5, height=4, bg = 'transparent')


