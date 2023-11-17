library(ggplot2)
culturing <- read.csv("Desktop/Projects/MDRO_UTI/181021_CulturingData5.csv", header = TRUE)
#x<-culturing[which(culturing$Categorie=="Non-recurrent")]

culturing$Patient_1<-factor(culturing$Patient, levels=c('4','6','5','13','14','16','18'))

ggplot(culturing[which(culturing$Categorie=="Recurrent"),], aes(x=Time, y=Detail, size=CFU.ml, col=Source)) + 
  geom_point() + theme_bw()+ 
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),strip.text = element_text(face="bold", size=12),panel.grid.major = element_blank(), panel.grid.minor = element_line(colour="black", size=0.5), panel.border = element_rect(colour = "black"), axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.text.y= element_text(size=12,face="bold"))  + 
  scale_size(name="CFU/ml",range = c(0, 20),trans="sqrt",breaks=c(0,10, 100,1000,10000,100000,1000000,5000000,10000000)) +  facet_grid(~Patient_1, scale="free_y", drop=TRUE) +
  scale_y_continuous(name="",labels=c("Episode 1","Episode 2","Episode 3"), breaks=c(1.5, 3.5, 5.5),minor_breaks = c(2.5,4.5)) +
  scale_x_continuous(name="",labels=c("DxU","Enrollment","0 days pAT","7 days pAT","30 days pAT","180 days pAT"),minor_breaks = c(0))+
  scale_color_manual(name="",values=c("#A97942","#FFD300"))+
  guides(colour = guide_legend(override.aes = list(size = 5)))



ggplot(culturing[which(culturing$Categorie=="Non-Recurrent"),], aes(x=Time, y=Detail, size=CFU.ml, col=Source)) + 
  geom_point() + theme_bw()+ 
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),strip.text = element_text(face="bold", size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.text.y= element_text(size=12,face="bold"))  + 
  scale_size(name="CFU/ml",range = c(0, 33),trans="sqrt",breaks=c(0,10, 100,1000,10000,100000,1000000,5000000,10000000)) +  facet_grid(~Patient, scale="free_y", drop=TRUE) +
  scale_y_continuous(name="",labels=c("Episode 1"), breaks=c(2.5)) +
  scale_x_continuous(name="",labels=c("DxU","Enrollment","0 days pAT","7 days pAT","30 days pAT","180 days pAT"))+
  scale_color_manual(name="",values=c("#A97942","#FFD300"))+
  guides(colour = guide_legend(override.aes = list(size = 5)))
