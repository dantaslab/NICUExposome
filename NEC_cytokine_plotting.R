NEC_cytokines <- read.delim("Desktop/Projects/200701_NEC_cytokines.txt", header = TRUE)
NEC_cytokines_metadata <- read.delim("Desktop/Projects/200701_cytokine_metadata.txt", header =TRUE)
NEC_merge<-merge(NEC_cytokines, NEC_cytokines_metadata,by.x = "Sample_ID", by.y = "Sample_ID")


ggplot(NEC_merge, aes(x=Days_preonset, y=log(concentration+1),col=as.factor(NEC_status))) +
  geom_point() +
  theme_bw()+
  facet_wrap(~Cytokine, scales="free")+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold"))  +
  ylab("Concentration")+
  xlim(0,10)+
  geom_smooth( method = "loess")

NEC_chemokines <- read.delim("Desktop/Projects/200701_chemokines.txt", header = TRUE)
NEC_merge_chemo<-merge(NEC_chemokines, NEC_cytokines_metadata,by.x = "Sample_ID", by.y = "Sample_ID")


ggplot(NEC_merge_chemo, aes(x=Days_preonset, y=concentration,col=as.factor(NEC_status))) +
  geom_point() +
  theme_bw()+
  facet_wrap(~Cytokine, scales="free")+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold"))  +
  ylab("Concentration")+
  xlim(0,10)+
  geom_smooth( method = "loess")
