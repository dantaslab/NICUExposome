library("ggplot2")
library("reshape")

ARG_uniqueARG <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/shortbred/shift_vs_permutation_uniqueARGs.txt", header = TRUE,check.names=FALSE)

boxplot_ARGs<-ggplot(ARG_uniqueARG, aes(x=Case, y=UniqueARGs,fill=Case))+
  geom_boxplot(alpha=0.75)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black",size=0.5, fill = NA),
        axis.text.x= element_blank(),axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  ylab('')+
  xlab('')+
  scale_fill_manual(values=c("#0D1137","#ed3572"))+scale_y_continuous(limits=c(0,100),breaks=c(0,30,60,90))

boxplot_ARGs

ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_tax_shift_analysis/uniqueARGs.png", boxplot_ARGs, width=1, height=2, bg = "transparent")


wilcox.test(ARG_uniqueARG$UniqueARGs ~ ARG_uniqueARG$Case,data=ARG_uniqueARG)


####variance around shifts
finalpermutationtaxadataset_controlsamples <-read.delim("Box/2022_RobertEric_NECMultiOmics/EK R scripts/Microbiome shifts in preterm IGM (Metaphlan and ARGs)/211201_finalpermutationtaxadataset_controlsamples.txt", header = TRUE, sep = "\t")

gp<-ggplot(finalpermutationtaxadataset_controlsamples, aes(x=LogFC, y=0, fill=Type)) +
  geom_density_ridges2(scale = 10, alpha = 0.6) +
  facet_grid(Variance2 ~ ., scales="free") +
  #scale_fill_manual(values=c("#4077e6","#0D1137")) +
  #scale_y_continuous(breaks=c(0,8))+
  scale_fill_manual(values=c("#0D1137","#ed3572")) +
  labs(y = "Density", x ="Fold change in consecutive-sample species abundance (log)") +
  theme(legend.title = element_blank()) +
  xlim(-6,6) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #panel.background = element_blank(), axis.line = element_line(colour = "black"))
  theme(legend.position="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line = element_line(size=0.25),
        panel.background = element_blank(),axis.ticks.y=element_blank())
gp

ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_tax_shift_analysis/changing_taxa.png", gp, width=3, height=5, bg = "transparent")



shifts <-read.delim("Box/2022_RobertEric_NECMultiOmics/EK R scripts/Microbiome shifts in preterm IGM (Metaphlan and ARGs)/211202_shiftinvasionhistogram.txt", header = TRUE, sep = "\t")

intro_plots<-ggplot(shifts, aes(x=Permutation_value)) +
  geom_histogram(binwidth=1, boundary=0) +
  #facet_grid(Order ~ ., scales="free") +
  facet_grid(Order ~ .) +
  geom_vline(aes(xintercept=Shift_vertline),
             color="#ed3572", size=1)+
  #scale_fill_manual(values=c("#4077e6","#0D1137")) +
  scale_y_continuous(breaks=c(0,300), position="right")+
  #scale_fill_manual(values=c("#0D1137","#4077e6")) +
  labs(y = "Frequency", x ="Acquisition events") +
  theme(legend.title = element_blank()) +
  #xlim(0,25) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #panel.background = element_blank(), axis.line = element_line(colour = "black"))
  theme(legend.position="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line = element_line(size=0.25),
        panel.background = element_blank(),axis.ticks.y=element_blank(),
        strip.text.y = element_blank()) 


ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_tax_shift_analysis/introduced_taxa.png", intro_plots, width=3, height=5, bg = "transparent")

#######ARG_persistence
ARG_uniqueARG_persist <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/shortbred/uniqueARG_persistence.txt", header = TRUE,check.names=FALSE)

gg<-ggplot(ARG_uniqueARG_persist, aes(x=as.numeric(Shift_DOL.1), y=as.numeric(percentage_persistence))) +
  geom_point(size=2, pch=21,alpha=0.2)+
  #scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,60)) +
  #scale_x_reverse(limits=c(30,0)) +
  labs(y = "", x ="") +
  geom_smooth(method = "loess",color="black",size=1) +
  theme()
gg


#######Taxa_persistence
Taxa_persist <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_tax_shift_analysis/220424_shift_taxa_persistence.txt", header = TRUE,check.names=FALSE)
Taxa_persist_1<-subset(Taxa_persist,Taxa_persist$Species!="Klebsiella_pneumoniae")


gg<-ggplot(Taxa_persist_1, aes(x=as.numeric(DOL_to_shift), y=as.numeric(Presence),col=Species)) +
  facet_grid(Species ~ .) +
  geom_point(size=2, pch=21,alpha=0)+
  scale_y_continuous(limits=c(0,1.3),breaks=c(0,0.5,1)) +
  scale_x_continuous(limits=c(0,10), breaks=c(0,5,10)) +
  #scale_x_reverse(limits=c(30,0)) +
  labs(y = "", x ="") +
  geom_smooth(method = "loess",size=1,fill="lightgrey") +
  theme()+
  coord_cartesian(ylim = c(0,1))+
  scale_color_manual(values=c("#aed6dc","#81b7d2","#4a536b","#B85042"))+
  theme(legend.position="none",
        panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour="#D9D9D9"),panel.grid.major.x =element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line = element_line(size=0.25),
        panel.background = element_blank(),axis.ticks.y=element_blank(),
        strip.text.y = element_blank(),panel.spacing = unit(1, "lines")) 
gg
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_tax_shift_analysis/introduced_taxa_persistence.png", gg, width=3, height=2, bg = "transparent")


#######Shift persistence
shift_persist <- read.csv("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_tax_shift_analysis/BC_shift.csv", header = TRUE,check.names=FALSE)


gg<-ggplot(shift_persist, aes(x=as.numeric(DOL_post_shift_normalized), y=as.numeric(value))) +
  geom_point(size=2, pch=21,alpha=0)+
  #scale_x_reverse(limits=c(30,0)) +
  labs(y = "", x ="") +
  geom_smooth(method = "loess",size=1,fill="lightgrey") +
  theme()+
  coord_cartesian(ylim = c(0,1))+
  theme() 
gg






ARG_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/shortbred/shortbred_matrix.txt", header = TRUE,check.names=FALSE,row.names = 1)
ARG_DF_richness<-as.data.frame(rowSums(ARG_DF[,1:ncol(ARG_DF)]>0.0))
colnames(ARG_DF_richness)[1] <-"Richness"

metadata <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_tax_shift_analysis/220422_metadata_shift_1.txt", header = TRUE,check.names=FALSE,row.names = 1)

ARG_richness_meta<-merge(ARG_DF_richness, metadata,by.x = "row.names", by.y = "row.names")
ARG_DF_meta_shifts<-subset(ARG_richness_meta,ARG_richness_meta$Shift_DOL!="NA")
write.csv(ARG_DF_meta_shifts,"./Desktop/Projects/NEC/final_analysis/control_analysis/shortbred/richness_shifts.csv")

gg<-ggplot(ARG_DF_meta_shifts, aes(x=as.numeric(Shift_DOL.1), y=as.numeric(Richness))) +
  geom_point(size=2, pch=21,alpha=0.2)+
  #scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(-2,2)) +
  #scale_x_reverse(limits=c(30,0)) +
  labs(y = "", x ="") +
  geom_smooth(method = "loess",color="#A7BEAE",size=1) +
  theme()
gg


metadata <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_tax_shift_analysis/220422_metadata_permutations_1.txt", header = TRUE,check.names=FALSE,row.names = 1)

ARG_richness_meta<-merge(ARG_DF_richness, metadata,by.x = "row.names", by.y = "row.names")
ARG_DF_meta_shifts<-subset(ARG_richness_meta,ARG_richness_meta$Shift_DOL!="NA")
write.csv(ARG_DF_meta_shifts,"./Desktop/Projects/NEC/final_analysis/control_analysis/shortbred/richness_permutations.csv")

gg<-ggplot(ARG_DF_meta_shifts, aes(x=as.numeric(Shift_DOL.1), y=as.numeric(Richness))) +
  geom_point(size=2, pch=21,alpha=0.2)+
  #scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(-2,2)) +
  #scale_x_reverse(limits=c(30,0)) +
  labs(y = "", x ="") +
  geom_smooth(method = "loess",color="#A7BEAE",size=1) +
  theme()
gg
