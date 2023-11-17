
library(nlme)
library("multcomp")

iRep_averagepersample <- read.delim("iRep/iRep_averagepersample.txt", header = TRUE)

iRep_averagepersample_week<-iRep_averagepersample[which(iRep_averagepersample$dayspreonset<=7),]

iRep_averagepersample_late<-iRep_averagepersample[which(iRep_averagepersample$Onset=='late'&iRep_averagepersample$dayspreonset<=7),]
iRep_averagepersample_early<-iRep_averagepersample[which(iRep_averagepersample$Onset=='early'&iRep_averagepersample$dayspreonset<=7),]


model.b = lme(iRepaveragefortimepoint ~ NEC*dayspreonset,
              random = ~1|iRep_Patient, data=iRep_averagepersample_week)
summary(model.b)
model.null = lme(iRepaveragefortimepoint ~ dayspreonset,
                 random = ~1|iRep_Patient, data=iRep_averagepersample_week)
anova(update(model.b, . ~ ., method = 'ML'),
      update(model.null, . ~ ., method = 'ML'))

summary(glht(model.b, lincfit=mcp(NEC='Tukey')))

model.b$terms

iRep_graph<-ggplot(iRep_averagepersample_week, aes(x = dayspreonset, y = iRepaveragefortimepoint, col=NEC)) + 
  geom_point(pch='.',alpha=0.6)+
  geom_smooth(method="loess",size=1)+
  scale_x_continuous(limit=c(0,7),breaks=c(0,2,4,6))+
  scale_y_continuous(limit=c(0,2.5),breaks=c(0,1,2))+
  theme_classic()+
  ylab("iRep")+
  xlab("days_pre_onset")+
  theme_classic()+
  theme(legend.text=element_text(size=12, face='bold'),legend.title=element_blank(),legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.line=element_line(size=0.5),
        axis.text.x= element_blank(),axis.ticks.y=element_line(size=0.5),axis.ticks.x=element_line(size=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  scale_color_manual(values=c("#3b4d61","#ef9d10"))
iRep_graph
ggsave("figures/iRep.png", iRep_graph, width=5, height=3, bg = "transparent")



iRep_averagepersample <- read.delim("Desktop/Projects/NEC/final_analysis/iRep/iRep_averagepersample_Enterobacteriaceae.txt", header = TRUE)

iRep_Entero_7days<-iRep_averagepersample[which(iRep_averagepersample$dayspreonset<=7),]
iRep_averagepersample_late<-iRep_averagepersample[which(iRep_averagepersample$Onset=='Late'&iRep_averagepersample$dayspreonset<=7),]
iRep_averagepersample_early<-iRep_averagepersample[which(iRep_averagepersample$Onset=='Early'&iRep_averagepersample$dayspreonset<=7),]


model.b = lme(iRepaveragefortimepoint ~ NEC*dayspreonset,
              random = ~1|iRep_Patient, data=iRep_Entero_7days)
summary(model.b)
model.null = lme(iRepaveragefortimepoint ~ dayspreonset,
                 random = ~1|iRep_Patient, data=iRep_Entero_7days)
anova(update(model.b, . ~ ., method = 'ML'),
      update(model.null, . ~ ., method = 'ML'))

summary(glht(model.b, lincfit=mcp(NEC_status='Tukey')))

model.b$terms

iRep_graph<-ggplot(iRep_Entero_7days, aes(x = dayspreonset, y = iRepaveragefortimepoint, col=NEC)) + 
  geom_point(pch='.',alpha=0.6)+
  geom_smooth(method="loess",size=1.5)+
  theme_classic()+
  ylab("iRep")+
  xlab("days_pre_onset")+
  theme_classic()+
  theme(legend.text=element_text(size=12, face='bold'),legend.title=element_blank(),legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.line=element_line(size=1),
        axis.text.x= element_blank(),axis.ticks=element_line(size=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  scale_color_manual(values=c("#3b4d61","#ef9d10"))+
  scale_x_continuous(limit=c(0,7),breaks=c(0,2,4,6))+
  scale_y_continuous(limit=c(0,8))
iRep_graph
ggsave("./Desktop/Projects/NEC/final_analysis/figures/iRep_Enterobacteriaceae.png", iRep_graph, width=10, height=5, bg = "transparent")







