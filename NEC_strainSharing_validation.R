###################################
## NEC strain sharing validation ##
###################################
library("ggplot2")

C_diff_hist <- read.delim("./Desktop/inStrain/plotting/C_diff_hist.txt", header = TRUE,row.names = 1, na.strings = "NA")
hist(C_diff_hist$popANI,ylim=c(0,20),col="#A7BEAE",breaks=100,xlim=c(0.9965,1))
hist(C_diff_hist$popANI,ylim=c(0,20),col="#A7BEAE",breaks=10000,xlim=c(0.99999,1))
hist(C_diff_hist$percent_compared,col="#A7BEAE",breaks=40,xlim=c(0.2,1))


C_diff_plot <- read.delim("~/Desktop/inStrain/plotting/C_diff_plot.txt", header = TRUE, na.strings = "NA")

gplot_points<-ggplot(C_diff_plot, aes(x=x_SNPs_log, y=y_value,fill=type)) +
  geom_point(size=4,pch=21,stroke=1)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = "transparent"),
        axis.text.x= element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  geom_vline(xintercept=-0.379231294)+scale_fill_manual(values=c("#1e847f","#1e847f","#ecc19c","#ecc19c"))

gplot_points
ggsave("./Desktop/inStrain/plotting/C_diff_plot.png", gplot_points, width=4, height=6, bg = "transparent")


S_epidermidis_hist <- read.delim("~/Desktop/inStrain/plotting/S_epidermidis_hist.txt", header = TRUE,row.names = 1, na.strings = "NA")
hist(S_epidermidis_hist$popANI,ylim=c(0,50),col="#B85042",breaks=100)
hist(S_epidermidis_hist$popANI,ylim=c(0,30),col="#B85042",breaks=50000,xlim=c(0.99999,1))
hist(S_epidermidis_hist$percent_compared,ylim=c(0,30),col="#B85042",breaks=40,xlim=c(0,1))


S_epi_plot <- read.delim("./Desktop/inStrain/plotting/S_epi_plot_filtered.txt", header = TRUE, na.strings = "NA")

gplot_points<-ggplot(S_epi_plot, aes(x=x_SNPs_log, y=y_value,fill=type)) +
  geom_point(size=4,pch=21,stroke=1)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = "transparent"),
        axis.text.x= element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  geom_vline(xintercept=-0.379231294)+scale_fill_manual(values=c("#1e847f","#1e847f","#ecc19c","#ecc19c"))

gplot_points
ggsave("./Desktop/inStrain/plotting/S_epi_plot.png", gplot_points, width=4, height=12, bg = "transparent")
