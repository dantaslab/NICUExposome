SM_DF_shift <- read.delim("./Desktop/Projects/MDRO_UTI_COHORT/Metagenomics/200824_metaphlan_shift_data.txt", header = TRUE)
row.names(SM_DF_shift)<-SM_DF_shift[,1]
SM_DF_shift<-SM_DF_shift[,-1]

###### Load metadata ######
SM_DF_shift_meta<-read.delim("./Desktop/Projects/MDRO_UTI_COHORT/Metagenomics/200824_shift_metadata.txt", header = TRUE,check.names = FALSE)

SM_DF_shift_BC <-vegdist(SM_DF_shift, method = "bray")

SM_shift_pco.bray<-pco(SM_DF_shift_BC,k=4)
SM_shift_pco.bray_DF<-as.data.frame(SM_shift_pco.bray$points)
SM_BC_shift_meta<-merge(SM_shift_pco.bray_DF, SM_DF_shift_meta,by.x = 'row.names', by.y = "Sample_ID")

#########HISTORY############
g<-ggplot(SM_BC_shift_meta, aes(x=V3, y=V4, fill=Group)) +geom_point(size=5, pch=21)+
  theme_bw()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"), 
        axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1,face="bold"),
        axis.title.x = element_text(size=14, face="bold",margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=14, face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y= element_text(size=12,face="bold"))  +
  ylab("PCO2")+
  xlab("PCO1")+
  coord_fixed((SM_shift_pco.bray$eig[3]/sum(SM_shift_pco.bray$eig))/(SM_shift_pco.bray$eig[4]/sum(SM_shift_pco.bray$eig)))
g

betadispresults<-betadisper(SM_DF_shift_BC,SM_DF_shift_meta$Patient, type = "centroid")
centroids.btwn<-betadispresults$centroids

btw_grps<-read.delim("./Desktop/treatmentbtw_grps.txt", header = TRUE,check.names = FALSE)
perMOVbtwn<-adonis(centroids.btwn~treatment+patient,btw_grps, method = "euclidean",permutations=999)
perMOVbtwn

adonis(SM_DF_shift_BC~Group+Patient,data=SM_DF_shift_meta,method="bray",permutations=999, add="cailliez")

subject<-SM_DF_shift_meta$Patient
subject_data<-read.delim("./Desktop/Projects/MDRO_UTI_COHORT/Metagenomics/subject_metadata.txt", header = TRUE,check.names = FALSE, row.names=1)
sample_data<-read.delim("./Desktop/Projects/MDRO_UTI_COHORT/Metagenomics/sample_metadata.txt", header = TRUE,check.names = FALSE, row.names=1)


PERMANOVA_repeat_measures(
  SM_DF_shift_BC,
  subject, subject_data = subject_data,
  sample_data = sample_data,
  metadata_order = c(names(subject_data), names(sample_data)),
  permutations=999, ncores=1)

