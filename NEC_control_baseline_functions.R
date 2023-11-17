#########################################
######### Baseline development ##########
#########################################


SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base_functions/ALL_control_human_community_matchingRNA.txt", header = TRUE,check.names=FALSE,row.names = 1)

metadata <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base_functions/metadata_STL.txt", header = TRUE,check.names=FALSE,row.names = 1)

SM_DF_richness<-as.data.frame(rowSums(SM_DF[,1:ncol(SM_DF)]>0.0))
colnames(SM_DF_richness)[1] <-"Richness"
SM_DF_richness_meta<-merge(SM_DF_richness, metadata,by.x = 'row.names', by.y = "row.names")

rare_R_SB<-ggplot(SM_DF_richness_meta, aes(x=DOL, y=Richness))+
  geom_point(size=2, pch=21,alpha=0)+
  theme_classic()+
  theme(legend.text=element_text(size=14, face='bold'),legend.title=element_text(size=14, face='bold'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  ylab('Richness')+
  xlab('DOL')+geom_smooth(method="loess",col="#ecc19c")+
  scale_y_continuous(limits=c(0,500),expand = c(0.01,0.01))+
  scale_x_continuous(limits=c(0,80),expand = c(0.01,0.01))
rare_R_SB
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base/Richness_functions.png", rare_R_SB, width=5, height=3, bg = "transparent")

####################################
####### Species richness/PWY #######
####################################

SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/RNA_average_pathway_species_richness.txt", header = TRUE,check.names=FALSE,row.names = 1)

metadata <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base_functions/metadata_STL.txt", header = TRUE,check.names=FALSE,row.names = 1)

SM_DF_richness_meta<-merge(SM_DF, metadata,by.x = 'row.names', by.y = "row.names")

rare_R_SB<-ggplot(SM_DF_richness_meta, aes(x=DOL, y=average_richness))+
  geom_point(size=2, pch=21,alpha=0)+
  theme_classic()+
  theme(legend.text=element_text(size=14, face='bold'),legend.title=element_text(size=14, face='bold'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),panel.background = element_rect(fill = "transparent"),
        axis.text.x= element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  ylab('Richness')+
  xlab('DOL')+geom_smooth(method="loess",col="#1e847f")+
  scale_x_continuous(limits=c(0,80),expand = c(0.01,0.01))+
  scale_y_continuous(limits=c(0,10),expand = c(0.01,0.01))
rare_R_SB
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/RNA_average_species_richness_functions.png", rare_R_SB, width=5, height=3, bg = "transparent")


SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/DNA_average_pathway_species_richness.txt", header = TRUE,check.names=FALSE,row.names = 1)

metadata <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base_functions/metadata_STL.txt", header = TRUE,check.names=FALSE,row.names = 1)

SM_DF_richness_meta<-merge(SM_DF, metadata,by.x = 'row.names', by.y = "row.names")

rare_R_SB<-ggplot(SM_DF_richness_meta, aes(x=DOL, y=Average_richness))+
  geom_point(size=2, pch=21,alpha=0)+
  theme_classic()+
  theme(legend.text=element_text(size=14, face='bold'),legend.title=element_text(size=14, face='bold'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  ylab('Richness')+
  xlab('DOL')+geom_smooth(method="loess",col="#ecc19c")+
  scale_x_continuous(limits=c(0,80),expand = c(0.01,0.01))+
  scale_y_continuous(limits=c(0,10),expand = c(0.01,0.01))
rare_R_SB
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/DNA_average_pathway_richness_functions.png", rare_R_SB, width=5, height=3, bg = "transparent")

######################################
####### RNA overrepresentation #######
######################################

overrepresentation <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/RNA_overrepresentation.txt", header = TRUE,check.names=FALSE,row.names = 1)

rare_R_SB<-ggplot(overrepresentation, aes(x=DNA_median, y=RNA_median, size=Prevalence,fill=Genus))+
  geom_point(pch=21)+
  theme_classic()+
  theme(legend.text=element_text(size=14, face='bold'),legend.title=element_text(size=14, face='bold'),legend.position="right",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.line=element_line(size=1),axis.ticks=element_line(size=1),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())+
  ylab('')+xlim(0,25)+ylim(0,25)+
  xlab('')+ geom_abline(intercept = 0, slope = 1,size=1)+
  scale_size_continuous(range=c(2.5,22),breaks=c(0.01,0.25,0.5,0.75,1),limits=c(0,1))+
  scale_fill_manual(values=c("#d72631","#f3ca20","#d9a5b3","#aed6dc","#81b7d2","#4a536b","#ff6e40","#77c593","#B85042","#ff6e40","#500472"))+
  guides(fill = guide_legend(override.aes = list(size = 5)))
rare_R_SB
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/RNA_overrepresentation.png", rare_R_SB, width=9, height=9, bg = "transparent")



SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/RNA_overrepresentation_species.txt", header = TRUE,check.names=FALSE)

SM_DF$Species <- factor(SM_DF$Species, levels=c("Escherichia_coli","Enterobacter_cloacae_complex","Enterococcus_faecalis","Staphylococcus_epidermidis","Klebsiella_pneumoniae",
                                                "Clostridium_perfringens","Veillonella_parvula"))

rare_R_SB<-ggplot(SM_DF, aes(x=DOL, y=value,color=metric))+
  geom_point(size=2, pch=21,alpha=0)+
  theme_classic()+
  theme(legend.text=element_text(size=14, face='bold'),legend.title=element_text(size=14, face='bold'),legend.position="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.line=element_line(size=1),axis.ticks=element_line(size=1),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        strip.text.y = element_blank())+
  ylab('Richness')+
  xlab('DOL')+geom_smooth(method="loess")+facet_grid(Species ~ ., scales='free_y')+
  scale_color_manual(values=c("#ecc19c","#1e847f"))+
  coord_cartesian(ylim=c(-10,110))+scale_y_continuous(breaks=c(0,100))+scale_x_continuous(limits=c(0,80),breaks=c(0,20,40,60,80))
rare_R_SB

ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/RNA_overrepresentation_by_species.png", rare_R_SB, width=9, height=9, bg = "transparent")


############################################
### Bubble plot change over time summary ###
############## Expression ##################
############################################
library("ggthemes")

bubble_plot_functions <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/Cat2_fig2A.txt", header = TRUE,check.names=FALSE)
bubble_plot_functions2 <- bubble_plot_functions[which(bubble_plot_functions$Total_pathways>3),]

xx = ggplot(bubble_plot_functions2, aes(x = Percentage_flex, y = Prevalence)) + 
  geom_point(aes(size = Normalized_expression, fill = Cat1), alpha = 0.75, shape = 21) + 
  labs( x= "", y = "", size = "Normalized Expression (%)", fill = "")  + 
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="right",legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line = element_line(size=0.75),
        panel.background = element_blank()) +
  scale_size_continuous(range=c(2.5,22),breaks=c(0.1,1,5,10,15))+
  geom_rangeframe()+
  scale_x_continuous(limits=c(0,100), breaks=c(0,25,50,75,100))+
  scale_y_continuous(limits=c(0,110), breaks=c(0,25,50,75,100))+
  scale_fill_manual(values=c("#ecc19c","#1e847f","#000000"))+
  coord_cartesian(ylim = c(0,100))
xx
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/fig2A.png", xx, width=6, height=6, bg = "transparent")

richness <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/pathway_taxa_richness.txt", header = TRUE,check.names=FALSE)

shapiro.test(richness$Richness)
wilcox.test(Richness ~ Grp, data = richness)
density_plot<-ggplot(richness, aes(x=Richness, fill=Grp)) +
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),axis.line = element_line(size=0.6),
        panel.background = element_blank())+
  geom_density(alpha=0.5)+scale_x_continuous(limits=c(1,15), breaks=c(1,5,10,15))+scale_y_continuous(limits=c(0,0.35), breaks=c(0,0.1,0.2,0.3))+
  scale_fill_manual(values=c("#aed6dc","#ff9a8d"))
density_plot
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/fig2B.png", density_plot, width=4, height=1.5, bg = "transparent")

#Underrepresented taxa core/flex
stool_col_type_s_GO_ORA <- read.delim('./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/Species_enrichment_core_flex.txt', check.names = FALSE, header = TRUE, row.names = 1)
s_col_s<-as.data.frame(apply(stool_col_type_s_GO_ORA, 1, 
                             function(x) {
                               tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                               fisher.test(tbl, alternative="two.sided")$p.value
                             }))
s_col_s$p.adjust<-sapply(s_col_s,p.adjust,method="fdr")
write.table(s_col_s, "./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/species_enrichment_results.txt", sep="\t")


volcano <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/volcano_plot_drivers.txt", header = TRUE,check.names=FALSE)


volcano_plot<-ggplot(volcano, aes(x=FC, y=FDR, fill=Genus)) +
  geom_point(pch=21,size=5)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line = element_line(size=1.5),
        panel.background = element_blank())+
  scale_fill_manual(values=c("#d72631","#d9a5b3","#f3ca20","#aed6dc","#81b7d2","#4a536b","#77c593","#FFFFFF","#B85042","#ff6e40"))+
  scale_x_continuous(limits=c(-13,5), breaks=c(-10,-5,0,5))+geom_hline(yintercept = 1.301029996)
volcano_plot
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/fig2B_volcano.png", volcano_plot, width=6, height=5, bg = "transparent")



library("scatterpie")

fig2c <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/figure_2C.txt", header = TRUE,check.names=FALSE)
fig2c$x1<-as.numeric(fig2c$x1_fake)

pie<-ggplot() + geom_scatterpie(aes(x=x1_fake, y=y1_fake, group=Group, r=log), data=fig2c, cols=LETTERS[1:10])+
                           coord_equal()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="right",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),panel.background = element_rect(fill = "transparent"),
        axis.text.x= element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.line = element_line(size=1.1))+
  scale_fill_manual(values=c("#d9a5b3","#ff6e40","#81b7d2","#77c593","#aed6dc","#f3ca20","#4a536b","#B85042","#d72631","#FFFFFF"))+ 
  geom_hline(yintercept = 0)+scale_x_continuous(limits=c(1.5,30.75),breaks=c(3,6,9,12,15,18,21,24,27,30))+
  scale_y_continuous(limits=c(-4,6))+geom_scatterpie_legend(fig2c$log, x=27, y=-1, n=3, labeller=function(x) x)
pie
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/fig2C.png", pie, width=9, height=7, bg = "transparent")



fig2c_error <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/figure_2C_range.txt", header = TRUE,check.names=FALSE)

pie_range<-ggplot(fig2c_error, aes(x=x1_fake))+
  geom_linerange(aes(ymin=min_fake,ymax=max_fake))+
  coord_equal()+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks = element_blank(),
        axis.title.y = element_blank(),axis.line = element_blank(),
        panel.background = element_blank())+
  scale_x_continuous(limits=c(1.5,30.75))+scale_y_continuous(limits=c(-4,6))
pie_range
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/fig2C_background.png", pie_range, width=9, height=7, bg = "transparent")


###########################################
######### Heatmap Transcription ###########
###########################################
heatmap_enriched <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/enrichment_metadata.txt", header = TRUE,check.names=FALSE)

heatmap_enriched$metadata <- factor(heatmap_enriched$metadata, levels=c("Other_ANTIM_recent","Recent_lincosamide","Recent_amingoglycoside","Recent_carbapenem","Recent_3rd_gen_cephalosporin","Recent_1st_gen_cephalosporin","Recent_penicillin",
                                                                      "recent_Dopamine","recent_Hydrocortisone","recent_Dexamethasone","recent_Indomethacin","recent_Midazolam","recent_Phenobarbitol","recent_Morphine","recent_Fentanyl",
                                                                      "recent_Chlorathiazide","recent_Calcium_Gluconate","recent_PotassiumChloride","recent_SodiumBicarbonate","recent_SodiumChloride","recent_Iron","recent_MViron","recent_Cholecalciferol",
                                                                      "recent_Surfactant","formula_cumulative","milk_cumulative"))


ggplot(heatmap_enriched, aes(x =variable , y = metadata, fill =log2_value)) +
  geom_tile(color = "#818181",
            lwd = 0.5,
            linetype = 1)+
  facet_grid(group ~ ., scales='free_y', space='free')+ scale_fill_gradientn(colors=c("white", "#408ec6", "#1e2761","#7a2048"))
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks = element_blank(),
        axis.title.y = element_blank(),axis.line = element_blank(),
        panel.background = element_blank())



  heatmap_depletion <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/depletion_metadata.txt", header = TRUE,check.names=FALSE)
  
  heatmap_depletion$metadata <- factor(heatmap_depletion$metadata, levels=c("Other_ANTIM_recent","Recent_lincosamide","Recent_amingoglycoside","Recent_carbapenem","Recent_3rd_gen_cephalosporin","Recent_1st_gen_cephalosporin","Recent_penicillin",
                                                                          "recent_Dopamine","recent_Hydrocortisone","recent_Dexamethasone","recent_Indomethacin","recent_Midazolam","recent_Phenobarbitol","recent_Morphine","recent_Fentanyl",
                                                                          "recent_Chlorathiazide","recent_Calcium_Gluconate","recent_PotassiumChloride","recent_SodiumBicarbonate","recent_SodiumChloride","recent_Iron","recent_MViron","recent_Cholecalciferol",
                                                                          "recent_Surfactant","formula_cumulative","milk_cumulative"))
  
  
  ggplot(heatmap_depletion, aes(x =variable , y = metadata, fill =log2_value)) +
    geom_tile(color = "#818181",
              lwd = 0.5,
              linetype = 1)+
    facet_grid(group ~ ., scales='free_y', space='free')+ scale_fill_gradientn(colors=c("#F2F2F2", "#d2601a","#d71b3b","#000000"))



  
fig_2E <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/Fig_2E.txt", header = TRUE,check.names=FALSE)
  
fig_2E$metadata <- factor(fig_2E$metadata, levels=c("DOL","dummy","dummy2","dummy3","dummy4","dummy5","dummy6","Recent_3rd_gen_cephalosporin","Recent_carbapenem","Other_ANTIM_recent",
                                                    "recent_Calcium_Gluconate","recent_SodiumChloride","recent_Chlorathiazide","recent_Dexamethasone","recent_Hydrocortisone","recent_Cholecalciferol","formula_cumulative"))
  
fig_2E$pwy <- factor(fig_2E$pwy, levels=c("Cofactor, Prosthetic Group, Electron Carrier, and Vitamin Biosynthesis ","Fatty Acid and Lipid Biosynthesis ","Amino Acid Biosynthesis ",
                                          "Nucleoside and Nucleotide Biosynthesis ","Aromatic Compound Degradation ","Cell Structure Biosynthesis ","Fermentation ",
                                          "Secondary Metabolite Degradation ","Carbohydrate Biosynthesis ","Amino Acid Degradation ","Carbohydrate Degradation ","Carboxylate Degradation ",
                                          "Other Generation of Precursor Metabolite and Energy ","TCA cycle ","Nucleoside and Nucleotide Degradation ","Secondary Metabolite Biosynthesis ",
                                                    "Inorganic Nutrient Metabolism ","Respiration ","C1 Compound Utilization and Assimilation ","Fatty Acid and Lipid Degradation "))

  
f_2e<-ggplot(fig_2E, aes(x = coeff, y = y_value, fill =metadata, size=pwy_number, color =metadata,alpha=metadata,shape=metadata)) +
  geom_point()+
  geom_vline(xintercept = 0,size=1,color="black")+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="right",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks.y = element_blank(),axis.ticks.x=element_line(size=1,color="black"),
        axis.title.y = element_blank(),axis.line.y = element_blank(),axis.line.x=element_line(size=1,color="black"),
        strip.background = element_blank(),
        strip.text.y = element_blank(),panel.spacing = unit(1, "lines"))+
  facet_grid(grp ~ ., scales='free_y', space='free')+
  scale_size_continuous(range=c(2,10), breaks=c(1,10,20,30,40,50))+scale_x_continuous(limits=c(-1.1,1), breaks=c(-1,-0.5,0,0.5,1))+
  scale_fill_manual(values=c("white","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF",
                             "#e75874","#be1558","#fbcbc9",
                             "#aed6dc","#6883bc","#3a6b35","#e3b448",
                             "#cbd18f","#1868ae","#322514"))+
  scale_color_manual(values=c("black","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF",
                             "black","black","black",
                             "black","black","black","black",
                             "black","black","black"))+
  scale_alpha_manual(values=c(1,0,0,0,0,0,0,
                              1,1,1,1,1,1,1,1,1,1))+
  scale_shape_manual(values=c(23,21,21,21,21,21,
                              21,21,21,21,21,
                              21,21,21,21,21,
                              21))
  


f_2e
  
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/fig2e.png", f_2e, width=9, height=7, bg = "transparent")


#scale_y_continuous(breaks=c(1,2,3,4,5,6,7,10,11,12,13,14,15,16,17,18,21,22,23),labels=c("Cofactor, Prosthetic Group, Electron Carrier, and Vitamin Biosynthesis ","Fatty Acid and Lipid Biosynthesis ","Amino Acid Biosynthesis ",
#                                                                                        "Nucleoside and Nucleotide Biosynthesis ","Cell Structure Biosynthesis ","Carbohydrate Biosynthesis ","Secondary Metabolite Biosynthesis ",
#                                                                                        "Aromatic Compound Degradation ","Secondary Metabolite Degradation ","Amino Acid Degradation ","Carbohydrate Degradation ","Carboxylate Degradation ",
#                                                                                        "Nucleoside and Nucleotide Degradation ","Inorganic Nutrient Metabolism ","C1 Compound Utilization and Assimilation ","Fatty Acid and Lipid Degradation ",
#                                                                                        "Fermentation ","TCA cycle ","Respiration "))



heatmap_DNARNA <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/Fig_2F.txt", header = TRUE,check.names=FALSE)

heatmap_DNARNA$metadata <- factor(heatmap_DNARNA$metadata, levels=c("Other_ANTIM_recent","Recent_carbapenem","Recent_3rd_gen_cephalosporin",
                                                                          "recent_Dopamine","recent_Hydrocortisone","recent_Dexamethasone","recent_Indomethacin","recent_Midazolam","recent_Phenobarbitol","recent_Morphine","recent_Fentanyl",
                                                                          "recent_Chlorathiazide","recent_Calcium_Gluconate","recent_PotassiumChloride","recent_SodiumBicarbonate","recent_SodiumChloride","recent_Iron","recent_MViron","recent_Cholecalciferol",
                                                                          "recent_Surfactant","formula_cumulative","milk_cumulative"))

heatmap_DNARNA$med_grp<-factor(heatmap_DNARNA$med_grp, levels=c("diet","med","antimicrobial"))

ggplot(heatmap_DNARNA, aes(x =metadata , y = grp, fill =percentage)) +
  geom_tile(color = "#818181",
            lwd = 0.5,
            linetype = 1)+
  facet_grid(cat1 ~ med_grp, scales='free', space='free')+ 
  scale_fill_gradientn(colors=c("#7a2048", "#ffcce7","#322514"),na.value = "white")+
  theme()




################################
######### PROCRUSTES ###########
################################

#### DNA_functions
SM_DF_functions <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base_functions/ALL_control_human_community_matchingRNA.txt", header = TRUE,row.names = 1)

#CLR transform functional data for PCA
DNA_F_DNA_clr<-clr(SM_DF_functions)

DNA_F_DNA_clr.rda <- prcomp(DNA_F_DNA_clr)

#### RNA_functions
SM_DF_functions_RNA <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base_functions/RNA_humann3_community_relAbund.txt", header = TRUE,row.names = 1)

#CLR transform functional data for PCA
DNA_F_RNA_clr<-clr(SM_DF_functions_RNA)

DNA_F_RNA_clr.rda <- prcomp(DNA_F_RNA_clr)

pro <- procrustes(X = DNA_F_DNA_clr.rda, Y = DNA_F_RNA_clr.rda,symmetric=TRUE)
pro

protest(X = DNA_F_DNA_clr.rda, Y = DNA_F_RNA_clr.rda, scores = "sites", permutations = 999)

ctest <- data.frame(rda1=pro$Yrot[,1],
                    rda2=pro$Yrot[,2],xrda1=pro$X[,1],
                    xrda2=pro$X[,2])

ctest_DF<-as.data.frame(ctest)
#ctest_DF_meta<-merge(ctest_DF, SM_meta,by.x = "row.names", by.y = "row.names")

procrustes<-ggplot(ctest_DF) +
  geom_point(aes(x=xrda1, y=xrda2),pch=21,size=3,fill="#ecc19c") +
  geom_point(aes(x=rda1, y=rda2),pch=21,size=3,fill="#1e847f",) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2),lwd=0.05)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black",size=0.5, fill = NA),
        axis.text.x= element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())
procrustes


ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/figs2b.png", procrustes, width=8, height=8, bg = "transparent")


#### DNA_metaphlan
SM_DF_functions <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_PERMANOVA/ALL_community_metaphlan_matchingRNA.txt", header = TRUE,row.names = 1)

#CLR transform functional data for PCA
DNA_F_DNA_clr<-clr(SM_DF_functions)

DNA_F_DNA_clr.rda <- prcomp(DNA_F_DNA_clr)

#### RNA_functions
SM_DF_functions_RNA <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base_functions/RNA_humann3_community_relAbund.txt", header = TRUE,row.names = 1)

#CLR transform functional data for PCA
DNA_F_RNA_clr<-clr(SM_DF_functions_RNA)

DNA_F_RNA_clr.rda <- prcomp(DNA_F_RNA_clr)

pro <- procrustes(X = DNA_F_DNA_clr.rda, Y = DNA_F_RNA_clr.rda,symmetric=TRUE)
pro

protest(X = DNA_F_DNA_clr.rda, Y = DNA_F_RNA_clr.rda, scores = "sites", permutations = 999)

ctest <- data.frame(rda1=pro$Yrot[,1],
                    rda2=pro$Yrot[,2],xrda1=pro$X[,1],
                    xrda2=pro$X[,2])

ctest_DF<-as.data.frame(ctest)
#ctest_DF_meta<-merge(ctest_DF, SM_meta,by.x = "row.names", by.y = "row.names")

procrustes<-ggplot(ctest_DF) +
  geom_point(aes(x=xrda1, y=xrda2),pch=21,size=3,fill="black",alpha=0.2) +
  geom_point(aes(x=rda1, y=rda2),pch=21,size=3,fill="#1e847f",) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2),lwd=0.05)+
  theme_classic()+
  theme(legend.position ="none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black",size=0.5, fill = NA),
        axis.text.x= element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())
procrustes


ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/figs2c.png", procrustes, width=8, height=8, bg = "transparent")



####################################
############ Maaslin2 ##############
####################################

library("Maaslin2")

SM_meta <- read.delim("./Desktop/metadata_CG.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("./Desktop/test1.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.01,
  output = "./Desktop/test_maaslin_CG", 
  fixed_effects = c("DOL","Group"
  ),
  random_effects = c("Patient"))



                 
                 
