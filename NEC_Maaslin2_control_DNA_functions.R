####################################
############ Maaslin2 ##############
####################################

library("Maaslin2")

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/maaslin_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_base_functions/ALL_control_human_community_matchingRNA.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.01,
  output = "./Desktop/Projects/NEC/final_analysis/control_analysis/RNA_Maaslin2/results_DNA_DOL_log", 
  fixed_effects = c("DOL","Recent_1st_gen_cephalosporin","Recent_3rd_gen_cephalosporin","Recent_amingoglycoside","Recent_carbapenem",
                    "Recent_penicillin","Recent_lincosamide","Other_ANTIM_recent","formula_cumulative","milk_cumulative",
                    "recent_Calcium_Gluconate","recent_Chlorathiazide","recent_Cholecalciferol","recent_Dexamethasone","recent_Dopamine",
                    "recent_Fentanyl","recent_Hydrocortisone","recent_Indomethacin","recent_Iron","recent_Midazolam","recent_Morphine",
                    "recent_MViron","recent_Phenobarbitol","recent_Surfactant","recent_SodiumChloride","recent_SodiumBicarbonate","recent_PotassiumChloride"
                    ),
  random_effects = c("Patient"))


####################################
############### PCA ################
####################################
library("compositions")

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/211117_Maaslin_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/function_PERMANOVA/211128_control_human_community.txt", header = TRUE,row.names = 1)



DNA_F_clr<-clr(SM_DF)

head(DNA_F_clr)





####################################
############# Figure ###############
####################################

library('pheatmap')

species_fig <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/Maaslin2_cases_DOL_prev_all/211127_species_figure.txt", header = TRUE, na.strings = "NA")
species_DF <- cast(species_fig, category~species, value='value',fun.aggregate =mean)
row.names(species_DF)<-species_DF$category
species_DF<-species_DF[,-1]
species_DF[is.na(species_DF)]<- 0


paletteLength <- 50
myColor <- colorRampPalette(c("#1d3c45", "white", "#d2601a"))(paletteLength)

myBreaks <- c(seq(min(species_DF), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(species_DF)/paletteLength, max(species_DF), length.out=floor(paletteLength/2)))

# Plot the heatmap
pHEAT<-pheatmap(species_DF, color=myColor, breaks=myBreaks, cluster_rows = F, gaps_row = c(1,2,11))
ggsave('./Desktop/Projects/NEC/final_analysis/taxonomy/Maaslin2/heatmap.png', pHEAT, width=9, height=9, bg = 'transparent')

