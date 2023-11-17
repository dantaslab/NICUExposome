####################################
############ Maaslin2 ##############
####################################

library("Maaslin2")

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/Maaslin/210923_maaslin_metadata_log.txt", header = TRUE,row.names = 1, na.strings = "NA")

SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/211006_rel_pathabundances.txt", header = TRUE,row.names = 1,check.names=FALSE)
SM_DF2<-t(SM_DF)
metadata<-read.delim('Desktop/Projects/NEC/final_analysis/taxonomy/210429_nec_metadata.txt', header = TRUE)
SM_DF_selection <- SM_DF2[row.names(SM_DF2) %in% metadata$SampleID, ]

fit_data = Maaslin2(
  input_data = SM_DF_selection, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  output = "./Desktop/Projects/NEC/final_analysis/humann/Maaslin/species_log", 
  fixed_effects = c("GA_birth","NEC_status","DOL","ABX_recent","Recent_penicillin","Recent_amingoglycoside",
                    "Recent_carbapenem","Recent_macrolide","Recent_1st_gen_cephalosporin","Recent_3rd_gen_cephalosporin",
                    "Recent_4th_gen_cephalosporin","Other_ANTIM_recent","Abx_cumulative","ANTIM_cumulative","mat_preeclampsia","inf",
                    "feeding_current","milk_current","milk_cumulative","formula_cumulative"),
  random_effects = c("Patient_ID"))


#######Genus
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/210201_metaphlan_table_genus.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  output = "./Desktop/Projects/NEC/final_analysis/taxonomy/Maaslin2/genus_log", 
  fixed_effects = c("GA_birth","NEC_status","DOL","ABX_recent","Recent_penicillin","Recent_amingoglycoside","Recent_lincosamide","Recent_glycopeptide",
                    "Recent_carbapenem","Recent_macrolide","Recent_1st_gen_cephalosporin","Recent_2nd_gen_cephalosporin","Recent_3rd_gen_cephalosporin",
                    "Recent_4th_gen_cephalosporin","Other_ANTIM_recent","Abx_cumulative","ANTIM_cumulative","mat_preeclampsia","route","intraventricular_hemorrhage",
                    "feeding_current","milk_current","milk_cumulative","formula_cumulative"),
  random_effects = c("Patient_ID"))

#######Family
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/210201_metaphlan_table_family.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  output = "./Desktop/Projects/NEC/final_analysis/taxonomy/Maaslin2/family_log", 
  fixed_effects = c("GA_birth","NEC_status","DOL","ABX_recent","Recent_penicillin","Recent_amingoglycoside","Recent_lincosamide","Recent_glycopeptide",
                    "Recent_carbapenem","Recent_macrolide","Recent_1st_gen_cephalosporin","Recent_2nd_gen_cephalosporin","Recent_3rd_gen_cephalosporin",
                    "Recent_4th_gen_cephalosporin","Other_ANTIM_recent","Abx_cumulative","ANTIM_cumulative","mat_preeclampsia","route","intraventricular_hemorrhage",
                    "feeding_current","milk_current","milk_cumulative","formula_cumulative"),
  random_effects = c("Patient_ID"))

#######Order
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/210201_metaphlan_table_order.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  output = "./Desktop/Projects/NEC/final_analysis/taxonomy/Maaslin2/order_log", 
  fixed_effects = c("GA_birth","NEC_status","DOL","ABX_recent","Recent_penicillin","Recent_amingoglycoside","Recent_lincosamide","Recent_glycopeptide",
                    "Recent_carbapenem","Recent_macrolide","Recent_1st_gen_cephalosporin","Recent_2nd_gen_cephalosporin","Recent_3rd_gen_cephalosporin",
                    "Recent_4th_gen_cephalosporin","Other_ANTIM_recent","Abx_cumulative","ANTIM_cumulative","mat_preeclampsia","route","intraventricular_hemorrhage",
                    "feeding_current","milk_current","milk_cumulative","formula_cumulative"),
  random_effects = c("Patient_ID"))

#######Class
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/210201_metaphlan_table_class.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  output = "./Desktop/Projects/NEC/final_analysis/taxonomy/Maaslin2/class_log", 
  fixed_effects = c("GA_birth","NEC_status","DOL","ABX_recent","Recent_penicillin","Recent_amingoglycoside","Recent_lincosamide","Recent_glycopeptide",
                    "Recent_carbapenem","Recent_macrolide","Recent_1st_gen_cephalosporin","Recent_2nd_gen_cephalosporin","Recent_3rd_gen_cephalosporin",
                    "Recent_4th_gen_cephalosporin","Other_ANTIM_recent","Abx_cumulative","ANTIM_cumulative","mat_preeclampsia","route","intraventricular_hemorrhage",
                    "feeding_current","milk_current","milk_cumulative","formula_cumulative"),
  random_effects = c("Patient_ID"))

#######Phylum
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/210201_metaphlan_table_phylum.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  output = "./Desktop/Projects/NEC/final_analysis/taxonomy/Maaslin2/phylum_log", 
  fixed_effects = c("GA_birth","NEC_status","DOL","ABX_recent","Recent_penicillin","Recent_amingoglycoside","Recent_lincosamide","Recent_glycopeptide",
                    "Recent_carbapenem","Recent_macrolide","Recent_1st_gen_cephalosporin","Recent_2nd_gen_cephalosporin","Recent_3rd_gen_cephalosporin",
                    "Recent_4th_gen_cephalosporin","Other_ANTIM_recent","Abx_cumulative","ANTIM_cumulative","mat_preeclampsia","route","intraventricular_hemorrhage",
                    "feeding_current","milk_current","milk_cumulative","formula_cumulative"),
  random_effects = c("Patient_ID"))


####################################
############# Figure ###############
####################################

library('pheatmap')

species_fig <- read.delim("./Desktop/Projects/NEC/final_analysis/taxonomy/Maaslin2/species_figure.txt", header = TRUE, na.strings = "NA")
species_DF <- cast(species_fig, category~species, value='value',fun.aggregate =mean)
row.names(species_DF)<-species_DF$category
species_DF<-species_DF[,-1]
species_DF[is.na(species_DF)]<- 0

paletteLength <- 50
myColor <- colorRampPalette(c("#1d3c45", "white", "#d2601a"))(paletteLength)

myBreaks <- c(seq(min(species_DF), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(species_DF)/paletteLength, max(species_DF), length.out=floor(paletteLength/2)))

# Plot the heatmap
pHEAT<-pheatmap(species_DF, color=myColor, breaks=myBreaks, cluster_rows = F, gaps_row = c(1,11,13,17,18))
ggsave('./Desktop/Projects/NEC/final_analysis/taxonomy/Maaslin2/heatmap.png', pHEAT, width=9, height=9, bg = 'transparent')

