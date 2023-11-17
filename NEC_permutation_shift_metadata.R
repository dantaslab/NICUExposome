shift_permutations <- read.table("Desktop/Projects/NEC/final_analysis/control_analysis/tax_shift_analysis/211130_permutation_IDs_for_taxa.txt", header = TRUE, row.names = 1, check.names=FALSE)
metadata_NEC <- read.table("Desktop/Projects/NEC/final_analysis/control_analysis/211117_metaphlan_table_controls.txt", header = TRUE, check.names=FALSE)


merged_metadata_frame_permutations<-merge(shift_permutations, metadata_NEC,by.x = "sample_I", by.y = "SampleID")
merged_metadata_frame_permutationsII<-merge(merged_metadata_frame_permutations, metadata_NEC,by.x = "sample_II", by.y = "SampleID")
write.table(merged_metadata_frame_permutationsII, "Desktop/Projects/NEC/final_analysis/control_analysis/tax_shift_analysis/211130_permutations_taxa_all.txt", sep="\t")

metadata_NEC <- read.table("Desktop/Projects/NEC/final_analysis/control_analysis/tax_shift_analysis/211130_permutations_metadata.txt", header = TRUE, row.names = 1, check.names=FALSE)
feeding_current<-aggregate(metadata_NEC$feeding_current, by=list(Category=metadata_NEC$permutation), FUN=sum)
milk_current<-aggregate(metadata_NEC$milk_current, by=list(Category=metadata_NEC$permutation), FUN=sum)
ABX_recent<-aggregate(metadata_NEC$ABX_recent, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_1st_gen_cephalosporin<-aggregate(metadata_NEC$Recent_1st_gen_cephalosporin, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_2nd_gen_cephalosporin<-aggregate(metadata_NEC$Recent_2nd_gen_cephalosporin, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_3rd_gen_cephalosporin<-aggregate(metadata_NEC$Recent_3rd_gen_cephalosporin, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_4th_gen_cephalosporin<-aggregate(metadata_NEC$Recent_4th_gen_cephalosporin, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_amingoglycoside<-aggregate(metadata_NEC$Recent_amingoglycoside, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_carbapenem<-aggregate(metadata_NEC$Recent_carbapenem, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_glycopeptide<-aggregate(metadata_NEC$Recent_glycopeptide, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_lincosamide<-aggregate(metadata_NEC$Recent_lincosamide, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_macrolide<-aggregate(metadata_NEC$Recent_macrolide, by=list(Category=metadata_NEC$permutation), FUN=sum)
Recent_penicillin<-aggregate(metadata_NEC$Recent_penicillin, by=list(Category=metadata_NEC$permutation), FUN=sum)
ANTIM_cumulative<-aggregate(metadata_NEC$ANTIM_cumulative, by=list(Category=metadata_NEC$permutation), FUN=mean)
Abx_cumulative<-aggregate(metadata_NEC$Abx_cumulative, by=list(Category=metadata_NEC$permutation), FUN=mean)
formula_cumulative<-aggregate(metadata_NEC$formula_cumulative, by=list(Category=metadata_NEC$permutation), FUN=mean)
milk_cumulative<-aggregate(metadata_NEC$milk_cumulative, by=list(Category=metadata_NEC$permutation), FUN=mean)
metadata_summary<-cbind(feeding_current,milk_current,ABX_recent,Recent_1st_gen_cephalosporin,Recent_2nd_gen_cephalosporin,
      Recent_3rd_gen_cephalosporin,Recent_4th_gen_cephalosporin,Recent_amingoglycoside,Recent_carbapenem,
      Recent_glycopeptide,Recent_lincosamide,Recent_macrolide,Recent_penicillin,ANTIM_cumulative,
      Abx_cumulative,formula_cumulative,milk_cumulative)

write.table(metadata_summary, "Desktop/Projects/NEC/final_analysis/control_analysis/tax_shift_analysis/211130_permutations_metadata_summary.txt", sep="\t")

p_value <- read.table("Desktop/test.txt", check.names=FALSE)
p.adjust(p_value$V1,method="fdr")
