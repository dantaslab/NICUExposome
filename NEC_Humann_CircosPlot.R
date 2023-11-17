###################################################
######## Humann3 - Maaslin visualization ##########
###################################################

library("circlize")

circos_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/211026_circos_humann_late_2.txt", header = TRUE, na.strings = "NA")

col_species_trans_1 = t_col("#FEB24C",perc=40)
col_species_trans_2 = t_col("#88419D",perc=40)
col_species_trans_3 = t_col("#0570B0",perc=40)
col_species_trans_4 = t_col("#74A9CF",perc=40)
col_species_trans_5 = t_col("#deebf7",perc=40)
col_species = c(col_species_trans_1,col_species_trans_2,col_species_trans_3,col_species_trans_4,col_species_trans_5)
segment_col=c("black")

circos.par("track.height" = 0.6,cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(circos_DF$species, x = circos_DF$pseudo_x)
circos.track(circos_DF$species, y = circos_DF$coeff,ylim=c(-1.3,1.3),bg.col=col_species)
col = c("#0D1137", "#E52165", "#E52165", "#E52165", "#E52165")
circos.segments(circos_DF$pseudo_x,circos_DF$stdev_low,circos_DF$pseudo_x,circos_DF$stdev_bar,sector.index = "Enterococcus faecalis")
circos.segments(circos_DF$pseudo_x_bar1,circos_DF$stdev_low1,circos_DF$pseudo_x_bar1,circos_DF$stdev_bar1,sector.index = "Escherichia coli")
circos.segments(circos_DF$pseudo_x_bar2,circos_DF$stdev_low2,circos_DF$pseudo_x_bar2,circos_DF$stdev_bar2,sector.index = "Klebsiella pneumoniae")
circos.segments(circos_DF$pseudo_x_bar3,circos_DF$stdev_low3,circos_DF$pseudo_x_bar3,circos_DF$stdev_bar3,sector.index = "Klebsiella quasipneumoniae")
circos.segments(circos_DF$pseudo_x_bar4,circos_DF$stdev_low4,circos_DF$pseudo_x_bar4,circos_DF$stdev_bar4,sector.index = "Klebsiella variicola")
circos.trackLines(circos_DF$species, circos_DF$pseudo_x,circos_DF$pseudo_y_line)
circos.trackPoints(circos_DF$species, circos_DF$pseudo_x, bg=col,circos_DF$coeff, pch = 21, cex = 1)
circos.segments(52.5,-1.3,52.5,1.3,sector.index = "Enterococcus faecalis",col=segment_col)
circos.segments(60.5,-1.3,60.5,1.3,sector.index = "Enterococcus faecalis",col=segment_col)
circos.segments(63.5,-1.3,63.5,1.3,sector.index = "Enterococcus faecalis",col=segment_col)
circos.segments(69.5,-1.3,69.5,1.3,sector.index = "Klebsiella pneumoniae",col=segment_col)
circos.segments(101.5,-1.3,101.5,1.3,sector.index = "Klebsiella pneumoniae",col=segment_col)
circos.segments(117.5,-1.3,117.5,1.3,sector.index = "Klebsiella pneumoniae",col=segment_col)
circos.segments(118.5,-1.3,118.5,1.3,sector.index = "Klebsiella pneumoniae",col=segment_col)
circos.segments(15.5,-1.3,15.5,1.3,sector.index = "Klebsiella quasipneumoniae",col=segment_col)
circos.segments(16.5,-1.3,16.5,1.3,sector.index = "Klebsiella quasipneumoniae",col=segment_col)
circos.segments(8.5,-1.3,8.5,1.3,sector.index = "Klebsiella variicola",col=segment_col)
circos.segments(14.5,-1.3,14.5,1.3,sector.index = "Klebsiella variicola",col=segment_col)




###################################
############ Early ################
###################################

circos_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/humann/211026_circos_humann_early.txt", header = TRUE, na.strings = "NA")

col_species_trans_1 = t_col("#8C96C6",perc=40)
col_species_trans_2 = t_col("#E31A1C",perc=40)
col_species_trans_3 = t_col("#818181",perc=40)

col_species = c(col_species_trans_1,col_species_trans_2,col_species_trans_3)
segment_col=c("black")

circos.par("track.height" = 0.6,cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(circos_DF$species, x = circos_DF$pseudo_x)
circos.track(circos_DF$species, y = circos_DF$coeff,ylim=c(-1.3,1.3),bg.col=col_species)
col = c("#E52165", "#0D1137", "#0D1137")
circos.segments(circos_DF$pseudo_x,circos_DF$stdev_low,circos_DF$pseudo_x,circos_DF$stdev_bar,sector.index = "Enterobacter cloacae complex")
circos.segments(circos_DF$pseudo_x_bar1,circos_DF$stdev_low1,circos_DF$pseudo_x_bar1,circos_DF$stdev_bar1,sector.index = "Staphylococcus epidermidis")
circos.segments(circos_DF$pseudo_x_bar2,circos_DF$stdev_low2,circos_DF$pseudo_x_bar2,circos_DF$stdev_bar2,sector.index = "unclassified")
circos.trackLines(circos_DF$species, circos_DF$pseudo_x,circos_DF$pseudo_y_line)
circos.trackPoints(circos_DF$species, circos_DF$pseudo_x, bg=col,circos_DF$coeff, pch = 21, cex = 1)
circos.segments(2.5,-1.3,2.5,1.3,sector.index = "unclassified",col=segment_col)












