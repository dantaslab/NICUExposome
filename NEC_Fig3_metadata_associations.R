######################################
####### Associations over DOL ########
######################################
library(scales)
library(ggplot2)

# Taxa
fisher <- matrix(c(22, 7, 17,11,4,7,22,39,30,35,36,36), nrow = 6)
fisher.test(fisher,alternative="two.sided",simulate.p.value=TRUE)

# Functions
fisher <- matrix(c(16,7,16,10,4,7,24,39,32,36,36,36), nrow = 6)
fisher.test(fisher,alternative="two.sided",simulate.p.value=TRUE)

####################################
######### Maaslin2 acute ###########
####################################

library("Maaslin2")

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_functions/DOL_10_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_functions/10_control_human_community.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.01,
  output = "./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_functions/10_acute", 
  fixed_effects = c("DOL_log","Caffeine_recent","Dopamine_recent","Famotidine_recent","Sufactant",
                    "MV_iron_recent","VitaminA_recent","Recent_1st_gen_cephalosporin","Recent_3rd_gen_cephalosporin",
                    "Recent_amingoglycoside","Recent_carbapenem","Recent_lincosamide","Recent_penicillin",
                    "Other_ANTIM_recent","formula_cumulative","milk_cumulative","GA_birth"
  ),
  random_effects = c("Patient_ID","Study"))


###################################
###### Fig 2B acute effects #######
###################################


fig_2B <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/220523_acute_effects_selection.txt", header = TRUE,check.names=FALSE)

fig_2B$metadata <- factor(fig_2B$metadata, levels=c("Sufactant","famotidine_recent","Dexamethasone_recent","fentanyl","Midazolam_recent","Dopamine_recent","PotassiumChloride_recent",
                                                    "VitA_recent","MV_iron_recent","iron_recent","Recent_lincosamide","Recent_glycopeptide","Recent_amingoglycoside","Recent_penicillin",
                                                    "Carbapenem","Recent_4th_gen_cephalosporin","Recent_3rd_gen_cephalosporin","formula","milk_cumulative","GA_birth"
                                                    ))

fig_2B_plot<-ggplot(fig_2B, aes(x = coef, y = metadata, fill=genus)) +
  geom_point(pch=21,size=5,alpha=0.75)+
  geom_vline(xintercept = 0,size=1,color="black")+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="top",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank(),
        axis.text= element_blank(),
        axis.title.x = element_blank(),axis.ticks.y = element_blank(),axis.ticks.x=element_line(size=1,color="black"),
        axis.title.y = element_blank(),axis.line.y = element_blank(),axis.line.x=element_line(size=1,color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),panel.spacing = unit(1, "lines"))+
  facet_grid(. ~ timeframe)+scale_x_continuous(limits=c(-2,4))+
  scale_fill_manual(values=c("#4a536b","#aed6dc","#408ec6","#d72631","#f3ca20","#77c593",
                             "#B85042","#ff6e40","#500472"))

fig_2B_plot

ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/figures/Fig2b_Acute_effects.png", fig_2B_plot, width=4, height=5, bg = "transparent")



####################################
####### Maaslin2 persistent ########
####################################

library("Maaslin2")

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/10_persistent_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/10_control_metaphlan_community.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.02,
  output = "./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/10_persistent_02", 
  fixed_effects = c("X10_penicillin","X10_1st_gen_cephalosporin","X10_3rd_gen_cephalosporin","X10_carbapenem",
                    "X10_lincosamide","X10_amingoglycoside","X10_VitaminA","X10_MV_iron","X10_Lasix",
                    "X10_Famotidine","X10_Dopamine","X10_Caffeine",
                    "DOL_log","PotassiumChloride_recent","Phenobarbitol_recent","Morphine_recent","Midazolam_recent",
                    "Lasix_recent","Iron_recent","Insulin_recent","Famotidine_recent",
                    "Dopamine_recent","Dexamethasone_recent","Caffeine_recent","Acetaminophen_recent",
                    "VitaminA_recent","Recent_3rd_gen_cephalosporin","Recent_4th_gen_cephalosporin",
                    "Recent_glycopeptide","GA","Recent_1st_gen_cephalosporin","Recent_2nd_gen_cephalosporin",
                    "Recent_amingoglycoside","Recent_carbapenem","Recent_lincosamide","Recent_penicillin","Recent_macrolide",
                    "Other_ANTIM_recent","formula_cumulative","milk_cumulative"
  ),
  random_effects = c("Patients","Study"))

##################################################
###### Bifidobacterium_breve/VitaminA/DOL10 ######
##################################################

SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/10_Bifidobacterium_breve.txt", header = TRUE,row.names = 1)

B_breve<-ggplot(SM_DF, aes(x=DOL, y=Bifidobacterium_breve,col=X10_VitaminA))+
  geom_point(size=1, pch=20,alpha=0)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.border = element_blank(),panel.grid.major.y = element_line(size=0.25,color="#818181"),panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks =  element_line(size=0.75),
        axis.title.y = element_blank(),axis.line =  element_line(size=0.75),
        panel.background = element_blank())+
  ylab('Bifidobacterium breve')+
  xlab('DOL')+geom_smooth(method="loess")+
  scale_color_manual(values = c("#f57e7e","#1e847f"))+scale_y_continuous(breaks=c(0,40,80))
B_breve
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/10_B_breve.png", B_breve, width=4.5, height=2, bg = "transparent")



SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/10_Klebsiella_pneumoniae.txt", header = TRUE,row.names = 1)

K_pneumo<-ggplot(SM_DF, aes(x=DOL, y=Klebsiella_pneumoniae,col=X10_Caffeine))+
  geom_point(size=1, pch=20,alpha=0)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.border = element_blank(),panel.grid.major.y = element_line(size=0.25,color="#818181"),panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks.y =  element_line(size=0.75),axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),axis.line.y =  element_line(size=0.75),
        panel.background = element_blank())+
  ylab('Klebsiella_pneumoniae')+
  xlab('DOL')+geom_smooth(method="loess")+
  scale_color_manual(values = c("#f57e7e","#1e847f"))+scale_y_continuous(breaks=c(0,40,80),limits=c(-5,80))
K_pneumo
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/10_K_pneumo.png", K_pneumo, width=4.5, height=2, bg = "transparent")


#########################
### DOL 30 persistent ###
#########################

library("Maaslin2")

SM_meta <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/30_persistent_metadata.txt", header = TRUE,row.names = 1, na.strings = "NA")

#####Species
SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/30_control_metaphlan_community.txt", header = TRUE,row.names = 1)

fit_data = Maaslin2(
  input_data = SM_DF, 
  input_metadata = SM_meta,
  transform = "LOG",
  normalization = "NONE",
  min_abundance = 0.0001,
  min_prevalence = 0.05,
  output = "./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/30_persistent_05", 
  fixed_effects = c("DOL_log","Recent_4th_gen_cephalosporin","Recent_glycopeptide","GA","Other_ANTIM_recent","milk_cumulative",
                    "formula_cumulative","Morphine_recent","Midazolam_recent","Lasix_recent","Acetaminophen_recent",
                    "Recent_1st_gen_cephalosporin","Recent_macrolide","Recent_lincosamide",
                    "X30_Dexamethasone","X30_Iron","X30_Insulin","X30_Famotidine","X30_Midazolam","X30_Potassium_Chloride",
                    "X30_VitaminA","X30_Formula","X30_4thGenCep","X30_Carbapenem","X30_Lincosamide","X30_Glycopeptide"
  ),
  random_effects = c("Patients","Study"))


SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/30_Eubacterium_rectale.txt", header = TRUE,row.names = 1)

Eubacterium_rectale<-ggplot(SM_DF, aes(x=DOL, y=Eubacterium_rectale,col=X30_VitaminA))+
  geom_point(size=1, pch=20,alpha=0)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.border = element_blank(),panel.grid.major.y = element_line(size=0.25,color="#818181"),panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks =  element_line(size=0.75),
        axis.title.y = element_blank(),axis.line =  element_line(size=0.75),
        panel.background = element_blank())+
  ylab('Erysipelatoclostridium_ramosum')+
  xlab('DOL')+geom_smooth(method="loess")+
  scale_color_manual(values = c("#f57e7e","#1e847f"))+scale_y_continuous(breaks=c(0,20,40),limits=c(-5,40),oob=rescale_none)+
  scale_x_continuous(breaks=c(40,60,80),limits=c(30,80))
Eubacterium_rectale
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/30_Eubacterium_rectale.png", Eubacterium_rectale, width=4.5, height=2, bg = "transparent")


SM_DF <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/30_Enterococcus_faecalis.txt", header = TRUE,row.names = 1)

E_feacalis<-ggplot(SM_DF, aes(x=DOL, y=Enterococcus_faecalis,col=X30_Glycopeptide))+
  geom_point(size=1, pch=20,alpha=0)+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.border = element_blank(),panel.grid.major.y = element_line(size=0.25,color="#818181"),panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x = element_blank(),axis.ticks.y =  element_line(size=0.75),axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),axis.line.y =  element_line(size=0.75),
        panel.background = element_blank())+
  ylab('Enterococcus_faecalis')+
  xlab('DOL')+geom_smooth(method="loess")+
  scale_color_manual(values = c("#f57e7e","#1e847f"))+scale_y_continuous(breaks=c(0,50,100),limits=c(-5,100),oob=rescale_none)+
  scale_x_continuous(breaks=c(40,60,80),limits=c(30,80))
E_feacalis
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/30_E_feacalis.png", E_feacalis, width=4.5, height=2, bg = "transparent")

##########################
##### Summary Figure #####
##########################

library(tidyverse)

#mydat <- structure(list(Species = structure(c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L), .Label = c("setosa", "versicolor", "virginica"), class = "factor"), flower_att = c("Sepal.Length", "Sepal.Length", "Sepal.Length", "Sepal.Width", "Sepal.Width", "Sepal.Width", "Petal.Length", "Petal.Length", "Petal.Length", "Petal.Width", "Petal.Width", "Petal.Width"), measurement = c(5.1, 7, 6.3, 3.5, 3.2, 3.3, 1.4, 4.7, 6, 0.2, 1.4, 2.5), month = c("January", "February", "January", "February", "January", "February", "January", "February", "January", "February", "January", "February")),
#                   row.names = c(NA, -12L), class = "data.frame"
#)

make_triangles <- function(x, y, point = "up") {
  x <- as.integer(as.factor((x)))
  y <- as.integer(as.factor((y)))
  
  if (point == "up") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x - 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y + 0.5, y + 0.5)
    }, simplify = FALSE)
  } else if (point == "down") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x + 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y - 0.5, y + 0.5)
    }, simplify = FALSE)
  }
  data.frame(x = unlist(newx), y = unlist(newy))
}

#####################
######## 10 #########
#####################

mydat <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/10_persistent_plot_2C.txt", header = TRUE)
month.name<-c("Staphylococcus","Klebsiella","Streptococcus","Proteus","Lactobacillus","Clostridium","Actinomyces","Bacteroides","Bifidobacterium")

# required, otherwise you cannot use the values as fill
mydat_wide <- mydat %>% pivot_wider(names_from = "flower_att", values_from = "measurement")
# making your ordered months factor
mydat_wide$month <- droplevels(factor(mydat_wide$month, levels = month.name))
# The actual triangle computation
newcoord_up <- make_triangles(mydat_wide$month, mydat_wide$Species)
newcoord_down <- make_triangles(mydat_wide$month, mydat_wide$Species, point = "down")
# just a dirty trick for renaming
newcoord_down <- newcoord_down %>% select(xdown = x, ydown = y)
# you need to repeat each row of your previous data frame 3 times
repdata <- map_df(1:nrow(mydat_wide), function(i) mydat_wide[rep(i, 3), ])
newdata <- bind_cols(repdata, newcoord_up, newcoord_down)
newdata[is.na(newdata)] <- 0

test<-ggplot(newdata) +
  geom_polygon(aes(x = x, y = y, fill = enriched.real, group = interaction(Species, month)), color = "black") +
  scale_fill_gradientn(colors=c("#FFFFFF","#1e847f","#1e3d59"), limits = c(0, 5)) +
  ggnewscale::new_scale_fill() +
  geom_polygon(aes(x = xdown, y = ydown, fill = depleted.real, group = interaction(Species, month)), color = "black") +
  scale_fill_gradientn(colors = c("#FFFFFF", "#d2601a","#d72631"), limits = c(0, 5)) +
  scale_x_continuous(breaks = seq_along(unique(mydat_wide$month)), 
                     labels = unique(levels(mydat_wide$month))) +
  scale_y_continuous(breaks = seq_along(unique(mydat_wide$Species)),
                     labels = unique(mydat_wide$Species))+
  coord_equal()+theme_classic()
test
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/10_summary.png", test, width=9, height=8, bg = "transparent")

#####################
######## 30 #########
#####################

mydat <- read.delim("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/30_persistent_plot_2F.txt", header = TRUE)

# required, otherwise you cannot use the values as fill
mydat_wide <- mydat %>% pivot_wider(names_from = "flower_att", values_from = "measurement")
# making your ordered months factor
mydat_wide$month <- droplevels(factor(mydat_wide$month, levels = month.name))
# The actual triangle computation
newcoord_up <- make_triangles(mydat_wide$month, mydat_wide$Species)
newcoord_down <- make_triangles(mydat_wide$month, mydat_wide$Species, point = "down")
# just a dirty trick for renaming
newcoord_down <- newcoord_down %>% select(xdown = x, ydown = y)
# you need to repeat each row of your previous data frame 3 times
repdata <- map_df(1:nrow(mydat_wide), function(i) mydat_wide[rep(i, 3), ])
newdata <- bind_cols(repdata, newcoord_up, newcoord_down)
newdata[is.na(newdata)] <- 0

test<-ggplot(newdata) +
  geom_polygon(aes(x = x, y = y, fill = enriched.real, group = interaction(Species, month)), color = "black") +
  scale_fill_gradientn(colors=c("#FFFFFF","#1e847f","#1e3d59"), limits = c(0, 5)) +
  ggnewscale::new_scale_fill() +
  geom_polygon(aes(x = xdown, y = ydown, fill = depleted.real, group = interaction(Species, month)), color = "black") +
  scale_fill_gradientn(colors = c("#FFFFFF", "#d2601a","#d72631"), limits = c(0, 5)) +
  scale_x_continuous(breaks = seq_along(unique(mydat_wide$month)), 
                     labels = unique(levels(mydat_wide$month))) +
  scale_y_continuous(breaks = seq_along(unique(mydat_wide$Species)),
                     labels = unique(mydat_wide$Species))+
  coord_equal()+theme_classic()
test
ggsave("./Desktop/Projects/NEC/final_analysis/control_analysis/DNA_Maaslin2/DNA_Maaslin2_taxa/persistent_effects/30_summary.png", test, width=9, height=8, bg = "transparent")
























