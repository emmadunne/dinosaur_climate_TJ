# ******************************************************
#
#   Early Dinosaurs & Climate
#
#   E. M. Dunne (2022)
# ______________________________________________________
#
#   Evolutionary model fitting analysis
# 
# ******************************************************


## Load packages:
library(tidyverse)
library(reshape)
library(ape)
library(OUwie)
library(phytools)
library(beepr)
library(cowplot)
library(ggpubr)
library(RColorBrewer)



# Prepare climate data for species on supertree ----------------------------------------------------

## Load dataset where climate variables have been assigned to each species on the supertree using the coordinates of their occurrences
tree_climate <- read_csv("./data/supertree_occs_climate.csv") 

## Extract the standard deviations of climate variables
sd.MAT.temp <- aggregate( tree_climate$MAT_degC , list( tree_climate$taxon_name) , sd )
sd.MAP.temp <- aggregate( tree_climate$MAP_cm , list( tree_climate$taxon_name) , sd )
sd.s_temp.temp <- aggregate( tree_climate$s_temp , list( tree_climate$taxon_name) , sd )
sd.s_precip.temp <- aggregate( tree_climate$s_precip , list( tree_climate$taxon_name) , sd )

## Extract the means of climate variables
mean.MAT.temp <- aggregate( tree_climate$MAT_degC , list( tree_climate$taxon_name) , mean )
mean.MAP.temp <- aggregate( tree_climate$MAP_cm , list( tree_climate$taxon_name) , mean )
mean.s_temp.temp <- aggregate( tree_climate$s_temp , list( tree_climate$taxon_name) , mean )
mean.s_precip.temp <- aggregate( tree_climate$s_precip , list( tree_climate$taxon_name) , mean )

occurrences.temp <- aggregate( tree_climate$MAT_degC , list( tree_climate$taxon_name) , length )

## Combine into single dataset: (previously climate.data.table.temp)
climate_table <- cbind( sd.MAT.temp , mean.MAT.temp[,2] , sd.MAP.temp[,2] ,
                        mean.MAP.temp[,2] , sd.s_temp.temp[,2] , mean.s_temp.temp[,2] ,
                        sd.s_precip.temp[,2] , mean.s_precip.temp[,2] ,
                        occurrences.temp[,2] )[ order(occurrences.temp[,2] , decreasing = TRUE ),]

## Change col and row names
colnames( climate_table ) <- c( "taxon_name" , "sd.MAT" , "mean.MAT" , "sd.MAP" , "mean.MAP" , "sd.s_temp" , "mean.s_temp" , "sd.s_precip" , "mean.s_precip" , "occurrence.count" )
rownames( climate_table ) <- climate_table[,"taxon_name"]


## Save a copy of this data for ease of plotting in script 04:
write_csv(climate_table , file = "./data/climate_data_table.csv")




# Clean tree data ---------------------------------------------------------------


## Load time-scaled trees (100 randomly resolved topologies for two alternative sauropodomorph trees)
## (NB: The rest of the code in this script follows the Chapelle sauropodomorph topology for ease)
dino_Chapelle_100 <- read.tree("./trees/100_ts_trees_Chapelle.tre")
# dino_McPhee_100 <- read.tree(file = "./trees/100_ts_trees_McPhee.tre") 

## List of taxa to remove from supertree
taxa_to_remove <- c("Arizonasaurus", "Ctenosauriscus", "Dongusuchus", 
                    "Teleocrater", "Xilousuchus", "Yarasuchus", # Outgroup taxa
                    "Adeopapposaurus_mognai", # Early Jurassic
                    "Asilisaurus_kongwe", "Lutungutali_sitwensis") # Anisian

## Remove from tree
trees_Chapelle <- lapply(dino_Chapelle_100, drop.tip, tip = taxa_to_remove)
## ...and return objects to correct class
class(trees_Chapelle)<-"multiPhylo"

## Check for differences between tree and climate data:
dtaxa.Chapelle <- trees_Chapelle[[1]]$tip.label[ !trees_Chapelle[[1]]$tip.label %in% rownames(climate_table) ] # on Chapelle supertree but not in temp. data
dtaxa.Chapelle

## Drop the taxa that have missing data
trees_Chapelle <- lapply( trees_Chapelle , drop.tip , tip = dtaxa.Chapelle )
## ...and return object to correct class
class(trees_Chapelle) <- "multiPhylo"


## Save ccleaned copies for plotting in script 04
writeNexus(trees_Chapelle, "./trees/trees_Chapelle_cleaned.nex")
#writeNexus(trees_McPhee, "./datasets/evo_rates/trees_McPhee_cleaned.nex")



# Model set-up ------------------------------------------------------------

## Find the tip height of Laquintasaurus (because it comes immediately after the T/J)
Laquintasaurus_heights <- list()
for (i in 1:length(trees_Chapelle)) {
  Laquintasaurus_heights[[i]] <- vcv.phylo(trees_Chapelle[[i]])[ "Laquintasaura_venezuelae" , "Laquintasaura_venezuelae" ]
}
	
## Make an 'era' simmap tree (i.e. a tree with regimes painted onto the branches according to time)
era_trees <- list()
for (i in 1:length(trees_Chapelle)) {
  era_trees[[i]] <- make.era.map(trees_Chapelle[[i]] , c( 0 , Laquintasaurus_heights[[i]] - 1 )) # age of Laquintasaurus minus 1.0 Ma
}
class(era_trees) <- "multiPhylo" # return it to the correct class
plotSimmap(era_trees[[1]]) # see what it looks like

## Subset the tree to a specific group of interest. In this example, Sauropodomorpha
clade.specifiers <- c( "Saturnalia_tupiniquim" , "Isanosaurus_attavipachi" )	# MRCA of these two defines Sauropodomorpha

subclade_era_trees <- list()
for (i in 1:length(era_trees)) {
  subclade_era_trees[[i]] <- extract.clade.simmap(era_trees[[i]] , getMRCA(era_trees[[i]] , clade.specifiers))
}
class(subclade_era_trees) <- "multiPhylo" # return to the correct class
plotSimmap(subclade_era_trees[[1]])	# see what it looks like


## The other part of the tree i.e. non-sauropodomorph dinosaurs:
antisubclade_era_trees <- lapply(era_trees, drop.tip.simmap , tip = subclade_era_trees[[1]]$tip.label)
plotSimmap(antisubclade_era_trees[[1]])	# see what it looks like



# Sauropodomorpha (subclade) - run models ---------------------------------------------


## Subset the data for OUwie: 
data_subclades <- list()
for (i in 1:length(subclade_era_trees)) {
  data_subclades[[i]] <- data.frame(
    taxon = subclade_era_trees[[i]]$tip.label ,
    regime = names( unlist( subclade_era_trees[[i]]$maps ) )[ 1:length( subclade_era_trees[[i]]$tip.label ) ] ,
    trait = climate_table[ subclade_era_trees[[i]]$tip.label , "mean.MAT" ] ,
    stdev = rep( 3.1652468, length( subclade_era_trees[[i]]$tip.label ) ) 
  )
}

## Set the vector of models that you want to fit
models = c( "BM1" , "BMS" , "OU1" , "OUM" , "OUMV" , "OUMA" , "OUMVA" )

## Run model fitting analysis
OUwie_results_subclade <- list()
for(i in 1:length(subclade_era_trees)){
  OUwie_results_subclade[[i]] <- lapply(models, OUwie,
                                        phy = subclade_era_trees[[i]],
                                        data = data_subclades[[i]],
                                        simmap.tree = TRUE,
                                        mserr = "known" )
  names( OUwie_results_subclade[[i]] ) <- models
}
beep("coin") # to indicate when run is complete :)



# Sauropodomorpha (subclade) - model output -------------------------------

### Extract AICc scores
AICc_scores_subclade <- list()
for(i in 1:length(OUwie_results_subclade)) {
  AICc_scores_subclade[[i]] <- unlist( lapply( OUwie_results_subclade[[i]] , function( X ){ X$AICc } ) )
}
AICc_subclade_df <- as.data.frame(do.call("rbind", AICc_scores_subclade))

### Plot in ggplot
AICc_sub_melt <- melt(AICc_subclade_df) # 'melt' dataset so it can be plotted in ggplot
AICc_sub_melt$variable <- factor(AICc_sub_melt$variable, levels = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), ordered = TRUE) # arrange variables

## Colour palette:
orange_gradient <- colorRampPalette( c("#DA7209","#FFA402", "#FEECC4") )

AICc_sub_boxplot <- ggplot(AICc_sub_melt, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  scale_fill_manual(values = orange_gradient(7)) +
  xlab("") + ylab("AICc") + ggtitle("Sauropodomorpha") +
  theme_bw(base_size = 12) + theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) 
AICc_sub_boxplot


### AICc weights 
AICc_weights_subclade <- list()
for(i in 1:length(AICc_scores_subclade)) {
  AICc_weights_subclade[[i]] <- round( exp(-0.5 * (AICc_scores_subclade[[i]] - min( AICc_scores_subclade[[i]] )) ) / sum(exp(-0.5 * (AICc_scores_subclade[[i]] - min( AICc_scores_subclade[[i]] )) )) , 2 )
}
AICc_weights_subclade_df <- as.data.frame(do.call("rbind", AICc_weights_subclade))

### Plot in ggplot
AICc_weights_sub_melt <- melt(AICc_weights_subclade_df)
AICc_weights_sub_melt$variable <- factor(AICc_weights_sub_melt$variable, levels = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), ordered = TRUE)

AICc_w_sub_boxplot <- ggplot(AICc_weights_sub_melt, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  scale_fill_manual(values = orange_gradient(7)) +
  xlab("") + ylab("AICc weight") + ggtitle("Sauropodomorpha") +
  theme_bw(base_size = 12) + theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  scale_y_continuous(limits = c(0, 0.55))
AICc_w_sub_boxplot


### Alpha scores:
## (Similar code can be used to extract theta values)
alpha_subclade <- list()
for(i in 1:length(OUwie_results_subclade)) {
  alpha_subclade[[i]] <- unlist( lapply( OUwie_results_subclade[[i]] , function( X ){ X$solution[1,1] } ) )
}
alpha_subclade_df <- as.data.frame(do.call("rbind", alpha_subclade))
alpha_subclade_df



# Other dinosaurs (anti-subclade) - run models --------------------------

## Next is the same thing for the tree -excluding- the clade of interest

## Subset the data for OUwie: 
data_antisubclades <- list()
for (i in 1:length(antisubclade_era_trees)) {
  data_antisubclades[[i]] <- data.frame(
    taxon = antisubclade_era_trees[[i]]$tip.label ,
    regime = names( unlist( antisubclade_era_trees[[i]]$maps ) )[ 1:length( antisubclade_era_trees[[i]]$tip.label ) ] ,
    trait = climate_table[ antisubclade_era_trees[[i]]$tip.label , "mean.MAT" ] ,
    stdev = rep( 3.1652468 , length( antisubclade_era_trees[[i]]$tip.label ) )
  )
}


## Set the vector of models that you want to fit
models = c( "BM1" , "BMS" , "OU1" , "OUM" , "OUMV" , "OUMA" , "OUMVA" )

## Run model fitting analysis
OUwie_results_antisubclade <- list() 
for(i in 1:length(antisubclade_era_trees)){
  OUwie_results_antisubclade[[i]] <- lapply(models, OUwie, 
                                            phy = antisubclade_era_trees[[i]], 
                                            data = data_antisubclades[[i]], 
                                            simmap.tree = TRUE, 
                                            mserr = "known" )
  names( OUwie_results_antisubclade[[i]] ) <- models
}
beep("coin")


# Other dinosaurs (anti-subclade) - model output --------------------------

## Extract AICc, alpha, theta values, etc. as for subclade above

### Extract AICc scores 
AICc_scores_antisubclade <- list()
for(i in 1:length(OUwie_results_antisubclade)) {
  AICc_scores_antisubclade[[i]] <- unlist( lapply( OUwie_results_antisubclade[[i]] , function( X ){ X$AICc } ) )
}
AICc_antisubclade_df <- as.data.frame(do.call("rbind", AICc_scores_antisubclade))


### Plot in ggplot
AICc_anti_melt <- melt(AICc_antisubclade_df)
AICc_anti_melt$variable <- factor(AICc_anti_melt$variable, levels = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), ordered = TRUE)

## Colour palette:
teal_gradient <- colorRampPalette( c("#025F87","#0094A6", "#D8FAFE") )

AICc_anti_boxplot <- ggplot(AICc_anti_melt, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  scale_fill_manual(values = teal_gradient(7)) +
  xlab("") + ylab("AICc") + ggtitle("Non-sauropodomorpha") +
  theme_bw(base_size = 12) + theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  scale_y_continuous(limits = c(435, 500))
AICc_anti_boxplot


#### AICc weights (more useful)
AICc_weights_antisubclade <- list()
for(i in 1:length(AICc_scores_antisubclade)) {
  AICc_weights_antisubclade[[i]] <- round( exp(-0.5 * (AICc_scores_antisubclade[[i]] - min( AICc_scores_antisubclade[[i]] )) ) / sum(exp(-0.5 * (AICc_scores_antisubclade[[i]] - min( AICc_scores_antisubclade[[i]] )) )) , 2 )
}
AICc_weights_antisubclade_df <- as.data.frame(do.call("rbind", AICc_weights_antisubclade))

### Plot in ggplot
AICc_weights_anti_melt <- melt(AICc_weights_antisubclade_df)
AICc_weights_anti_melt$variable <- factor(AICc_weights_anti_melt$variable, levels = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), ordered = TRUE)

AICc_w_anti_boxplot <- ggplot(AICc_weights_anti_melt, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  scale_fill_manual(values = teal_gradient(7)) +
  xlab("") + ylab("AICc weight") + ggtitle("Non-sauropodomorpha") +
  theme_bw(base_size = 12) + theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  scale_y_continuous(limits = c(0, 1.0))
AICc_w_anti_boxplot
