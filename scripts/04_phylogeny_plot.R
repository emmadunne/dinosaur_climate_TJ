# ******************************************************
#
#   Early Dinosaurs & Climate
#
#   E. M. Dunne (2022)
# ______________________________________________________
#
#   Phylogeny plotting using phytools
# 
# ******************************************************

## Load packages:
library(tidyverse)
library(phytools)
library(viridis)
library(strap)
library(paleotree)



## Load time-scaled trees that were cleaned in script 03:
trees_Chapelle <- read.nexus("./trees/trees_Chapelle_cleaned.nex")
#trees_McPhee <- read.nexus("./trees/trees_McPhee_cleaned.nex")

## Randomly pick a single tree for plotting:
random_tree <- sample(1:100, 1)
tree <- trees_Chapelle[[random_tree]]


## Organise the continuous data:

## Import mean MAT from script 03:
mean_MAT_values <- read.csv("./data/climate_data_table.csv", row.names = 1)

## remove taxa that are in this dataset but not on the tree:
taxa_to_remove <- rownames(mean_MAT_values)[ !rownames(mean_MAT_values) %in% tree$tip.label ] # in MAT data but not on tree
mean_MAT_values <- mean_MAT_values[!rownames(mean_MAT_values) %in% taxa_to_remove , ]
rownames(mean_MAT_values)[ !rownames(mean_MAT_values) %in% tree$tip.label ] #check

## Convert to matrix:
MAT_matrix <- as.matrix(mean_MAT_values) [,2] # 2nd column

## contMap()
MATmapped <- contMap(tree, MAT_matrix, plot = FALSE)
MATmapped <- setMap(MATmapped, invert = TRUE)
n <- length(MATmapped$cols)
MATmapped$cols[1:n] <- plasma(n)
plot(MATmapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAT (°C)")


#### Plot tree on a time scale:
tipAges <- cbind(c("Guaibasaurus_candelariensis", "Emausaurus_ernsti", "Lewisuchus_admixtus"), c(224.7, 182.5, 235.35)) # constrained taxa plus their mean stat ages
tree2 <- setRootAge(tree, tipAges)
geoscalePhylo(tree = tree2, boxes = "Age", cex.tip = 0.4)
