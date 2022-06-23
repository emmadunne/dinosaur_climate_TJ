# ******************************************************
#
#   Early Dinosaurs & Climate
#
#   E. M. Dunne & P. L. Godoy (2022)
# ______________________________________________________
#
#   Preparing tree data to be time-scaled using 
#     the FBD process in MrBayes
# 
# ******************************************************

## Load packages
library(paleotree)
library(phytools)
library(geiger)



## ______________________________________________________
## 1. Creating NEXUS file for MrBayes -------------------
## ______________________________________________________


## Importing the files

# Tree
treeConstraints <- read.nexus("./trees/dino_supertree_Chapelle.nex") # Switch to other trees where necessary
str(treeConstraints)
treeConstraints$tip.label

# Ages
tipTimes_full <- read.csv("./trees/dino_ages.csv", header=T) # import full dataset
tipTimes <- subset(tipTimes_full, select=c(Taxon, FAD, LAD)) # truncate to only necessary columns
tipTimes$Taxon <- gsub(" ", "_", tipTimes$Taxon) # replace spaces with underscores to match the tree
rownames(tipTimes) <- tipTimes$Taxon # rownames as taxa names
tipTimes <- tipTimes[,-1] # remove Taxon column
str(tipTimes) # check if everything looks right
rownames(tipTimes)


## Check that the tree and the age data match
name.check(treeConstraints, tipTimes)

## Drop tips from tree that have no data
treeConstraints2 <- drop.tip(treeConstraints, treeConstraints$tip.label[!(treeConstraints$tip.label %in% rownames(tipTimes))])
treeConstraints <- treeConstraints2

## Drop taxa from data that aren't on tree
tipTimes2 <- tipTimes[(rownames(tipTimes) %in% treeConstraints$tip.label), ]
str(tipTimes2)
tipTimes <- tipTimes2

## Check again:
name.check(treeConstraints, tipTimes) # yay! :)


## Save results (the command file) as:
file <- "./trees/dino_tree_Chapelle_FBD_MrB.nex"


## Creating the NEXUS file for MrBayes

## *** Below, be sure to change: ***
## anchorTaxon = one with 'surest' dates, doesn't matter what age, but best to be percise
## runName = dino_FBD
## ngen = 30 million 
createMrBayesTipDatingNexus(tipTimes, outgroupTaxa = NULL, treeConstraints = treeConstraints, ageCalibrationType = "uniformRange",
                            whichAppearance = "first", treeAgeOffset = 0.1, minTreeAge = NULL, collapseUniform = TRUE,
                            anchorTaxon = "Leyesaurus_marayensis", newFile = file, origNexusFile = NULL, parseOriginalNexus = TRUE,
                            createEmptyMorphMat = TRUE, morphModel = "strong", runName = "dino_FBD", ngen = "15000000",
                            doNotRun = FALSE, cleanNames = TRUE, printExecute = TRUE)

# Should get this message in the console:
#     Now go to MrBayes and paste in this line:  
#     Execute "/Users/emmadunne/Projects/Early dinosaurs/trees/dino_tree_FBD_MrB.nex";



## ______________________________________________________
## 2. Edit NEXUS file in TextEdit -----------------------
## ______________________________________________________

#   Open the .nex file that was just produced (stored in the same older as specified above)
#     a. Scroll down towards the bottom
#       - find prset treeagepr - change offsetexp to uniform
#       - These ages should correspond to what you think should be the 
#          oldest age of the root of the entire tree (out- and ingroup)
#     b. Go to prset sampleprob (the proportion of extant taxa)
#       - Set to 0.2
#     c. Go to burninfrac
#       - Set to 0.25
#     d. Add sumt; to the final section [RUN] after sump;
#       - This gives you out a summary tree
#     e. Save file & close




## ______________________________________________________
## 3. Run MrBayes analysis in the terminal --------------
## ______________________________________________________

#     a. Set the working directory - I copied the .nex file from above into a new folder on the desktop for ease
#       - cd changes the directory 
#       - e.g. cd ./Desktop/MrBayes_
#     b. Check the directory using pwd (print working directory)
#     c. Activate MrBayes - type mb
#     d. Run the analysis using:
#       - Execute dino_tree_Chapelle_FBD_MrB.nex
#     e. Watch the “average standard number”
#      - should be below 0.01




## ______________________________________________________
## 4. Getting trees from FBD analyses -------------------
## ______________________________________________________


# set work directory
setwd("~/Desktop/MrBayes_Chapelle")


## You'll need to make sure the files from the MrBayes analysis are in the same folder
#   to be more specific, you'll need a total of 5 files:
#   first, you need the .nex file you used to run the analysis in MrBayes
#   then, you'll also need 4 other files: two .p files and two .t files (e.g. "dino_tree_FBD.run1.p")
# **** make sure all these files have the same name as the stem for the .nex file (e.g. "dino_tree_Chapelle_FBD_MrB") ****


### Getting single MCC tree
dino_MCCT_tree <- obtainDatedPosteriorTreesMrB(runFile = "dino_tree_Chapelle_FBD_MrB.run1.t", 
                                               nRuns = 2, burnin = 0.25, outputTrees = "MCCT",
                                               labelPostProb = FALSE, getFixedTimes = TRUE,
                                               originalNexusFile = NULL, file = NULL) 
#this can take some time (~15 minutes?) - depends on the number of generations used


# great, let's check if it looks right
dino_MCCT_tree$root.time # this checks if the tree has an element for the age of the root

# if it shows a number, that's good
# but let's be extra cautious and plot the tree 
dino_MCCT_tree <- ladderize(dino_MCCT_tree,right=FALSE)
plot(dino_MCCT_tree,direction="right",cex=0.3)
axisPhylo()
# if it looks right (the right ages), you can save the tree


# saving files
# there are different formats you can use to save the files
# I recommend using two of those: .nex and .txt
# the .nex is important if you want to open it in Mesquite or other software (such as FigTree) and for sharing the tree(s) as a supplementary file
# but the .txt file is the best for importing the tree(s) into R again for running extra analyses, since it keeps the information about the age of the root ($root.time)

write.nexus(dino_MCCT_tree, file="dino_Chapelle_MCCT_tree.nex")
dput(dino_MCCT_tree, file = "dino_Chapelle_MCCT_tree.txt") #for keeping the root age


### Getting 100 trees
dino_100_trees <- obtainDatedPosteriorTreesMrB(runFile = "dino_tree_Chapelle_FBD_MrB.run1.t", 
                                               nRuns = 2, burnin = 0.25, outputTrees = 100,
                                               labelPostProb = FALSE, getFixedTimes = TRUE,
                                               originalNexusFile = NULL, file = NULL)

## Check tree number 1:
dino_100_trees[[1]]$root.time

## Ladderizing and plotting:
for (i in 1: length(dino_100_trees)) {
  dino_100_trees[[i]]<- ladderize(dino_100_trees[[i]],right=FALSE)} 
plot(dino_100_trees[[73]],direction="right",cex=0.3)
axisPhylo()

## Saving:
write.nexus(dino_100_trees, file="dino_Chapelle_100_trees.nex")
dput(dino_100_trees, file = "dino_Chapelle_100_trees.txt")  #for keeping the root age
