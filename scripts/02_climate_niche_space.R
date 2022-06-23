# ******************************************************
#
#   Early Dinosaurs & Climate
#
#   E. M. Dunne (2022)
# ______________________________________________________
#
#   Part 1:
#     Running the PCA analysis to illustrate 'climate 
#     niche space' for early dinosaurs (Figure 1)
#
#   Part 2:
#     Constructing raincloud plots to explore climate 
#     ranges (Figure 2 and Figure S1)
# 
# ******************************************************

## Load package(s):
library(tidyverse)
library(factoextra)
library(gridExtra)
library(ggalt)
library(ordr)
library(RVAideMemoire)
library(beepr)
library(vegan)
library(ggdist)


# PART 1: PCA ===============================================================================================

# Organise data -----------------------------------------------------------


## Load occurrence data where climate variables have been assigned to each occurrence based on the coordinates:
occs_plus_climate <- read_csv("./data/tetrapod_occs_climate.csv.csv")

## Truncate to just the data necessary for a PCA:
PCA_data <- select(occs_plus_climate, occurrence_no, MAT_degC, MAP_cm, s_temp, s_precip, group, stage) %>% distinct()
## Remove duplicated occs
PCA_data <- PCA_data[!duplicated(PCA_data$occurrence_no),]

## Divide into Late Triassic and Early Jurassic
PCA_data_LT <- filter(PCA_data, stage == "Carnian" | stage == "Norian" | stage == "Rhaetian") %>% na.omit()  %>% distinct()
PCA_data_EJ <- filter(PCA_data, stage == "Hettangian" | stage == "Sinemurian" | 
                        stage == "Pliensbachian" | stage == "Toarcian") %>% na.omit()  %>% distinct()

## Change the group names to be shorter for plotting
PCA_data_LT$group <- gsub("Dinosauromorpha", "Dino.", PCA_data_LT$group )
PCA_data_LT$group <- gsub("Sauropodomorpha", "Sauro.", PCA_data_LT$group )
PCA_data_EJ$group <- gsub("Dinosauromorpha", "Dino.", PCA_data_EJ$group )
PCA_data_EJ$group <- gsub("Sauropodomorpha", "Sauro.", PCA_data_EJ$group )



# PCA ---------------------------------------------------------------------

## Biplot confidence ellipses
## See here for more info: https://corybrunson.github.io/ordr/reference/stat-biplot-ellipse.html


### Late Triassic
LT_pca <- PCA_data_LT[, 2:5] %>%
  prcomp(center = TRUE, scale. = TRUE) %>% 
  as_tbl_ord() %>%
  mutate_rows(group = PCA_data_LT$group)

confellip_LT <- LT_pca %>%
  ggbiplot(aes(color = group)) +
  theme_bw() + theme(legend.position = "none") +
  geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  scale_colour_manual(values = c("#0891A3", "#FFA93D", "#B00B69")) +
  scale_fill_manual(values = c("#0891A3", "#FFA93D", "#B00B69"))
confellip_LT

## Example code to explore the output:
loading_scores <- LT_pca$rotation[,1]
climate_scores <- abs(loading_scores) ## get the magnitudes
climate_score_ranked <- sort(climate_scores, decreasing=TRUE)
LT_pca$rotation[top_4,1] ## show the scores (and +/- sign)


### Early Jurassic:
EJ_pca <- PCA_data_EJ[, 2:5] %>%
  prcomp(scale = TRUE) %>%
  as_tbl_ord() %>%
  mutate_rows(group = PCA_data_EJ$group)

confellip_EJ <- EJ_pca %>%
  ggbiplot(aes(color = group)) +
  theme_bw() + theme(legend.position = "none") +
  geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  scale_colour_manual(values = c("#0891A3", "#FFA93D", "#B00B69")) +
  scale_fill_manual(values = c("#0891A3", "#FFA93D", "#B00B69"))
confellip_EJ



# npMANOVA ----------------------------------------------------------------


## More info here: https://www.rdocumentation.org/packages/RVAideMemoire/versions/0.9-80/topics/pairwise.perm.manova

# The "fact" argument will be your groups (a two-column data frame; first column with species names 
# and second column with the group they belong to [i.e., non-dinosaur tetrapods, non-sauropodomorph 
# dinosaurs, or sauropodomorphs]). 

LT_PCgroups <- PCA_data_LT$group
head(LT_PCgroups)

EJ_PCgroups <- PCA_data_EJ$group
head(EJ_PCgroups)

# The "resp" argument will be the euclidian distance of your PC scores (you can use all PCs instead 
# of just PC1 and PC2). Again, you'll need the PC scores as a data frame. First column being the taxa; 
# other columns the PCs.

LT_PCscores <- as.data.frame(LT_pca$x) # PC scores from above
LT_PCscores <- cbind(LT_PCscores, LT_PCgroups)

EJ_PCscores <- as.data.frame(EJ_pca$x) # PC scores from above
EJ_PCscores <- cbind(EJ_PCscores, EJ_PCgroups)


## npMANOVA:
npmanova_LT <- pairwise.perm.manova(dist(LT_PCscores[,1:2], "euclidian"), LT_PCscores$LT_PCgroups, nperm = 10000, p.method = "BH")
npmanova_EJ <- pairwise.perm.manova(dist(EJ_PCscores[,1:2], "euclidian"), EJ_PCscores$EJ_PCgroups, nperm = 10000, p.method = "BH")
beep("fanfare") # to indicate when above lines have run
# The results will be a table of p-values for pairwise comparisons between the groups. 
# If p is lower than 0.05, the groups are significantly different (i.e., the PC scores of them are significantly different)



# PART 2: Raincloud plots ===============================================================================================

## For more info see: https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/

## Taking the data imported above...
## Divide into Late Triassic and Early Jurassic and omit NAs
boxplot_data_LT <- occs_plus_climate %>% filter(stage == "Carnian" | stage == "Norian" | stage == "Rhaetian") %>% na.omit()
boxplot_data_EJ <- occs_plus_climate %>% filter(stage == "Hettangian" | stage == "Sinemurian" | 
                                                  stage == "Pliensbachian" | stage == "Toarcian") %>% na.omit()

## Change the group names to be shorter for plotting
boxplot_data_LT$group <- gsub("Dinosauromorpha", "Dino.", boxplot_data_LT$group )
boxplot_data_LT$group <- gsub("Sauropodomorpha", "Sauro.", boxplot_data_LT$group )
boxplot_data_EJ$group <- gsub("Dinosauromorpha", "Dino.", boxplot_data_EJ$group )
boxplot_data_EJ$group <- gsub("Sauropodomorpha", "Sauro.", boxplot_data_EJ$group )

## Arrange the group column before plotting
boxplot_data_LT$group <- factor(boxplot_data_LT$group,
                                levels = c("Tetrapoda", "Dino.", "Sauro."),ordered = TRUE)
boxplot_data_EJ$group <- factor(boxplot_data_EJ$group,
                                levels = c("Tetrapoda", "Dino.", "Sauro."),ordered = TRUE)


## Set theme:
theme_set(theme_minimal())
theme_update(
  panel.grid.major = element_line(color = "grey92", size = .4),
  panel.grid.minor = element_line(color = "grey92", size = .2),
  axis.title.x = element_text(color = "grey30", margin = margin(t = 7), size=12),
  axis.title.y = element_text(color = "grey30", margin = margin(r = 7), size=12),
  axis.text = element_text(color = "grey50", size = 11),
  plot.margin = margin(rep(15, 4)),
  legend.position = "none")


### Plots for each of the 4 variables in both the Late Triassic and Early Jurassic 

## Late Triassic
rain_LT_MAT <- ggplot(boxplot_data_LT, aes(x = group, y = MAT_degC, fill = group)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA, alpha = .6) + 
  geom_boxplot(aes(fill = group, colour = group), width = .25, outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(colour = group), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1)) + 
  scale_color_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  scale_fill_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  labs(x = NULL, y = "MAT (°C)") + ggtitle("Late Triassic") +
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

rain_LT_MAP <- ggplot(boxplot_data_LT, aes(x = group, y = MAP_cm, fill = group)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA, alpha = .6) + 
  geom_boxplot(aes(fill = group, colour = group), width = .25, outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(colour = group), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1)) + 
  scale_color_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  scale_fill_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  labs(x = NULL, y = "MAP (mm/day)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

rain_LT_ST <- ggplot(boxplot_data_LT, aes(x = group, y = s_temp, fill = group)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA, alpha = .6) + 
  geom_boxplot(aes(fill = group, colour = group), width = .25, outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(colour = group), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1)) + 
  scale_color_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  scale_fill_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  labs(x = NULL, y = "Seasonal range in temp. (°C)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

rain_LT_SP <- ggplot(boxplot_data_LT, aes(x = group, y = s_precip, fill = group)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA, alpha = .6) + 
  geom_boxplot(aes(fill = group, colour = group), width = .25, outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(colour = group), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1)) + 
  scale_color_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  scale_fill_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  labs(x = NULL, y = "Seasonal range in precip. (mm/day)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")



## Early Jurassic
rain_EJ_MAT <- ggplot(boxplot_data_EJ, aes(x = group, y = MAT_degC, fill = group)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA, alpha = .6) + 
  geom_boxplot(aes(fill = group, colour = group), width = .25, outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(colour = group), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1)) + 
  scale_color_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  scale_fill_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  labs(x = NULL, y = "MAT (°C)") + ggtitle("Early Jurassic") +
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

rain_EJ_MAP <- ggplot(boxplot_data_EJ, aes(x = group, y = MAP_cm, fill = group)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA, alpha = .6) + 
  geom_boxplot(aes(fill = group, colour = group), width = .25, outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(colour = group), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1)) + 
  scale_color_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  scale_fill_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  labs(x = NULL, y = "MAP (mm/day)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

rain_EJ_ST <- ggplot(boxplot_data_EJ, aes(x = group, y = s_temp, fill = group)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA, alpha = .6) + 
  geom_boxplot(aes(fill = group, colour = group), width = .25, outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(colour = group), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1)) + 
  scale_color_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  scale_fill_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  labs(x = NULL, y = "Seasonal range in temp. (°C)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

rain_EJ_SP <- ggplot(boxplot_data_EJ, aes(x = group, y = s_precip, fill = group)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA, alpha = .6) + 
  geom_boxplot(aes(fill = group, colour = group), width = .25, outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(colour = group), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1)) + 
  scale_color_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  scale_fill_manual(values = c("#B00B69", "#0891A3", "#FFA93D")) +
  labs(x = NULL, y = "Seasonal range in precip. (mm/day)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")


## All 4 climate variables:
rainplots <- ggarrange(rain_LT_MAT, rain_EJ_MAT,
                       rain_LT_MAP, rain_EJ_MAP,
                       rain_LT_ST, rain_EJ_ST,
                       rain_LT_SP, rain_EJ_SP,
                       ncol = 2, nrow = 4,
                       labels = c("(a)", "(e)", "(b)", "(f)", "(c)", "(g)", "(d)", "(h)"))
rainplots
