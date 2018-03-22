### Plotting PCA
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last updated: March 22, 2018

setwd("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\")
library(tidyverse)
library(vegan)
library(ggthemes)

source("C:\\Users\\wilco\\Desktop\\Git projects\\C2E\\Community Paper\\calculate glass delta.R")

glass_means <- change_glass_d %>%
  group_by(site_project_comm, treatment.x, plot_mani) %>%
  summarise_at(vars(abs_richness_glass:losses_glass), funs(mean)) %>%
  replace_at()

replace_
?replace_all
?prcomp
pca_glass <- prcomp(dplyr::select(glass_means, abs_richness_glass:losses_glass))
