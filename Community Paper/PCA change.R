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
  mutate(abs_richness_glass = replace(abs_richness_glass, abs_richness_glass == "Inf", NA),
         abs_evenness_glass = replace(abs_evenness_glass, abs_evenness_glass == "Inf", NA),
         rank_glass = replace(rank_glass, rank_glass == "Inf", NA),
         gains_glass = replace(gains_glass, gains_glass == "Inf", NA),
         losses_glass = replace(losses_glass, losses_glass == "Inf", NA)) %>%
  ungroup() %>%
  na.omit()

filter(glass_means, is.na(abs_evenness_glass))

pca_glass <- prcomp(dplyr::select(glass_means, abs_richness_glass:losses_glass), scale.=T)

pca_glass$x
