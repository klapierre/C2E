### Plotting PCA
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last updated: March 22, 2018

setwd("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\")
setwd("/Users/meghanavolio/Dropbox/C2E/Products/CommunityChange/March2018 WG/")
library(tidyverse)
library(vegan)
library(ggthemes)

source("C:\\Users\\wilco\\Desktop\\Git projects\\C2E\\Community Paper\\calculate glass delta.R")
source("/Users/meghanavolio/Documents/R Projects/C2E/Community Paper/calculate glass delta.R")


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


pca_glass <- prcomp(dplyr::select(glass_means, abs_richness_glass:losses_glass), scale.=T)

### combine treatment and site information
trt_info <- read.csv("ExperimentInformation_Nov2017.csv") %>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_"))
site_info <- read.csv("SiteExperimentDetails_Dec2016.csv") %>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_"))

###
pca_df <- data.frame(dplyr::select(glass_means, site_project_comm:plot_mani),
                     pca_glass$x) %>%
  rename(treatment=treatment.x) %>%
  left_join(trt_info, by=c("site_project_comm","treatment")) %>%
  left_join(site_info, by=c("site_project_comm"))

pca_df_long <- pca_df %>%
  gather(key=site_var, value=site_value, -(site_project_comm:project_name.y), -public.y)

ggplot(pca_df_long, aes(x=PC1, y=PC2, col=site_value)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~site_var)+
  theme(legend.position="none")

a <- ggplot(pca_df, aes(x=PC1, y=PC2, col=MAP)) +
  geom_point() +
  theme_bw() +
  xlim(-5,5) +
  ylim(-5,2.5)
#  theme(legend.position="none")

b <-ggplot(pca_df, aes(x=PC1, y=PC2, col=MAT)) +
  geom_point() +
  theme_bw() +
  xlim(-5,5) +
  ylim(-5,2.5)

c <- ggplot(pca_df, aes(x=PC1, y=PC2, col=experiment_length)) +
  geom_point() +
  theme_bw() +
  xlim(-5,5) +
  ylim(-5,2.5)

d <- ggplot(pca_df, aes(x=PC1, y=PC2, col=rrich)) +
  geom_point() +
  theme_bw() +
  xlim(-5,5) +
  ylim(-5,2.5)
install.packages("gridExtra")
library(gridExtra)
grid.arrange(a,b,c,d)
