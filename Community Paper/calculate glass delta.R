### Calculating Glass's delta
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last updated: March 22, 2018

### Set up workspace
setwd("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\")
library(tidyverse)
library(vegan)
library(ggthemes)

### Read in data 
change_metrics_perm <- read.csv("CORRE_RACS_Subset_Perm.csv") %>%
  mutate(abs_richness_change = abs(richness_change),
         abs_evenness_change = abs(evenness_change))

### Control data
change_control <- change_metrics_perm %>%
  filter(plot_mani==0) %>%
  dplyr::select(treatment_year, treatment_year2, abs_richness_change, abs_evenness_change, 
                rank_change, gains, losses, site_project_comm, treatment, plot_mani) %>%
  rename(abs_richness_change_ctrl = abs_richness_change,
         abs_evenness_change_ctrl = abs_evenness_change,
         rank_change_ctrl = rank_change,
         gains_ctrl = gains,
         losses_ctrl = losses
  ) %>%
  group_by(site_project_comm, treatment, treatment_year2) %>%
  summarise_at(vars(abs_richness_change_ctrl:losses_ctrl), funs(mean, sd))
               
change_glass_d <- change_metrics_perm %>%
  filter(plot_mani != 0) %>%
  group_by(site_project_comm, treatment, treatment_year, treatment_year2, plot_mani) %>%
  summarise(abs_richness_change = mean(abs_richness_change,na.rm=T),
            abs_evenness_change = mean(abs_evenness_change, na.rm=T),
            rank_change = mean(rank_change, na.rm=T),
            gains = mean(gains, na.rm=T),
            losses = mean(losses, na.rm=T)) %>%
  ungroup()%>%
  left_join(change_control, by=c("site_project_comm","treatment_year2")) %>%
  mutate(abs_richness_glass = (abs_richness_change-abs_richness_change_ctrl_mean)/abs_richness_change_ctrl_sd,
         abs_evenness_glass = (abs_evenness_change-abs_evenness_change_ctrl_mean)/abs_evenness_change_ctrl_sd,
         rank_glass = (rank_change-rank_change_ctrl_mean)/rank_change_ctrl_sd,
         gains_glass = (gains-gains_ctrl_mean)/gains_ctrl_sd,
         losses_glass = (losses-losses_ctrl_mean)/losses_ctrl_sd
         ) %>%
  dplyr::select(site_project_comm:plot_mani, abs_richness_glass:losses_glass)




