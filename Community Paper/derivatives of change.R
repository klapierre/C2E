### Calculating derivatives of change metrics
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last updated: March 23, 2018

### Set up workspace
setwd("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\")
library(tidyverse)
library(vegan)
library(ggthemes)

### Read in data 
change_metrics_bayes <- read.csv("CORRE_RACS_Subset_Perm.csv")

### Calculate cumulative sums of each metric
change_cumsum <- change_metrics_bayes %>%
  group_by(site_project_comm, treatment, plot_id) %>%
  mutate(richness_change_abs = abs(richness_change)) %>%
  mutate(evenness_change_abs = abs(evenness_change)) %>%
  mutate_at(vars(richness_change, richness_change_abs, evenness_change,evenness_change_abs, rank_change, gains, losses), funs(cumsum) ) %>%
  mutate(control = ifelse(plot_mani==0,"control","treatment"))

### Calculate derivatives for each timestep


df <- data.frame(a=rep(c(rep(1,5),rep(2,5)),50), b=rep(c(1:5,1:5),50))
ggplot(df, aes(x=a, y=b)) +
  geom_jitter(alpha=.5)

