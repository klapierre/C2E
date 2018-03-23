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
change_metrics_perm <- read.csv("CORRE_RACS_Subset_Perm.csv") %>%
  mutate(abs_richness_change = abs(richness_change),
         abs_evenness_change = abs(evenness_change))
