### Calculating cumulative sum of change metrics
### Author: Keivn Wilcox (wilcoxkr@gmail.com)
###
### Last updated March 21, 2018

### Set up workspace
library(tidyverse)
library(ggthemes)
library(RColorBrewer)

setwd("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\")

### Read in data 
change_metrics_bayes <- read.csv("CORRE_RACS_Subset_Perm.csv")

### Calculate cumulative sums of each metric
change_cumsum <- change_metrics_bayes %>%
  group_by(site_project_comm, treatment, plot_id) %>%
  mutate(richness_change_abs = abs(richness_change)) %>%
  mutate(evenness_change_abs = abs(evenness_change)) %>%
  mutate_at(vars(richness_change, richness_change_abs, evenness_change,evenness_change_abs, rank_change, gains, losses), funs(cumsum) ) %>%
  mutate(control = ifelse(plot_mani==0,"control","treatment"))

unique(subset(change_cumsum, site_project_comm=="Alberta_CCD_0")$treatment)

ggplot(subset(change_cumsum,site_project_comm=="CDR_e002_B"), aes(x=treatment_year2, y=rank_change))+
  geom_smooth()
### plot - faceted by site_project_comm and color by treatment
## absolute value of richness change
absrich_plot <- ggplot(subset(change_cumsum, site_project_comm=="Alberta_CCD_0"), aes(x=treatment_year2, y=richness_change_abs, group=treatment)) +
  geom_point(aes(col=control),pch=1) +
  geom_smooth(aes(col=control),se=F) +
#  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
##  richness change
rich_plot <- ggplot(change_cumsum, aes(x=treatment_year2, y=richness_change, group=treatment)) +
  geom_point(aes(col=control),pch=1) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
## evenness change
even_plot <- ggplot(change_cumsum, aes(x=treatment_year2, y=evenness_change, group=treatment)) +
  geom_point(aes(col=control),pch=1) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
## rank change
rank_plot <- ggplot(change_cumsum, aes(x=treatment_year2, y=rank_change, group=treatment)) +
  geom_point(aes(col=control),pch=1) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
## gains change
gains_plot <- ggplot(change_cumsum, aes(x=treatment_year2, y=gains, group=treatment)) +
  geom_point(aes(col=control),pch=1) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
## losses change
losses_plot <- ggplot(change_cumsum, aes(x=treatment_year2, y=losses, group=treatment)) +
  geom_point(aes(col=control),pch=1) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 

png(paste0("figures\\absrich cumsum plot_", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(absrich_plot)
dev.off()

png(paste0("figures\\rich cumsum plot_", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(rich_plot)
dev.off()

png(paste0("figures\\evenness cumsum plot_", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(even_plot)
dev.off()

png(paste0("figures\\rank cumsum plot_", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(rank_plot)
dev.off()

png(paste0("figures\\gains cumsum plot_", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(gains_plot)
dev.off()

png(paste0("figures\\losses cumsum plot_", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(losses_plot)
dev.off()


### plotting just trendlines
## absolute value of richness change
absrich_plot2 <- ggplot(change_cumsum, aes(x=treatment_year2, y=richness_change_abs, group=treatment)) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
##  richness change
rich_plot2 <- ggplot(change_cumsum, aes(x=treatment_year2, y=richness_change, group=treatment)) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
## evenness change
even_plot2 <- ggplot(change_cumsum, aes(x=treatment_year2, y=evenness_change_abs, group=treatment)) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
## rank change
rank_plot2 <- ggplot(change_cumsum, aes(x=treatment_year2, y=rank_change, group=treatment)) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
## gains change
gains_plot2 <- ggplot(change_cumsum, aes(x=treatment_year2, y=gains, group=treatment)) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 
## losses change
losses_plot2 <- ggplot(change_cumsum, aes(x=treatment_year2, y=losses, group=treatment)) +
  geom_smooth(aes(col=control)) +
  facet_wrap(~site_project_comm, scales="free") +
  theme_few() +
  theme(legend.position="none") 

png(paste0("figures\\absrich cumsum plot_trendlines only_perm", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(absrich_plot2)
dev.off()

png(paste0("figures\\rich cumsum plot_trendlines only_perm", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(rich_plot2)
dev.off()

png(paste0("figures\\evenness cumsum plot_trendlines only_perm", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(even_plot2)
dev.off()

png(paste0("figures\\rank cumsum plot_trendlines only_perm", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(rank_plot2)
dev.off()

png(paste0("figures\\gains cumsum plot_trendlines only_perm", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(gains_plot2)
dev.off()

png(paste0("figures\\losses cumsum plot_trendlines only_perm", Sys.Date(),".png"), width=11, height=8, units="in", res=600)
print(losses_plot2)
dev.off()

