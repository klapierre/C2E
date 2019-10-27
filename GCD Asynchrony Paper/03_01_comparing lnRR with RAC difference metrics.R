### Compare lnRR with difference metrics
###
### Author: Kevin wilcox (kevin.wilcox@uwyo.edu)
### Created: May 14th 2019, last updated: June 13th, 2019

### Set up workspace
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\C2E_May2019\\asynchrony\\data\\")
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\C2E\\GCD asynchrony\\data') #kim's laptop
library(codyn)
library(tidyverse)

### Read in synchrony metrics
synch_metrics_sub <- read.csv("..\\synchrony metrics with environmental_subset.csv") %>%
  filter(!metric_name %in% c("spp_asynch","spatial_asynch","pop_asynch"))

### Read in difference metrics
rac_diff_metrics <- read.csv("..\\corre_community differences_March2019.csv") %>%
  dplyr::select(-treatment) %>%
  rename(treatment=treatment2, treatment_year=time) %>%
  mutate(dispersion_diff = ifelse(trt_greater_disp=="Control", abs_dispersion_diff, -(abs_dispersion_diff) )) %>%
  group_by(site_code, project_name, community_type, treatment) %>%
  summarise_at(.vars=vars(composition_diff, richness_diff:species_diff, dispersion_diff), mean, na.rm=T)

# dispersion_df <- read.csv("..\\ForBayesianAnalysis_May2017.csv") %>%
#   group_by(site_code, project_name, community_type, treatment) %>%
#   summarize(dispersion_change=mean(dispersion_change, na.rm=T)) KIM'S OLD DISPERSION VALUES

synch_RR_diff_metrics <- synch_metrics_sub %>%
  dplyr::select(site_code:community_type, metric_name, lnRR, trt_type2) %>%
  left_join(rac_diff_metrics, by=c("site_code","project_name", "community_type", "treatment"))

# synch_RR_diff_metrics <- synch_metrics_sub %>%
#   dplyr::select(site_code:community_type, metric_name, lnRR, trt_type2) %>%
#   left_join(dispersion_df, by=c("site_code","project_name", "community_type", "treatment"))


###
### Visualize
###
synch_rac_for_plotting <- synch_RR_diff_metrics %>%
  filter(trt_type2 %in% c("CO2", "drought", "Irrigation", "Precip. Vari.", "Temperature", "N", "P", "Mult. Nuts."))

synch_rac_for_plotting$trt_type2 <- factor(synch_rac_for_plotting$trt_type2,
                                           levels=c("CO2", "drought", "Irrigation", "Precip. Vari.", "Temperature", "N", "P", "Mult. Nuts."))
  
  ## All treatments
comp_plot_all <- ggplot(synch_rac_for_plotting, aes(x=composition_diff, y=lnRR, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  geom_hline(yintercept=0) +
  theme_few() +
  facet_wrap(~metric_name)

disp_plot_all <- ggplot(synch_rac_for_plotting, aes(x=dispersion_diff, y=lnRR, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  geom_hline(yintercept=0) +
  theme_few() +
  facet_wrap(~metric_name)

spp_plot_all <- ggplot(synch_rac_for_plotting, aes(x=species_diff, y=lnRR, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  geom_hline(yintercept=0) +
  # theme_few() +
  facet_wrap(~metric_name)

pdf("..//figures//composition_fig_17July2019.pdf", width=9, height=5, useDingbats=F)
print(comp_plot_all)
dev.off()

pdf("..//figures//dispersion_fig_17July2019.pdf", width=9, height=5, useDingbats=F)
print(disp_plot_all)
dev.off()

pdf("..//figures//sppdiff_fig_27Oct2019.pdf", width=9, height=5, useDingbats=F)
print(spp_plot_all)
dev.off()

## N
ggplot(subset(synch_RR_diff_metrics, trt_type2=="N"), aes(x=composition_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)
 
ggplot(subset(synch_RR_diff_metrics, trt_type2=="N"), aes(x=dispersion_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)
 
## Mult. Nuts.
ggplot(subset(synch_RR_diff_metrics, trt_type2=="Mult. Nuts."), aes(x=composition_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)

ggplot(subset(synch_RR_diff_metrics, trt_type2=="Mult. Nuts."), aes(x=dispersion_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)

## irr
ggplot(subset(synch_RR_diff_metrics, trt_type2=="Irrigation"), aes(x=composition_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)

ggplot(subset(synch_RR_diff_metrics, trt_type2=="Irrigation"), aes(x=dispersion_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)

## drought
ggplot(subset(synch_RR_diff_metrics, trt_type2=="drought"), aes(x=composition_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)

ggplot(subset(synch_RR_diff_metrics, trt_type2=="drought"), aes(x=dispersion_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)

## temp
ggplot(subset(synch_RR_diff_metrics, trt_type2=="Temperature"), aes(x=composition_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)

ggplot(subset(synch_RR_diff_metrics, trt_type2=="Temperature"), aes(x=dispersion_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)


###
### look at metrics, facet by treatement type
###
trt_type_vec <- c("drought","Irrigation","Precip. Vari.", "Mult. Nuts.","N","P","Temperature")

ggplot(subset(synch_RR_diff_metrics, metric_name=="gamma_stab" & trt_type2 %in% trt_type_vec), aes(x=dispersion_diff, y=lnRR, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~trt_type2)

ggplot(subset(synch_RR_diff_metrics, trt_type2 %in% trt_type_vec), aes(x=dispersion_diff, y=lnRR, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~metric_name)

ggplot(subset(synch_RR_diff_metrics, trt_type2 %in% trt_type_vec), aes(x=richness_diff, y=dispersion_diff, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  geom_hline(yintercept=0) +
  theme_classic() +
  facet_wrap(~trt_type2)

