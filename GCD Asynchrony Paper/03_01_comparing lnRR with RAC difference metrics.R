### Compare lnRR with difference metrics
###
### Author: Kevin wilcox (kevin.wilcox@uwyo.edu)
### Created: May 14th 2019, last updated: May 21, 2020

### Set up workspace
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\GCD asynchrony\\data\\") # Kevin's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\C2E\\GCD asynchrony\\data') #kim's laptop
setwd("C:\\Users\\kwilcox4\\Dropbox\\Shared working groups\\C2E\\GCD asynchrony\\data\\") # Kevin's work comp

library(codyn)
library(tidyverse)
library(ggthemes)

### Read in synchrony metrics
synch_metrics_sub <- read.csv("..\\synchrony metrics with environmental_subset_22Apr2020.csv") %>%
  filter(!metric_name %in% c("spp_asynch","spatial_asynch","pop_asynch"))

### Read in difference metrics
rac_diff_metrics <- read.csv("..\\corre_community differences_March2019.csv") %>%
  dplyr::select(-treatment) %>%
  rename(treatment=treatment2, treatment_year=time) %>%
  mutate(dispersion_diff = ifelse(trt_greater_disp=="Control", abs_dispersion_diff, -(abs_dispersion_diff) )) %>%
  group_by(site_code, project_name, community_type, treatment) %>%
  filter(treatment_year == max(treatment_year)) %>%
  filter(treatment_year > 4)

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
  facet_wrap(~metric_name, scales="free")

comp_plot_all <- ggplot(synch_rac_for_plotting, aes(x=composition_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_few() +
  facet_wrap(~metric_name, scales="free")

disp_plot_all <- ggplot(synch_rac_for_plotting, aes(x=dispersion_diff, y=lnRR, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  geom_hline(yintercept=0) +
  theme_few() +
  facet_wrap(~metric_name, scales="free")

disp_plot_all <- ggplot(synch_rac_for_plotting, aes(x=dispersion_diff, y=lnRR)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0) +
  theme_few() +
  facet_wrap(~metric_name, scales="free")

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


## Beta diversity box plots
ggplot(filter(synch_rac_for_plotting, metric_name=="alpha_stab"), aes(x=trt_type2, y=dispersion_diff, col=trt_type2)) +
  geom_hline(yintercept=0) +
  geom_boxplot() +
  geom_jitter() +
  theme_few()


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

synch_metric_sub_wide <- synch_metrics_sub %>%
  dplyr::select(site_code:community_type, site_proj_comm, trt_type2, metric_name, lnRR) %>%
  spread(key=metric_name, value=lnRR)

### gamma_stability ~ spatial_synchrony
ggplot(synch_metric_sub_wide, aes(x=spatial_synch, y=gamma_stab, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme_few()

ggplot(synch_metric_sub_wide, aes(x=spatial_synch, y=gamma_stab)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few()

### gamma_stability ~ alpha_stability
ggplot(synch_metric_sub_wide, aes(x=alpha_stab, y=gamma_stab, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme_few()

ggplot(synch_metric_sub_wide, aes(x=alpha_stab, y=gamma_stab)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few()

### spatial_synchrony ~ population_synchrony
ggplot(synch_metric_sub_wide, aes(y=spatial_synch, x=pop_synch, col=trt_type2)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme_few()

ggplot(synch_metric_sub_wide, aes(y=spatial_synch, x=pop_synch)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few()

### spatial_synchrony ~ all difference metrics
ggplot(filter(synch_RR_diff_metrics,metric_name=="spatial_synch"), aes(x=composition_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="spatial_synch"), aes(x=abs_dispersion_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="spatial_synch"), aes(x=richness_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="spatial_synch"), aes(x=evenness_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="spatial_synch"), aes(x=rank_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="spatial_synch"), aes(x=species_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="spatial_synch"), aes(x=dispersion_diff, y=lnRR, col=trt_type2, label=project_name)) +
  geom_point() + theme_few() + geom_hline(yintercept=0) + geom_text(size=2, nudge_x=0.01)

### pop_synchrony ~ all difference metrics
ggplot(filter(synch_RR_diff_metrics,metric_name=="pop_synch"), aes(x=composition_diff, y=lnRR, col=trt_type2, label=project_name)) +
  geom_point() + theme_few() + geom_hline(yintercept=0) + geom_text(size=2, nudge_x=0.03)
ggplot(filter(synch_RR_diff_metrics,metric_name=="pop_synch"), aes(x=composition_diff, y=lnRR)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)+geom_smooth(method="lm")
ggplot(filter(synch_RR_diff_metrics,metric_name=="pop_synch"), aes(x=abs_dispersion_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="pop_synch"), aes(x=richness_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="pop_synch"), aes(x=evenness_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="pop_synch"), aes(x=rank_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="pop_synch"), aes(x=species_diff, y=lnRR, col=trt_type2)) +
  geom_point() + theme_few() + geom_hline(yintercept=0)
ggplot(filter(synch_RR_diff_metrics,metric_name=="pop_synch"), aes(x=dispersion_diff, y=lnRR, col=trt_type2, label=project_name)) +
  geom_point() + theme_few() + geom_hline(yintercept=0) + geom_text(size=2, nudge_x=0.01)

