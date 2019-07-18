### Compare lnRR with Site-level variables
###
### Author: Kevin wilcox (kevin.wilcox@uwyo.edu)
### Created: May 14th 2019, last updated: May 15th, 2019

### Set up workspace
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\C2E_May2019\\asynchrony\\data\\")
library(tidyverse)
library(car)
library(MASS)
library(lsmeans)

### Read in synchrony metrics
synch_metrics_sub <- read.csv("..\\synchrony metrics with environmental_subset.csv")

### Read in difference metrics
site_info <- read.csv("SiteExperimentDetails_March2019.csv") %>%
  dplyr::select(-treatment) %>%
  rename(treatment=treatment2, treatment_year=time) %>%
  mutate(dispersion_diff = ifelse(trt_greater_disp=="Control", abs_dispersion_diff, -(abs_dispersion_diff) )) %>%
  group_by(site_code, project_name, community_type, treatment) %>%
  summarise_at(.vars=vars(composition_diff, richness_diff:species_diff, dispersion_diff), mean, na.rm=T)

synch_RR_site_info <- synch_metrics_sub %>%
  dplyr::select(site_code:community_type, metric_name, lnRR, trt_type2) %>%
  left_join(site_info, by=c("site_code","project_name", "community_type")) %>%
  mutate(
    MAP_norm = as.vector(scale(MAP)),
    MAT_norm =  as.vector(scale(MAT)),
    rrich_norm = as.vector(scale(rrich)),
    experiment_length_norm =  as.vector(scale(experiment_length)),
    anpp_norm =  as.vector(scale(anpp))
  )
  


### Run model loop to generate slopes and se's for each life form and functional group
trt_type2_vec <- unique(synch_RR_site_info$trt_type2) ## vector of metrics to cycle through

model_out <- {} ## set up data frame to fill with model output

for(TRT in 1:length(trt_type2_vec)){ # loop through traits
  
  df_temp <- synch_RR_site_info %>% # create temporary data frame for one trait at a time
    filter(trt_type2==trt_type2_vec[TRT])
  
  # run model with all clusters together
  #  ppt_weights_temp <- 1/df_trait_temp$ppt_slope_se
  
  MAP_model_temp <- lm(lnRR ~ MAP_norm*metric_name, data = df_temp)
  MAT_model_temp <- lm(lnRR ~ MAT_norm*metric_name, data = df_temp)
  anpp_model_temp <- lm(lnRR ~ anpp_norm*metric_name, data = df_temp)
  rrich_model_temp <- lm(lnRR ~ rrich_norm*metric_name, data = df_temp)
  experiment_length_model_temp <- lm(lnRR ~ experiment_length_norm*metric_name, data = df_temp)

  MAP_lstrends <- summary(lstrends(MAP_model_temp,  "metric_name", var="MAP_norm")) %>%
    mutate(site_char = "MAP") %>%
    mutate(trt_type2 = as.character(trt_type2_vec[TRT])) %>%
    rename(slope = MAP_norm.trend)
  MAT_lstrends <- summary(lstrends(MAT_model_temp,  "metric_name", var="MAT_norm")) %>%
    mutate(site_char = "MAT") %>%
    mutate(trt_type2 = as.character(trt_type2_vec[TRT])) %>%
    rename(slope = MAT_norm.trend)
  anpp_lstrends <- summary(lstrends(anpp_model_temp,  "metric_name", var="anpp_norm")) %>%
    mutate(site_char = "anpp") %>%
    mutate(trt_type2 = as.character(trt_type2_vec[TRT])) %>%
    rename(slope = anpp_norm.trend)
  rrich_lstrends <- summary(lstrends(rrich_model_temp,  "metric_name", var="rrich_norm")) %>%
    mutate(site_char = "rrich") %>%
    mutate(trt_type2 = as.character(trt_type2_vec[TRT])) %>%
    rename(slope = rrich_norm.trend)
  experiment_length_lstrends <- summary(lstrends(experiment_length_model_temp,  "metric_name", var="experiment_length_norm")) %>%
    mutate(site_char = "experiment_length") %>%
    mutate(trt_type2 = as.character(trt_type2_vec[TRT])) %>%
    rename(slope = experiment_length_norm.trend)
  
  model_out_temp <- rbind(MAP_lstrends, MAT_lstrends, anpp_lstrends, rrich_lstrends, experiment_length_lstrends)
  model_out <- rbind(model_out, model_out_temp)
}

### Overall trends
MAP_model_temp <- lm(lnRR ~ MAP_norm*metric_name, data = synch_RR_site_info)
MAT_model_temp <- lm(lnRR ~ MAT_norm*metric_name, data = synch_RR_site_info)
anpp_model_temp <- lm(lnRR ~ anpp_norm*metric_name, data = synch_RR_site_info)
rrich_model_temp <- lm(lnRR ~ rrich_norm*metric_name, data = synch_RR_site_info)
experiment_length_model_temp <- lm(lnRR ~ experiment_length_norm*metric_name, data = synch_RR_site_info)

MAP_lstrends <- summary(lstrends(MAP_model_temp,  "metric_name", var="MAP_norm")) %>%
  mutate(site_char = "MAP") %>%
  mutate(trt_type2 = "All") %>%
  rename(slope = MAP_norm.trend)
MAT_lstrends <- summary(lstrends(MAT_model_temp,  "metric_name", var="MAT_norm")) %>%
  mutate(site_char = "MAT") %>%
  mutate(trt_type2 = "All") %>%
  rename(slope = MAT_norm.trend)
anpp_lstrends <- summary(lstrends(anpp_model_temp,  "metric_name", var="anpp_norm")) %>%
  mutate(site_char = "anpp") %>%
  mutate(trt_type2 = "All") %>%
  rename(slope = anpp_norm.trend)
rrich_lstrends <- summary(lstrends(rrich_model_temp,  "metric_name", var="rrich_norm")) %>%
  mutate(site_char = "rrich") %>%
  mutate(trt_type2 = "All") %>%
  rename(slope = rrich_norm.trend)
experiment_length_lstrends <- summary(lstrends(experiment_length_model_temp,  "metric_name", var="experiment_length_norm")) %>%
  mutate(site_char = "experiment_length") %>%
  mutate(trt_type2 = "All") %>%
  rename(slope = experiment_length_norm.trend)

model_out_temp <- rbind(MAP_lstrends, MAT_lstrends, anpp_lstrends, rrich_lstrends, experiment_length_lstrends)
model_out <- rbind(model_out, model_out_temp)
model_out_subset <- model_out %>%
  filter(trt_type2 %in% c(
    "drought",
    "Irrigation",
    "Precip. Vari.",
    "Temperature",
    "Mult. Nuts.",
    "N",
    "P",
    "All"
  ))
###
### Visualize slopes
###

ggplot(model_out_subset, aes(x=site_char, y=slope, ymin=lower.CL, ymax=upper.CL, col=trt_type2)) +
  geom_hline(yintercept=0) +
  geom_point(position=position_dodge(width=0.9)) +
  geom_errorbar(position=position_dodge(width=0.9), width=0) +
  facet_wrap(~metric_name, scales="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=60, hjust=1))





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

