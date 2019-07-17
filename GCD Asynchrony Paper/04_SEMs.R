### SEM to test beta diversity effects on synchrony metrics
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Script created: July 17, 2019; last updated: July 17, 2019

### Set up workspace
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\C2E_May2019\\asynchrony\\data\\")
library(tidyverse)
library(piecewiseSEM)

### Read in synchrony metrics
synch_metrics_sub <- read.csv("..\\synchrony metrics with environmental_subset.csv")

### Read in difference metrics and combine
rac_diff_metrics <- read.csv("..\\corre_community differences_March2019.csv") %>%
  dplyr::select(-treatment) %>%
  rename(treatment=treatment2, treatment_year=time) %>%
  mutate(dispersion_diff = ifelse(trt_greater_disp=="Control", abs_dispersion_diff, -(abs_dispersion_diff) )) %>%
  group_by(site_code, project_name, community_type, treatment) %>%
  summarise_at(.vars=vars(composition_diff, richness_diff:species_diff, dispersion_diff), mean, na.rm=T)

synch_RR_diff_wide <- synch_metrics_sub %>%
  filter(site_code != "NANT") %>% ## many plots have 1 species so synch metrics cannot be calculated
  dplyr::select(site_code:community_type, metric_name, lnRR, trt_type2) %>%
  left_join(rac_diff_metrics, by=c("site_code","project_name", "community_type", "treatment")) %>%
  filter(trt_type2 %in% c("CO2", "drought", "Irrigation", "Precip. Vari.", "Temperature", "N", "P", "Mult. Nuts.")) %>%
  spread(key=metric_name, value=lnRR)

### Run piecewise SEM
dispersion_model_all <- psem(
  lm(gamma_stab ~ spatial_asynch, data=synch_RR_diff_wide),
  lm(spatial_asynch ~ dispersion_diff + pop_asynch, data=synch_RR_diff_wide),
  lm(pop_asynch ~ dispersion_diff, data=synch_RR_diff_wide)
)

summary(dispersion_model_all)

coefs(dispersion_model_all, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)

# 
# 
# dispersion_model_all <- psem( ### flipping arrows
#   lm(gamma_stab ~ spatial_asynch, data=synch_RR_diff_wide),
#   lm(dispersion_diff ~ spatial_asynch + pop_asynch, data=synch_RR_diff_wide),
#   lm(spatial_asynch ~ pop_asynch, data=synch_RR_diff_wide)
# )
# 
# 
# dispersion_model_all_lavaan <- "
#   level:1
#   gamma_stab ~ spatial_asynch
#   spatial_asynch ~ dispersion_diff + pop_asynch
#   pop_asynch ~ dispersion_diff
#   level:2
#   gamma_stab ~~ spatial_asynch + pop_asynch + dispersion_diff
#   spatial_asynch ~~ pop_asynch + dispersion_diff
#   pop_asynch ~~ dispersion_diff
# "
# dispersion_model_all_output <- sem(dispersion_model_all_lavaan, data=synch_RR_diff_wide,
#                                    cluster="site_code", optim.method="em")
# 
# with(synch_RR_diff_wide, plot(gamma_stab, spatial_asynch))
# with(synch_RR_diff_wide, plot(pop_asynch, spatial_asynch))
# with(synch_RR_diff_wide, plot(dispersion_diff, spatial_asynch))
# with(synch_RR_diff_wide, plot(dispersion_diff, pop_asynch))
# 
# summary(dispersion_model_all)
# summary(compositionModel10 <- psem(
#   lm(anpp_pdiff ~ n + p + k + richness_difference + evenness_diff + rank_difference + species_difference, data=subset(allSEMdata, treatment_year==10)),
#   lm(richness_difference ~ n + p + k, data=subset(allSEMdata, treatment_year==10)),
#   lm(evenness_diff ~ n + p + k, data=subset(allSEMdata, treatment_year==10)),
#   lm(rank_difference ~ n + p + k, data=subset(allSEMdata, treatment_year==10)),
#   lm(species_difference ~ n + p + k, data=subset(allSEMdata, treatment_year==10))
# ))
# coefs10 <- coefs(compositionModel10, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
#   select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
#   mutate(treatment_year=10)
# 
# 
# spatial_model <- lm(spatial_asynch ~ dispersion_diff + pop_asynch, data=synch_RR_diff_wide)
# anova(spatial_model)

