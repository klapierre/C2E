### SEM to test beta diversity effects on synchrony metrics
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Script created: July 17, 2019; last updated: July 17, 2019

# Updated piecewiseSEM package
install.packages("devtools") # if devtools not installed
devtools::install_github("jslefche/piecewiseSEM@devel")

### Set up workspace
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\GCD asynchrony\\data\\")
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\C2E\\GCD asynchrony\\data') #kim's laptop
library(tidyverse)
library(piecewiseSEM)
library(nlme)
library(lme4)
library(lavaan)
library(PerformanceAnalytics)

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

#summary stats
summary(as.factor(synch_RR_diff_wide$site_code)) #33 sites
summary(as.factor(synch_RR_diff_wide$project_name)) #39 projects
summary(as.factor(synch_RR_diff_wide$community_type)) #51 project*community
length(synch_RR_diff_wide$treatment) #176 treatments


#correlation matrix
chart.Correlation(synch_RR_diff_wide[,c(11, 12, 13, 15, 17, 19, 20)])


# ### Run piecewise SEM
# ## original structure
# dispersion_model_all <- psem(
#   lmer(gamma_stab ~ spatial_synch + alpha_stab + (1|site_code), data=synch_RR_diff_wide),
#   lmer(spatial_synch ~ dispersion_diff + pop_synch + (1|site_code), data=synch_RR_diff_wide),
#   lmer(alpha_stab ~ spp_stab + spp_synch + (1|site_code), data=synch_RR_diff_wide),
#   lmer(pop_synch ~ dispersion_diff + (1|site_code), data=synch_RR_diff_wide),
#   spp_synch %~~% spatial_synch,
#   alpha_stab %~~% spatial_synch,
#   spp_synch %~~% dispersion_diff,
#   spp_stab %~~% spatial_synch,
#   spp_stab %~~% spp_synch,
#   spp_stab %~~% pop_synch,
#   spp_synch %~~% pop_synch
# )
# 
# summary(dispersion_model_all)
# 
# coefs(dispersion_model_all, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
#   select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)


## spatial only
test <- na.omit(synch_RR_diff_wide)

dispersion_model_all <- psem(
  lme(gamma_stab ~ spatial_synch, random = ~ 1 | site_code, na.action = na.omit, data = synch_RR_diff_wide),
  lme(spatial_synch ~ dispersion_diff + pop_synch, random = ~ 1 | site_code, na.action = na.omit, data = synch_RR_diff_wide),
  lme(pop_synch ~ dispersion_diff, random = ~ 1 | site_code, na.action = na.omit, data = synch_RR_diff_wide),
  gamma_stab %~~% dispersion_diff
)

summary(dispersion_model_all)

coefs(dispersion_model_all, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)


## no random effects, spatial only
dispersion_model_all <- psem(
  lm(gamma_stab ~ spatial_synch, data = synch_RR_diff_wide),
  lm(spatial_synch ~ dispersion_diff + pop_synch, data = synch_RR_diff_wide),
  lm(pop_synch ~ dispersion_diff, data = synch_RR_diff_wide),
  gamma_stab %~~% dispersion_diff
)

summary(dispersion_model_all)

coefs(dispersion_model_all, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)




# 
# ### Trying to figure out what happens when variables occur in different orders
# ### flipped pop and spatial synch structure
# dispersion_model_flipped <- psem(
#   lm(gamma_stab ~ pop_synch + alpha_stab, data=synch_RR_diff_wide),
#   lm(pop_synch ~ dispersion_diff + spatial_synch, data=synch_RR_diff_wide),
#   lm(alpha_stab ~ spp_stab + spp_synch, data=synch_RR_diff_wide),
#   lm(spatial_synch ~ dispersion_diff, data=synch_RR_diff_wide),
#   spp_synch %~~% dispersion_diff,
#   spp_synch %~~% spatial_synch,
#   spp_stab %~~% spatial_synch,
#   spp_stab %~~% spp_synch,
#   spp_stab %~~% pop_synch,
#   alpha_stab %~~% spatial_synch,
#   spp_synch %~~% pop_synch
# )
# 
# summary(dispersion_model_flipped)
# 
# coefs(dispersion_model_flipped, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
#   select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)
# 
# ## flipped pop and spatial synch structure
# dispersion_model_nopop <- psem(
#   lm(gamma_stab ~ spatial_synch + alpha_stab, data=synch_RR_diff_wide),
#   lm(alpha_stab ~ spp_stab + spp_synch, data=synch_RR_diff_wide),
#   lm(spatial_synch ~ dispersion_diff, data=synch_RR_diff_wide),
#   spp_synch %~~% dispersion_diff,
#   spp_synch %~~% spatial_synch,
#   spp_stab %~~% spatial_synch,
#   spp_stab %~~% spp_synch,
#   alpha_stab %~~% spatial_synch
# )
# 
# summary(dispersion_model_nopop)
# 
# coefs(dispersion_model_all, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
#   select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)
# 
# 
# dispersion_model_all3 <- psem(
#   lm(gamma_stab ~ spatial_synch + pop_synch, data=synch_RR_diff_wide),
#   lm(spatial_synch ~ dispersion_diff + pop_synch, data=synch_RR_diff_wide),
#   lm(pop_synch ~ dispersion_diff, data=synch_RR_diff_wide)
# )
# 
# 
# sem.fit(dispersion_model_all3)
# summary(dispersion_model_all3)
# str(dispersion_model_all3)
# with(synch_RR_diff_wide, plot(spatial_asynch,gamma_stab))
# 
# with(synch_RR_diff_wide, cor(pop_asynch, gamma_stab))
# 
# 
# summary(dispersion_model_all3)
# 
# ### Chekcing residuals
# gamma_resids <- residuals(lm(gamma_stab ~ pop_asynch, data=synch_RR_diff_wide))
# plot(gamma_resids, synch_RR_diff_wide$spatial_asynch)
# cor(gamma_resids, synch_RR_diff_wide$spatial_asynch)
# 
# #
# #
# # dispersion_model_all <- psem( ### flipping arrows
# #   lm(gamma_stab ~ spatial_asynch, data=synch_RR_diff_wide),
# #   lm(dispersion_diff ~ spatial_asynch + pop_asynch, data=synch_RR_diff_wide),
# #   lm(spatial_asynch ~ pop_asynch, data=synch_RR_diff_wide)
# # )
# #
# #


dispersion_model_spatial_lavaan <- "
  gamma_stab ~ spatial_synch
  spatial_synch ~ dispersion_diff + pop_synch
  pop_synch ~ dispersion_diff

  gamma_stab ~~ dispersion_diff
 "


dispersion_model_all_output <- sem(dispersion_model_spatial_lavaan, data=synch_RR_diff_wide)

summary(dispersion_model_all_output)






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

