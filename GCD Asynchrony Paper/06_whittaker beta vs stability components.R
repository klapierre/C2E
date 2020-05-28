### Comparing Whittaker beta diversity with response ratios of stability components
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### created Nov 15, 2019 (last updated Nov 15, 2019)



###
### TO DO, CALCULATE BETA WITH WHITTAKER AND GET DIFFERENCE, MERGE WITH STABILITY COMPONENTS RR, REDO ALL PLOTS AND SEM WITH WHITTAKER
###

### Set up workspace
rm(list=ls())
library(codyn)
library(tidyverse)
library(ggthemes)

setwd("C:\\Users\\kwilcox4\\Dropbox\\Shared working groups\\C2E\\GCD asynchrony\\data\\") # desktop
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\GCD asynchrony\\data\\") # laptop

### Read in synchrony metrics
synch_metrics_sub <- read.csv("..\\synchrony metrics with environmental_subset.csv") %>%
  filter(!metric_name %in% c("spp_asynch","spatial_asynch","pop_asynch"))

### Read in community data
cover_df <- read.csv('SpeciesRawAbundance_March2019.csv', header = T)
exp_info <- read.csv('ExperimentInformation_March2019.csv', header = T) %>%
  mutate(site_proj_comm_trt = paste(site_code, project_name, community_type, treatment, sep='_'))
site_info <- read.csv('SiteExperimentDetails_March2019.csv')
exp_info_short <- unique(dplyr::select(exp_info, -X, -calendar_year, -treatment_year))

cover_expInfo_df <- cover_df %>%
  full_join(exp_info, by=c("site_code", 
                           "project_name", 
                           "calendar_year",
                           "treatment_year",
                           "treatment",
                           "community_type"))

rm(cover_df, exp_info, exp_info_short, site_info)

### Create new treatment column and subset treatments and studies to those we use for synchrony metrics
old_metrics_vec <- read.csv("ecolett paper_site_proj_comm vector.csv") %>%
  .$site_proj_comm

cover_select_treatments <- cover_expInfo_df %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  filter(site_proj_comm %in% c(as.character(old_metrics_vec), "CHY_EDGE_0","HYS_EDGE_0","KNZ_EDGE_0","SGS_EDGE_0","SEV_EDGE_EB","SEV_EDGE_EG")) %>%
  filter(pulse == 0) %>% ## remove all pulse treatments
  filter(plant_trt == 0) %>% ## remove seeding treatments
  mutate(trt_type2=ifelse(trt_type=="N","N", 
                          ifelse(trt_type=="P", "P", 
                                 ifelse(trt_type=="CO2", "CO2",
                                        ifelse(trt_type=="irr", "Irrigation",
                                               ifelse(trt_type=="temp", "Temperature", 
                                                      ifelse(trt_type=="N*P"|trt_type=="mult_nutrient", "Mult. Nuts.", 
                                                             ifelse(trt_type=="drought", "drought", 
                                                                    ifelse(trt_type=="CO2*temp", "CO2*temp", 
                                                                           ifelse(trt_type=="drought*temp", "drought*temp", 
                                                                                  ifelse(trt_type=="irr*temp", "irr*temp",
                                                                                         ifelse(trt_type=="irr*CO2*temp"|trt_type=="N*CO2*temp"|trt_type=="N*irr*temp"|trt_type=="N*irr*CO2*temp", "mult_res*temp", 
                                                                                                ifelse(trt_type=="irr*herb_removal"|trt_type=="irr*plant_mani"|trt_type=="irr*plant_mani*herb_removal", "irr*NR", 
                                                                                                       ifelse(trt_type=="herb_removal"|trt_type=="till"|trt_type=="mow_clip"|trt_type=="burn"|trt_type=="plant_mani"|trt_type=="stone"|trt_type=="graze"|trt_type=="burn*graze"|trt_type=="fungicide"|trt_type=="plant_mani*herb_removal"|trt_type=="burn*mow_clip", "NR", 
                                                                                                              ifelse(trt_type=="precip_vari", "Precip. Vari.",  
                                                                                                                     ifelse(trt_type=="N*plant_mani"|trt_type=="N*burn"|trt_type=="N*mow_clip"|trt_type=="N*till"|trt_type=="N*stone"|trt_type=="N*burn*graze"|trt_type=="N*burn*mow_clip", "N*NR", 
                                                                                                                            ifelse(trt_type=="N*temp", "N*temp", 
                                                                                                                                   ifelse(trt_type=="N*CO2", "N*CO2",
                                                                                                                                          ifelse(trt_type=="irr*CO2", "irr*CO2",
                                                                                                                                                 ifelse(trt_type=="N*irr", "N*irr",
                                                                                                                                                        ifelse(trt_type=="mult_nutrient*herb_removal"|trt_type=="mult_nutrient*fungicide"|trt_type=="N*P*burn*graze"|trt_type=="N*P*burn"|trt_type=="*P*mow_clip"|trt_type=="N*P*burn*mow_clip"|trt_type=="N*P*mow_clip", "mult_nutrients*NR",
                                                                                                                                                               ifelse(trt_type=="P*mow_clip"|trt_type=="P*burn"|trt_type=="P*burn*graze"|trt_type=="P*burn*mow_clip", "P*NR", 
                                                                                                                                                                      ifelse(trt_type=="precip_vari*temp", "precip_vari*temp",
                                                                                                                                                                             ifelse(trt_type=="light","light",
                                                                                                                                                                                    ifelse(trt_type=="control","control",
                                                                                                                                                                                    ifelse(trt_type=="N*irr*CO2", "mult_res", 999))))))))))))))))))))))))))

### Vector of treatment types that have at least 5 sites where they are implemented
trt_type2_vector <- c(
  "control",
  "drought",
  "Irrigation",
  "Temperature",
  "Mult. Nuts.",
  "N",
  "P",
  "Precip. Vari."
)


cover_select_treatments <- cover_select_treatments %>%
  filter(trt_type2 %in% trt_type2_vector)

unique(cover_select_treatments$site_code)
alpha_richness_df <- cover_select_treatments %>%
  filter(!site_proj_comm %in% c("SERC_TMECE_MX", "SERC_TMECE_SP", "SERC_TMECE_SC")) %>%
  group_by(site_proj_comm, trt_type2, treatment, calendar_year, plot_id) %>%
  summarize(alpha_richness = specnumber(abundance)) %>%
  ungroup() %>%
  group_by(site_proj_comm, trt_type2, treatment) %>%
  summarize(alpha_richness = mean(alpha_richness)) %>%
  ungroup()

gamma_richness_df <- cover_select_treatments %>%
  filter(!site_proj_comm %in% c("SERC_TMECE_MX", "SERC_TMECE_SP", "SERC_TMECE_SC")) %>%
  group_by(site_proj_comm, trt_type2, treatment, calendar_year, genus_species) %>%
  summarize(abundance = mean(abundance, na.rm=T)) %>%
  ungroup() %>%
  group_by(site_proj_comm, trt_type2, treatment, calendar_year) %>%
  summarize(gamma_richness = specnumber(abundance)) %>%
  ungroup() %>%
  group_by(site_proj_comm, trt_type2, treatment) %>%
  summarize(gamma_richness = mean(gamma_richness)) %>%
  ungroup()

whittaker_beta <- alpha_richness_df %>%
  full_join(gamma_richness_df, by=c("site_proj_comm", "trt_type2", "treatment")) %>%
  mutate(whittaker_beta = gamma_richness/alpha_richness)

whittaker_control <- whittaker_beta %>%
  filter(trt_type2=="control") %>%
  dplyr::select(-trt_type2, -treatment, -alpha_richness, -gamma_richness) %>%
  rename(whittaker_beta_control = whittaker_beta)

whittaker_rr <- whittaker_beta %>%
  filter(trt_type2 != "control") %>%
  dplyr::select(-alpha_richness, -gamma_richness) %>%
  rename(whittaker_beta_trt = whittaker_beta) %>%
  full_join(whittaker_control, by=c("site_proj_comm")) %>%
  mutate(whittaker_rr = log(whittaker_beta_trt/whittaker_beta_control)) %>%
  filter(site_proj_comm != "ORNL_FACE_0")

write.csv(whittaker_rr, file="whittaker response ratios_2019Nov25.csv", row.names=F)
###
### Read in and process other RR values
###

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


test <- alpha_richness_df %>%
  full_join(synch_metrics_sub, by=c("site_proj_comm", "trt_type2"))

view(alpha_richness_df %>%
  filter(trt_type2!="control"))


view(synch_metrics_sub %>% filter(metric_name=="alpha_stab"))

view(alpha_richness)

metrics_test <- synch_metrics_sub %>%
  filter(metric_name=="alpha_stab") %>%
  select(site_proj_comm, trt_type2, lnRR)

whittaker_test <- alpha_richness_df %>%
  filter(trt_type2 != "control")

full_join <- metrics_test %>%
  full_join(whittaker_test, by=c("site_proj_comm", "trt_type2"))

view(synch_metrics_sub %>%
  filter(site_proj_comm == "Alberta_CCD_0" & trt_type2 == "NR"))

###
### Plot response ratios for pop synch, spatial sync, and gamma stability against whittaker beta
###
response_df <- read.csv("C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\GCD asynchrony\\data\\synchrony_community_withWhittaker_SEMdata.csv")

ggplot(response_df, aes(x=spatial_synch, y=dispersion_diff)) +
  geom_point()

pop_synch_whitt_plot <- ggplot(response_df, aes(x=pop_synch, y=whittaker_beta)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw()

spatial_synch_whitt_plot <- ggplot(response_df, aes(x=spatial_synch, y=whittaker_beta)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw()

gamma_stab_whitt_plot <- ggplot(response_df, aes(x=gamma_stab, y=whittaker_beta)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw()

ggsave(pop_synch_whitt_plot, 
       filename="C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\GCD asynchrony\\figures\\pop sync vs whittaker beta.png",
       height = 3, width = 3, dpi=300)
ggsave(spatial_synch_whitt_plot, 
       filename="C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\GCD asynchrony\\figures\\spatial sync vs whittaker beta.png",
       height = 3, width = 3, dpi=300)
ggsave(gamma_stab_whitt_plot, 
       filename="C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\GCD asynchrony\\figures\\gamma stability vs whittaker beta.png",
       height = 3, width = 3, dpi=300)

### Regression models of whittaker beta vs pop spatial synchrony and gamma stability
summary(lm(whittaker_beta ~ gamma_stab, data=response_df))
summary(lm(whittaker_beta ~ spatial_synch, data=response_df))
summary(lm(whittaker_beta ~ pop_synch, data=response_df))
summary(lm(dispersion_diff ~ gamma_stab, data=response_df))
