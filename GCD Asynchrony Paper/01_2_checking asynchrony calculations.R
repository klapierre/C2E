## Checking stability calculations
## April 22, 2020

setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\GCD asynchrony\\")
cover_df <- read.csv('data\\SpeciesRawAbundance_March2019.csv', header = T)
exp_info <- read.csv('data\\ExperimentInformation_March2019.csv', header = T) %>%
  mutate(site_proj_comm_trt = paste(site_code, project_name, community_type, treatment, sep='_'))
site_info <- read.csv('data\\SiteExperimentDetails_March2019.csv')
exp_info_short <- unique(dplyr::select(exp_info, -X, -calendar_year, -treatment_year))

full_df <- cover_df %>%
  full_join(exp_info, by=c("site_code", 
                           "project_name", 
                           "calendar_year",
                           "treatment_year",
                           "treatment",
                           "community_type"))

pplots_df <- full_df %>% filter(project_name=="pplots")


gamma_stability_pplots <- pplots_df %>%
  group_by(treatment, plot_id, calendar_year) %>%
  summarise(total_abund=sum(abundance)) %>%
  group_by(treatment) %>%
  summarise(mean_total_abund = mean(total_abund),
            sd_total_abund = sd(total_abund),
            gamma_stability = mean(total_abund)/sd(total_abund))


metacomm_stability <- species_cover %>%
  gather("species","cover", starts_with("sp")) %>%
  group_by(site_proj_comm, calendar_year) %>%
  summarise(total_abund = sum(cover)) %>%
  group_by(site_proj_comm) %>%
  summarise(gamma_stability = mean(total_abund)/sd(total_abund))

