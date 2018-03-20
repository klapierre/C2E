library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\C2E\\Products\\CommunityChange\\March2018 WG')

multMetrics <- read.csv('CORRE_Mult_Metrics_March2018.csv')
trt <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  mutate(site_project_comm=as.factor(paste(site_code, project_name, community_type, sep='_')))%>%
  select(-X)

multMetricsTrt <- multMetrics%>%
  left_join(trt)%>%
  select(treatment_year, treatment, composition_change, site_project_comm, plot_mani)
