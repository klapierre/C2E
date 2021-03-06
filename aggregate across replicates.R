library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(tidyr)

#kim's wd
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\to yang')

source('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\C2E\\C2E\\subsetting 7 yr N data.R')

relAbundAgg <- relAbundN%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, treatment_year, calendar_year, genus_species, plot_id, relcov, n, n_treatment)%>%
  group_by(site_code, project_name, community_type, treatment, treatment_year, calendar_year, genus_species, n, n_treatment)%>%
  summarise(relcov_agg=mean(relcov))