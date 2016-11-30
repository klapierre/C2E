library(plyr)
library(dplyr)
library(tidyr)
library(codyn)
library(ggplot2)
library(grid)

#kim's wd
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\to yang')

source('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\C2E\\C2E\\aggregate across replicates.R')

#calculate dominance: berger-parker method
dominance <- relAbundAgg%>%
  group_by(site_code, project_name, community_type, treatment, treatment_year, calendar_year)%>%
  summarise(bp_dominance=max(relcov_agg))

