library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

#kim's wd
if(exists("diffwd")==FALSE){setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\to yang')}
if(exists("diffwd")){ setwd(diffwd) }

#read in data
projectSummary <- read.csv('CORRE_project_summary.csv')
nExperiments <- read.csv('treatments_nitrogen_experiments.csv')
relAbund <- read.csv('CORRE_relative_abundance.csv')%>%
  left_join(projectSummary)%>%
  left_join(nExperiments, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  filter(experiment_length>6)%>%
  subset(!is.na(n))

#get rid of datasets that have fewer than 7 datapoints
dataLength <- relAbund%>%
  group_by(site_code, project_name, community_type, treatment, plot_id, treatment_year)%>%
  unique()%>%
  summarise(relcov=mean(relcov))%>%
  ungroup()%>%
  select(-relcov)%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  summarise(count=length(treatment_year))

relAbundLength <- relAbund%>%  
  ungroup()%>%
  left_join(dataLength)%>%
  filter(count>7)

#get controls
relAbundCtl <- relAbundLength%>%  
  ungroup()%>%
  filter(n==0)

#get highest N
relAbundN <- relAbundLength%>%
  group_by(site_code, project_name, community_type)%>%
  filter(n==max(n))%>%
  filter(n>0)%>%
  rbind(relAbundCtl)%>%
  ungroup()%>%
  mutate(n_treatment=ifelse(n==0, 0, 1))


