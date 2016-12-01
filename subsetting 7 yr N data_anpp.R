library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(tidyr)

#kim's wd
if(exists("diffwd")==FALSE){setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\to yang')}
if(exists("diffwd")){ setwd(diffwd) }

#read in data
projectSummary <- read.csv('CORRE_project_summary.csv')
nExperiments <- read.csv('treatments_nitrogen_experiments.csv')
relAbund <- read.csv('CORRE_anpp_raw.csv')%>%
  left_join(projectSummary)%>%
  left_join(nExperiments, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  filter(experiment_length>6)%>%
  subset(!is.na(n))

#get rid of datasets that have fewer than 7 datapoints
dataLength <- relAbund%>%
  group_by(site_code, project_name, community_type, treatment, plot_id, calendar_year)%>%
  unique()%>%
  summarise(anpp=mean(anpp))%>%
  ungroup()%>%
  select(-anpp)%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  summarise(count=length(calendar_year))

relAbundLength <- relAbund%>%  
  ungroup()%>%
  left_join(dataLength)%>%
  filter(count>2)

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
