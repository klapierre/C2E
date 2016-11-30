library(plyr)
library(dplyr)
library(tidyr)
library(codyn)
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

#get controls
relAbundCtl <- relAbund%>%
  filter(n==0)

#get highest N
relAbundN <- relAbund%>%
  group_by(site_code, project_name, community_type)%>%
  filter(n==max(n))%>%
  filter(n>0)%>%
  rbind(relAbundCtl)

