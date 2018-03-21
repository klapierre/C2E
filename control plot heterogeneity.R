library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\C2E\\Products\\CommunityChange\\March2018 WG')



#experiment information
expInfo <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  mutate(exp_year=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))

#diversity data
div <- read.csv('DiversityMetrics_Nov2017.csv')%>%
  separate(exp_year, c('site_code', 'project_name','community_type', 'calendar_year'), sep='::')%>%
  mutate(calendar_year=as.integer(calendar_year))%>%
  left_join(expInfo)%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  #create e^H metric
  mutate(expH=exp(H))

#subset out controls dispersion (as a metric of general heterogeneity at a site)
divControls <- subset(div, subset=(plot_mani==0))%>%
  select(exp_year, dispersion, expH, S, SimpEven, calendar_year, treatment_year)%>%
  separate(exp_year, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(site_project_comm, dispersion, treatment_year)%>%
  group_by(site_project_comm)%>%
  summarise(ctl_dispersion=mean(dispersion))%>%
  ungroup()

#merge control dispersion with bayesian output and determine which treatments that only have a significant intercept (i.e., no significant slopes) but that intercept is biologically meaningful (greater intercept than dispersion among control plots at the site)
bayesian <- read.csv('plot mani_equations_expinteractions_20yr_01172018.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  #join control dispersion data
  left_join(divControls)%>%
  #subset bayesian output to only look at compositional difference
  filter(variable=='mean')%>%
  select(site_project_comm, treatment, intercept, linear, quadratic, ctl_dispersion)%>%
  #generate a column to filter out anything that did not have a significant intercept, linear, or quadratic slope (keep), and generate a column to distinguish those experiments that had a significant linear or quadratic slope from those that only had a signficiant intercept
  mutate(keep=(intercept+linear+quadratic), slopes=(linear+quadratic))%>%
  #filter out anything that did not have a significant intercept, linear, or quadratic slope
  filter(keep!=0)%>%
  #back transform the intercepts by multiplying by the standard deviation and adding the mean compositional difference across all treatments
  mutate(intercept_trans=(intercept*0.165778+0.3160571))%>%
  #determine which are biologically meaningful intercepts for those experiments where intercept is the only factor
  mutate(intercept_bio=intercept_trans-ctl_dispersion)%>%
  filter(slopes>0|intercept_bio>0)

write.csv(bayesian, 'treatments_sig mult diff.csv', row.names=F)
