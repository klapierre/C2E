library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\C2E\\Products\\CommunityChange\\March2018 WG')

ttest <- read.csv('sig_ttest.csv')
regression <- read.csv('treatments_sig regression.csv')
bayesian_diff <- read.csv('treatments_sig mult diff.csv')%>%
  mutate(bayesian_intercept=intercept, bayesian_linear=linear, bayesian_quadratic=quadratic, bayesian=1)%>%
  select(site_project_comm, treatment, bayesian_intercept, bayesian_linear, bayesian_quadratic, bayesian)
permanova <- read.csv('permanova out.csv')%>%
  mutate(permanova=1)


comparison <- ttest%>%
  full_join(regression)%>%
  full_join(bayesian_diff)%>%
  full_join(permanova)%>%
  select(site_project_comm, treatment, ttest, regression, permanova, bayesian)%>%
  replace(is.na(.), 0)%>%
  mutate(num_tests=ttest+regression+permanova+bayesian, perm_bayes=permanova+bayesian)
