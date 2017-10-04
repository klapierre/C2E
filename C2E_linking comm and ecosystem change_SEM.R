library(tidyverse)
library(ggplot2)
library(grid)
library(PerformanceAnalytics)

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')


####data input
#community data
correCommChange <- read.csv('CORRE_ContTreat_Compare_OCT2017.csv')%>%
  select(-X)
#anpp data
correANPPchange <- read.csv('ForBayesianAnalysisANPP_Oct2017.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(-X)
#merge
correSEMdata <- correANPPchange%>%
  left_join(correCommChange)%>%
  na.omit()

###exploratory correlations and histograms (all variables compare treatment to control plots)
hist(sqrt(correSEMdata$mean_change)) #mean change  -- slightly left skewed; sqrt improves normality
hist(correSEMdata$PCSdiff) #richness percent change
hist(log(correSEMdata$PCEvendiff+1-min(correSEMdata$PCEvendiff))) #evenness percent change  -- slightly left skewed; improved by log following data translation
hist(correSEMdata$MRSc_diff) #mean rank shift change
hist(correSEMdata$spdiffc) #species number difference
hist(log(correSEMdata$anpp_PC+1-min(correSEMdata$anpp_PC))) #anpp percent change  -- slightly left skewed
dataVis <- correSEMdata%>%
  mutate(mean_change_transform=sqrt(mean_change), PCEvendiff_transform=log(PCEvendiff+1-min(PCEvendiff)), anpp_PC_transform=log(anpp_PC+1-min(anpp_PC)))%>%
  select(anpp_PC_transform, mean_change_transform, PCSdiff, PCEvendiff_transform, MRSc_diff, spdiffc) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)


