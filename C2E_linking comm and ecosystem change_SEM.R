library(tidyverse)
library(ggplot2)
library(grid)
library(PerformanceAnalytics)
library(lavaan)
library(lavaan.survey)

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')



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
  na.omit()%>%
  #transform to improve normality
  mutate(mean_change_transform=sqrt(mean_change), PCEvendiff_transform=log(PCEvendiff+1-min(PCEvendiff)),anpp_PC_transform=log(anpp_PC+1-min(anpp_PC)))
#all of these compare treatment to control!

###exploratory correlations and histograms (all variables compare treatment to control plots)
hist(sqrt(correSEMdata$mean_change)) #mean change  -- slightly left skewed; sqrt improves normality
hist(correSEMdata$PCSdiff) #richness percent change
hist(log(correSEMdata$PCEvendiff+1-min(correSEMdata$PCEvendiff))) #evenness percent change  -- slightly left skewed; improved by log following data translation
hist(correSEMdata$MRSc_diff) #mean rank shift change
hist(correSEMdata$spdiffc) #species number difference
hist(log(correSEMdata$anpp_PC+1-min(correSEMdata$anpp_PC))) #anpp percent change  -- slightly left skewed
dataVis <- correSEMdata%>%
  select(anpp_PC_transform, mean_change_transform, PCSdiff, PCEvendiff_transform, MRSc_diff, spdiffc) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)


#Model 1 - force everything through mean change
correModel1   <-  'anpp_PC_transform ~ mean_change_transform 
               mean_change_transform ~ PCSdiff + spdiffc + PCEvendiff_transform + MRSc_diff'
correModel1Fit <- sem(correModel1, data=correSEMdata, meanstructure=TRUE)
# adjusting for nested design
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=correSEMdata)
survey1Fit <- lavaan.survey(lavaan.fit=correModel1Fit, survey.design=surveyDesign)
survey1Fit  #gives chi-square
summary(survey1Fit)
modificationIndices1 <- modindices(survey1Fit) #modification indices
print(modificationIndices1[modificationIndices1$op == "~",])
print(modificationIndices1[modificationIndices1$op == "~~",])

#Model 2 - remove mean change, and look at direct effects of community metrics on ANPP change
correModel2   <-  'anpp_PC_transform ~ PCSdiff + spdiffc + PCEvendiff_transform + MRSc_diff'
correModel2Fit <- sem(correModel2, data=correSEMdata, meanstructure=TRUE)
# adjusting for nested design
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=correSEMdata)
survey2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=surveyDesign)
survey2Fit  #gives chi-square
summary(survey2Fit)
modificationIndices2 <- modindices(survey2Fit) #modification indices
print(modificationIndices2[modificationIndices2$op == "~",])
print(modificationIndices2[modificationIndices2$op == "~~",])
