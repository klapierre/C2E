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

#######OPTION 1: use all plots, looking at change from plot in yr 1 to yr x through time. compare treatment and control plot changes (vertical arrows); must group by site_project_comm!---------------------------------------
####data input
#treatment data
correTrt <- read.csv('ExperimentInformation_May2017.csv')%>%
  select(site_code, project_name, community_type, treatment_year, treatment, plot_mani)
#community data
correCommChange <- read.csv('CORRE_RAC_Metrics_Oct2017_allReplicates.csv')%>%
  select(-X)
#anpp data
correANPPchange <- read.csv('ANPP_Oct2017.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(-X)
#merge
correSEMdata <- correANPPchange%>%
  left_join(correCommChange)%>%
  na.omit()%>%
  #transform to improve normality
  mutate(mean_change_transform=sqrt(bc_dissim), E_q_transform=log(E_q), S_transform=sqrt(S), anpp_transform=log(anpp))%>%
  #make a binary treatment variable
  left_join(correTrt)%>%
  select(-treatment)%>%
  mutate(treatment=ifelse(plot_mani==0, 0, 1))

###exploratory correlations and histograms (all variables compare treatment to control plots)
dataVis <- correSEMdata%>%
  select(anpp_transform, mean_change_transform, appearance, disappearance, E_q_transform, MRSc, S_transform) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)


#Model 1 - force everything through mean change
#note that appearance/disappearance are count data, and must use a poisson model (need to figure out how to do this)
correModel1 <- 'anpp_transform ~ mean_change_transform + physiology_latent
                mean_change_transform ~ appearance + disappearance + E_q_transform + MRSc + S_transform
                physiology_latent =~ treatment
                appearance ~ treatment
                disappearance ~ treatment
                E_q_transform ~ treatment
                MRSc ~ treatment
                S_transform ~ treatment
                '
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
correModel2 <- 'anpp_transform ~ appearance + disappearance + E_q_transform + MRSc + S_transform + physiology_latent
                physiology_latent =~ treatment
                appearance ~ treatment
                disappearance ~ treatment
                E_q_transform ~ treatment
                MRSc ~ treatment
                S_transform ~ treatment'
correModel2Fit <- sem(correModel2, data=correSEMdata, meanstructure=TRUE)
# adjusting for nested design
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=correSEMdata)
survey2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=surveyDesign)
survey2Fit  #gives chi-square
summary(survey2Fit)
modificationIndices2 <- modindices(survey2Fit) #modification indices
print(modificationIndices2[modificationIndices2$op == "~",])
print(modificationIndices2[modificationIndices2$op == "~~",])




#######OPTION 2: use mean across plots within a treatment, looking at change from plot in yr 1 to yr x through time. compare treatment and control plot changes (vertical arrows)---------------------------------------
####data input
#treatment data
correTrt <- read.csv('ExperimentInformation_May2017.csv')%>%
  select(site_code, project_name, community_type, treatment_year, treatment, plot_mani)
#community data
correCommChange <- read.csv('CORRE_RAC_Metrics_Oct2017_allyears.csv')%>%
  select(-X)
#anpp data
correANPPchange <- read.csv('ANPP_Oct2017.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(-X)
#merge
correSEMdata <- correANPPchange%>%
  left_join(correCommChange)%>%
  na.omit()%>%
  #transform to improve normality
  mutate(mean_change_transform=sqrt(mean_change), even_transform=log(even), S_transform=sqrt(S), anpp_transform=sqrt(anpp))%>%
  #make a binary treatment variable
  left_join(correTrt)%>%
  select(-treatment)%>%
  mutate(treatment=ifelse(plot_mani==0, 0, 1))

###exploratory correlations and histograms (all variables compare treatment to control plots)
dataVis <- correSEMdata%>%
  select(anpp_transform, mean_change_transform, gain, loss, even_transform, MRSc, S_transform) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)


#Model 1 - force everything through mean change
#note that appearance/disappearance are count data, and must use a poisson model (need to figure out how to do this)
correModel1 <- 'anpp_transform ~ mean_change_transform + physiology_latent
                mean_change_transform ~ gain + loss + even_transform + MRSc + S_transform
                physiology_latent =~ treatment
                gain ~ treatment
                loss ~ treatment
                even_transform ~ treatment
                MRSc ~ treatment
                S_transform ~ treatment
                '
correModel1Fit <- sem(correModel1, data=correSEMdata, meanstructure=TRUE)
# adjusting for nested design
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata)
survey1Fit <- lavaan.survey(lavaan.fit=correModel1Fit, survey.design=surveyDesign)
survey1Fit  #gives chi-square
summary(survey1Fit)
modificationIndices1 <- modindices(survey1Fit) #modification indices
print(modificationIndices1[modificationIndices1$op == "~",])
print(modificationIndices1[modificationIndices1$op == "~~",])

#Model 2 - remove mean change, and look at direct effects of community metrics on ANPP change
correModel2 <- 'anpp_transform ~ gain + loss + even_transform + MRSc + S_transform + physiology_latent
                physiology_latent =~ treatment
                gain ~ treatment
                loss ~ treatment
                even_transform ~ treatment
                MRSc ~ treatment
                S_transform ~ treatment'
correModel2Fit <- sem(correModel2, data=correSEMdata, meanstructure=TRUE)
# adjusting for nested design
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata)
survey2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=surveyDesign)
survey2Fit  #gives chi-square
summary(survey2Fit)
modificationIndices2 <- modindices(survey2Fit) #modification indices
print(modificationIndices2[modificationIndices2$op == "~",])
print(modificationIndices2[modificationIndices2$op == "~~",])





#######OPTION 3: use change between treatment and control plots in each year (horizontal arrows); can't include treatment explicitly, as it is used to calculate difference at each time point---------------------------------------
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
dataVis <- correSEMdata%>%
  select(anpp_PC_transform, mean_change_transform, PCSdiff, PCEvendiff_transform, MRSc_diff, spdiffc) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)


##PROBLEM - how to include physiology as latent variable if no path goes into it?
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