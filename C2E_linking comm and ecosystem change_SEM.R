library(tidyverse)
library(ggplot2)
library(grid)
library(PerformanceAnalytics)
library(lavaan)
library(lavaan.survey)
library(semPlot)

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

# #######OPTION 1: use all plots, looking at change from plot in yr 1 to yr x through time. compare treatment and control plot changes (vertical arrows); must group by site_project_comm!---------------------------------------
####data input
#treatment data
correTrt <- read.csv('ExperimentInformation_May2017.csv')%>%
  select(site_code, project_name, community_type, treatment_year, treatment, plot_mani)
#community data
correCommChange <- read.csv('CORRE_RAC_Metrics_Oct2017_allReplicates.csv')%>%
  select(-X)
#anpp data
correANPPchange <- read.csv('ANPP_Oct2017.csv')%>%
  #need to fix data input issues that result in same plot having multiple anpp values
  group_by(site_code, project_name, community_type, treatment, treatment_year, calendar_year, plot_id)%>%
  summarise(anpp=mean(anpp))%>%
  ungroup()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))
#calculate experiment length
expLength <- correANPPchange%>%
  group_by(site_project_comm)%>%
  mutate(experiment_length=max(treatment_year))%>%
  ungroup()%>%
  select(site_project_comm, experiment_length)%>%
  unique() #mean experiment length is 10.36 across ANPP data
#merge
correSEMdata <- correANPPchange%>%
  left_join(correCommChange)%>%
  left_join(expLength)%>%
  na.omit()%>%
  #remove CDR e001 and e002 intermediate treatments (before CDR e001/e002 make up 33.8% of the data; after removal drops to 3 trts todal, to make up only 18.5% of total data)
  mutate(trt_drop=ifelse(project_name=='e001'&treatment==2, 1, ifelse(project_name=='e001'&treatment==3, 1, ifelse(project_name=='e001'&treatment==4, 1, ifelse(project_name=='e001'&treatment==5, 1, ifelse(project_name=='e001'&treatment==7, 1, ifelse(project_name=='e002'&treatment=='2_f_u_n', 1, ifelse(project_name=='e002'&treatment=='3_f_u_n', 1, ifelse(project_name=='e002'&treatment=='4_f_u_n', 1, ifelse(project_name=='e002'&treatment=='5_f_u_n', 1, ifelse(project_name=='e002'&treatment=='7_f_u_n', 1, 0)))))))))))%>%
  filter(trt_drop==0)%>%
  #transform to improve normality
  mutate(mean_change_transform=sqrt(bc_dissim), E_q_transform=log(E_q), S_transform=sqrt(S), anpp_change_transform=log(anpp))%>%
  #make a binary treatment variable
  left_join(correTrt)%>%
  mutate(treatment_binary=ifelse(plot_mani==0, 0, 1))

###exploratory correlations and histograms (all variables compare treatment to control plots)
dataVis <- correSEMdata%>%
  select(anpp_change_transform, mean_change_transform, appearance, disappearance, E_q_transform, MRSc, S_transform) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)


####3 models: all data (across all years); year 2 vs year 3 (29 experiments, must run 2 separate SEMs and compare); year 2 vx year 10 (12 experiments, must run 2 separate SEMs and compare)

#note that appearance/disappearance are count data, and must use a poisson model (need to figure out how to do this)
#model 1 - all data
correModel <- '
anpp_change_transform ~ mean_change_transform + treatment_binary
mean_change_transform ~ appearance + disappearance + E_q_transform + MRSc + S_transform
appearance ~ treatment_binary
disappearance ~ treatment_binary
E_q_transform ~ treatment_binary
MRSc ~ treatment_binary
S_transform ~ treatment_binary

#covariances
appearance~~disappearance
appearance~~E_q_transform
appearance~~MRSc
appearance~~S_transform
disappearance~~E_q_transform
disappearance~~MRSc
disappearance~~S_transform
E_q_transform~~MRSc
E_q_transform~~S_transform
MRSc~~S_transform'

correModelFit <- sem(correModel, data=correSEMdata, meanstructure=TRUE)
# adjusting for nested design
design <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata)
modelFit <- lavaan.survey(lavaan.fit=correModelFit, survey.design=design)
summary(modelFit, standardize=T)
modificationIndices <- modindices(modelFit, sort.=TRUE, minimum.value=3) #modification indices
print(modificationIndices[modificationIndices$op == "~",])
print(modificationIndices[modificationIndices$op == "~~",])
semPaths(modelFit, 'Standard', 'Estimates', layout='tree3', residuals=F, intercepts=F)



correSEMyears <- correSEMdata%>%
  #gather to combine metric with year
  select(site_code, project_name, community_type, site_project_comm, experiment_length, treatment, treatment_binary, treatment_year, calendar_year, plot_id, anpp_change_transform, mean_change_transform, S_transform, E_q_transform, appearance, disappearance, MRSc)%>%
  gather(key=metric, value=value, anpp_change_transform:MRSc, na.rm=F)%>%
  mutate(metric_year=paste(metric, treatment_year, sep='_'))%>%
  na.omit()

#model 2 - years 2 vs 3
correSEMdata3 <- correSEMyears%>%
  filter(treatment_year==2|treatment_year==10|treatment_year==6)%>%
  select(-metric, -treatment_year, -calendar_year)%>%
  spread(key=metric_year, value=value)%>%
  na.omit()

correModel2 <- '
#year 2
anpp_change_transform_2 ~ mean_change_transform_2 + treatment_binary
mean_change_transform_2 ~ appearance_2 + disappearance_2 + E_q_transform_2 + MRSc_2 + S_transform_2
appearance_2 ~ treatment_binary
disappearance_2 ~ treatment_binary
E_q_transform_2 ~ treatment_binary
MRSc_2 ~ treatment_binary
S_transform_2 ~ treatment_binary

#covariances
appearance_2~~disappearance_2
appearance_2~~E_q_transform_2
appearance_2~~MRSc_2
appearance_2~~S_transform_2
disappearance_2~~E_q_transform_2
disappearance_2~~MRSc_2
disappearance_2~~S_transform_2
E_q_transform_2~~MRSc_2
E_q_transform_2~~S_transform_2
MRSc_2~~S_transform_2


#year 5
anpp_change_transform_6 ~ mean_change_transform_6 + treatment_binary
mean_change_transform_6 ~ appearance_6 + disappearance_6 + E_q_transform_6 + MRSc_6 + S_transform_6
appearance_6 ~ treatment_binary
disappearance_6 ~ treatment_binary
E_q_transform_6 ~ treatment_binary
MRSc_6 ~ treatment_binary
S_transform_6 ~ treatment_binary

#covariances
appearance_6~~disappearance_6
appearance_6~~E_q_transform_6
appearance_6~~MRSc_6
appearance_6~~S_transform_6
disappearance_6~~E_q_transform_6
disappearance_6~~MRSc_6
disappearance_6~~S_transform_6
E_q_transform_6~~MRSc_6
E_q_transform_6~~S_transform_6
MRSc_6~~S_transform_6


#year 10
anpp_change_transform_10 ~ mean_change_transform_10 + treatment_binary
mean_change_transform_10 ~ appearance_10 + disappearance_10 + E_q_transform_10 + MRSc_10 + S_transform_10
appearance_10 ~ treatment_binary
disappearance_10 ~ treatment_binary
E_q_transform_10 ~ treatment_binary
MRSc_10 ~ treatment_binary
S_transform_10 ~ treatment_binary

#covariances
appearance_10~~disappearance_10
appearance_10~~E_q_transform_10
appearance_10~~MRSc_10
appearance_10~~S_transform_10
disappearance_10~~E_q_transform_10
disappearance_10~~MRSc_10
disappearance_10~~S_transform_10
E_q_transform_10~~MRSc_10
E_q_transform_10~~S_transform_10
MRSc_10~~S_transform_10


#temporal covariances
mean_change_transform_2~~mean_change_transform_6
anpp_change_transform_2~~anpp_change_transform_6
appearance_2~~appearance_6
disappearance_2~~disappearance_6
E_q_transform_2~~E_q_transform_6
MRSc_2~~MRSc_6
S_transform_2~~S_transform_6
mean_change_transform_2~~mean_change_transform_10
anpp_change_transform_2~~anpp_change_transform_10
appearance_2~~appearance_10
disappearance_2~~disappearance_10
E_q_transform_2~~E_q_transform_10
MRSc_2~~MRSc_10
S_transform_2~~S_transform_10
mean_change_transform_6~~mean_change_transform_10
anpp_change_transform_6~~anpp_change_transform_10
appearance_6~~appearance_10
disappearance_6~~disappearance_10
E_q_transform_6~~E_q_transform_10
MRSc_6~~MRSc_10
S_transform_6~~S_transform_10
'

correModel2Fit <- sem(correModel2, data=correSEMdata3, meanstructure=TRUE)
# adjusting for nested design
design <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata3)
model2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=design)
summary(model2Fit, standardize=T)
modificationIndices2 <- modindices(model2Fit, sort.=TRUE, minimum.value=3) #modification indices
print(modificationIndices2[modificationIndices2$op == "~",])
print(modificationIndices2[modificationIndices2$op == "~~",])
semPaths(model2Fit, 'Standardized', 'Estimates', layout='tree3', intercepts=F)









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
  select(-X)%>%
  group_by(site_code, project_name, treatment_year, calendar_year, treatment, community_type, site_project_comm)%>%
  summarise(anpp_change=mean(anpp))%>%
  ungroup()
#calculate experiment length
expLength <- correANPPchange%>%
  group_by(site_project_comm)%>%
  mutate(experiment_length=max(treatment_year))%>%
  ungroup()%>%
  select(site_project_comm, experiment_length)%>%
  unique() #mean experiment length is 10.36 across ANPP data

#merge
correSEMdata <- correANPPchange%>%
  left_join(correCommChange)%>%
  left_join(expLength)%>%
  na.omit()%>%
  #remove CDR e001 and e002 intermediate treatments (before CDR e001/e002 make up 33.8% of the data; after removal drops to 3 trts todal, to make up only 18.5% of total data)
  mutate(trt_drop=ifelse(project_name=='e001'&treatment==2, 1, ifelse(project_name=='e001'&treatment==3, 1, ifelse(project_name=='e001'&treatment==4, 1, ifelse(project_name=='e001'&treatment==5, 1, ifelse(project_name=='e001'&treatment==7, 1, ifelse(project_name=='e002'&treatment=='2_f_u_n', 1, ifelse(project_name=='e002'&treatment=='3_f_u_n', 1, ifelse(project_name=='e002'&treatment=='4_f_u_n', 1, ifelse(project_name=='e002'&treatment=='5_f_u_n', 1, ifelse(project_name=='e002'&treatment=='7_f_u_n', 1, 0)))))))))))%>%
  filter(trt_drop==0)%>%
  #transform to improve normality
  mutate(mean_change_transform=sqrt(mean_change), even_transform=log(even), S_transform=sqrt(S), anpp_change_transform=sqrt(anpp_change))%>%
  #make a site_code:community_type column
  mutate(site_comm=paste(site_code, community_type, sep='::'))%>%
  #make a binary treatment variable
  left_join(correTrt)%>%
  mutate(treatment_binary=ifelse(plot_mani==0, 0, 1))%>%
  #remove non-transformed variables
  select(site_code, project_name, community_type, site_project_comm, site_comm, treatment, treatment_binary, experiment_length, treatment_year, anpp_change_transform, mean_change_transform, gain, loss, S_transform, even_transform, MRSc)

###exploratory correlations and histograms (all variables compare treatment to control plots)
dataVis <- correSEMdata%>%
  select(anpp_change_transform, mean_change_transform, gain, loss, even_transform, MRSc, S_transform) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)

correSEMyears <- correSEMdata%>%
  #gather to combine metric with year
  gather(key=metric, value=value, anpp_change_transform:MRSc, na.rm=F)%>%
  mutate(metric_year=paste(metric, treatment_year, sep='_'))


####3 models: all data (across all years); year 2 vs year 3 (29 experiments, must run 2 separate SEMs and compare); year 2 vx year 10 (12 experiments, must run 2 separate SEMs and compare)

#note that appearance/disappearance are count data, and must use a poisson model (need to figure out how to do this)
#model 1 - all data
correModel <- '
anpp_change_transform ~ mean_change_transform + treatment_binary
mean_change_transform ~ gain + loss + even_transform + MRSc + S_transform
# physiology_latent =~ treatment_binary
gain ~ treatment_binary
loss ~ treatment_binary
even_transform ~ treatment_binary
MRSc ~ treatment_binary
S_transform ~ treatment_binary

#covariances
gain~~loss
gain~~even_transform
gain~~MRSc
gain~~S_transform
loss~~even_transform
loss~~MRSc
loss~~S_transform
even_transform~~MRSc
even_transform~~S_transform
MRSc~~S_transform'

correModelFit <- sem(correModel, data=correSEMdata, meanstructure=TRUE)
# adjusting for nested design
design <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata)
modelFit <- lavaan.survey(lavaan.fit=correModelFit, survey.design=design)
summary(modelFit, standardize=T)
modificationIndices <- modindices(modelFit, sort.=TRUE, minimum.value=3) #modification indices
print(modificationIndices[modificationIndices$op == "~",])
print(modificationIndices[modificationIndices$op == "~~",])
semPaths(modelFit, 'Standard', 'Estimates', layout='tree3', residuals=F, intercepts=F)




#model 2 - years 2 vs 3 - nothing matters
correSEMdata3 <- correSEMyears%>%
  filter(treatment_year!=1&treatment_year<4)%>%
  select(-metric, -treatment_year)%>%
  spread(key=metric_year, value=value)%>%
  na.omit()

correModel2 <- '
#year 2
anpp_change_transform_2 ~ mean_change_transform_2 + physiology_latent_2
mean_change_transform_2 ~ gain_2 + loss_2 + even_transform_2 + MRSc_2 + S_transform_2
physiology_latent_2 =~ treatment_binary
gain_2 ~ treatment_binary
loss_2 ~ treatment_binary
even_transform_2 ~ treatment_binary
MRSc_2 ~ treatment_binary
S_transform_2 ~ treatment_binary

#covariances
gain_2~~loss_2
gain_2~~even_transform_2
gain_2~~MRSc_2
gain_2~~S_transform_2
loss_2~~even_transform_2
loss_2~~MRSc_2
loss_2~~S_transform_2
even_transform_2~~MRSc_2
even_transform_2~~S_transform_2
MRSc_2~~S_transform_2'

correModel2Fit <- sem(correModel2, data=correSEMdata3, meanstructure=TRUE)
# adjusting for nested design
design <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata3)
model2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=design)
summary(model2Fit, standardize=T)
modificationIndices2 <- modindices(model2Fit, sort.=TRUE, minimum.value=3) #modification indices
print(modificationIndices2[modificationIndices2$op == "~",])
print(modificationIndices2[modificationIndices2$op == "~~",])
semPaths(model2Fit, 'Standardized', 'Estimates', layout='tree3')

correModel3 <- '
#year 3
anpp_change_transform_3 ~ mean_change_transform_3 + physiology_latent_3
mean_change_transform_3 ~ gain_3 + loss_3 + even_transform_3 + MRSc_3 + S_transform_3
physiology_latent_3 =~ treatment_binary
gain_3 ~ treatment_binary
loss_3 ~ treatment_binary
even_transform_3 ~ treatment_binary
MRSc_3 ~ treatment_binary
S_transform_3 ~ treatment_binary

#covariances
gain_3~~loss_3
gain_3~~even_transform_3
gain_3~~MRSc_3
gain_3~~S_transform_3
loss_3~~even_transform_3
loss_3~~MRSc_3
loss_3~~S_transform_3
even_transform_3~~MRSc_3
even_transform_3~~S_transform_3
MRSc_3~~S_transform_3'

correModel3Fit <- sem(correModel3, data=correSEMdata3, meanstructure=TRUE)
# adjusting for nested design
design <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata3)
model3Fit <- lavaan.survey(lavaan.fit=correModel3Fit, survey.design=design)
summary(model3Fit, standardize=T)
modificationIndices3 <- modindices(model3Fit, sort.=TRUE, minimum.value=3) #modification indices
print(modificationIndices3[modificationIndices3$op == "~",])
print(modificationIndices3[modificationIndices3$op == "~~",])
semPaths(model3Fit, 'Standardized', 'Estimates', layout='tree3')


#model 3 - years 2 vs 10
correSEMdata10 <- correSEMyears%>%
  filter(treatment_year==2|treatment_year==10)%>%
  select(-metric, -treatment_year)%>%
  spread(key=metric_year, value=value)%>%
  na.omit()

#year 10 of data - model significant (not good)
correModel10 <- '
#year 10
anpp_change_transform_10 ~ mean_change_transform_10 + physiology_latent_10
mean_change_transform_10 ~ gain_10 + loss_10 + even_transform_10 + MRSc_10 + S_transform_10
physiology_latent_10 =~ treatment_binary
gain_10 ~ treatment_binary
loss_10 ~ treatment_binary
even_transform_10 ~ treatment_binary
MRSc_10 ~ treatment_binary
S_transform_10 ~ treatment_binary

#covariances
gain_10~~loss_10
gain_10~~even_transform_10
gain_10~~MRSc_10
gain_10~~S_transform_10
loss_10~~even_transform_10
loss_10~~MRSc_10
loss_10~~S_transform_10
even_transform_10~~MRSc_10
even_transform_10~~S_transform_10
MRSc_10~~S_transform_10'

correModel10Fit <- sem(correModel10, data=correSEMdata10, meanstructure=TRUE)
# adjusting for nested design
design <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata10)
model10Fit <- lavaan.survey(lavaan.fit=correModel10Fit, survey.design=design)
summary(model10Fit, standardize=T)
modificationIndices10 <- modindices(model10Fit, sort.=TRUE, minimum.value=3) #modification indices
print(modificationIndices10[modificationIndices10$op == "~",])
print(modificationIndices10[modificationIndices10$op == "~~",])
semPaths(model10Fit, 'Standardized', 'Estimates', layout='tree3')


# #Model 2 - remove mean change, and look at direct effects of community metrics on ANPP change
# correModel2 <- 'anpp_transform ~ gain + loss + even_transform + MRSc + S_transform + physiology_latent
#                 physiology_latent =~ treatment
#                 gain ~ treatment
#                 loss ~ treatment
#                 even_transform ~ treatment
#                 MRSc ~ treatment
#                 S_transform ~ treatment'
# correModel2Fit <- sem(correModel2, data=correSEMdata, meanstructure=TRUE)
# # adjusting for nested design
# surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata)
# survey2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=surveyDesign)
# survey2Fit  #gives chi-square
# summary(survey2Fit)
# modificationIndices2 <- modindices(survey2Fit) #modification indices
# print(modificationIndices2[modificationIndices2$op == "~",])
# print(modificationIndices2[modificationIndices2$op == "~~",])





# #######OPTION 3: use change between treatment and control plots in each year (horizontal arrows); can't include treatment explicitly, as it is used to calculate difference at each time point---------------------------------------
# ####data input
# #community data
# correCommChange <- read.csv('CORRE_ContTreat_Compare_OCT2017.csv')%>%
#   select(-X)
# #anpp data
# correANPPchange <- read.csv('ForBayesianAnalysisANPP_Oct2017.csv')%>%
#   mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
#   select(-X)
# #merge
# correSEMdata <- correANPPchange%>%
#   left_join(correCommChange)%>%
#   na.omit()%>%
#   #transform to improve normality
#   mutate(mean_change_transform=sqrt(mean_change), PCEvendiff_transform=log(PCEvendiff+1-min(PCEvendiff)),anpp_PC_transform=log(anpp_PC+1-min(anpp_PC)))
# #all of these compare treatment to control!
# 
# ###exploratory correlations and histograms (all variables compare treatment to control plots)
# dataVis <- correSEMdata%>%
#   select(anpp_PC_transform, mean_change_transform, PCSdiff, PCEvendiff_transform, MRSc_diff, spdiffc) #make visualization dataframe
# chart.Correlation(dataVis, histogram=T, pch=19)
# 
# 
# ##PROBLEM - how to include physiology as latent variable if no path goes into it?
# #Model 1 - force everything through mean change
# correModel1   <-  'anpp_PC_transform ~ mean_change_transform 
# mean_change_transform ~ PCSdiff + spdiffc + PCEvendiff_transform + MRSc_diff'
# correModel1Fit <- sem(correModel1, data=correSEMdata, meanstructure=TRUE)
# # adjusting for nested design
# surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=correSEMdata)
# survey1Fit <- lavaan.survey(lavaan.fit=correModel1Fit, survey.design=surveyDesign)
# survey1Fit  #gives chi-square
# summary(survey1Fit)
# modificationIndices1 <- modindices(survey1Fit) #modification indices
# print(modificationIndices1[modificationIndices1$op == "~",])
# print(modificationIndices1[modificationIndices1$op == "~~",])
# 
# #Model 2 - remove mean change, and look at direct effects of community metrics on ANPP change
# correModel2   <-  'anpp_PC_transform ~ PCSdiff + spdiffc + PCEvendiff_transform + MRSc_diff'
# correModel2Fit <- sem(correModel2, data=correSEMdata, meanstructure=TRUE)
# # adjusting for nested design
# surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=correSEMdata)
# survey2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=surveyDesign)
# survey2Fit  #gives chi-square
# summary(survey2Fit)
# modificationIndices2 <- modindices(survey2Fit) #modification indices
# print(modificationIndices2[modificationIndices2$op == "~",])
# print(modificationIndices2[modificationIndices2$op == "~~",])