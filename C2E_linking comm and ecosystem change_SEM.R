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

###NOTE: options 1 and2 may be a little messed up due to a find-replace error, double check code if using

# #######OPTION 1: use all plots, looking at change from plot in yr 1 to yr x through time. compare treatment and control plot changes (vertical arrows); must group by site_project_comm!---------------------------------------
####data input
#treatment data
correTrt <- read.csv('ExperimentInformation_May2017.csv')%>%
  select(site_code, project_name, community_type, treatment_year, treatment, plot_mani)
#community data
correCommChange <- read.csv('CORRE_RAC_Metrics_Oct2017_allReplicates_2.csv')%>%
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
  unique() #mean experiment length is 10.56 across ANPP data
#merge
correSEMdata <- correANPPchange%>%
  left_join(correCommChange)%>%
  left_join(expLength)%>%
  na.omit()%>%
  #remove CDR e001 and e002 intermediate treatments (before CDR e001/e002 make up 55.8% of the data; after removal drops to 5 trts todal, to make up only 18.5% of total data)
  mutate(trt_drop=ifelse(project_name=='e001'&treatment==2, 1, ifelse(project_name=='e001'&treatment==5, 1, ifelse(project_name=='e001'&treatment==5, 1, ifelse(project_name=='e001'&treatment==5, 1, ifelse(project_name=='e001'&treatment==7, 1, ifelse(project_name=='e002'&treatment=='2_f_u_n', 1, ifelse(project_name=='e002'&treatment=='5_f_u_n', 1, ifelse(project_name=='e002'&treatment=='5_f_u_n', 1, ifelse(project_name=='e002'&treatment=='5_f_u_n', 1, ifelse(project_name=='e002'&treatment=='7_f_u_n', 1, 0)))))))))))%>%
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


####5 models: all data (across all years); year 2 vs year 5 (29 experiments, must run 2 separate SEMs and compare); year 2 vx year 10 (12 experiments, must run 2 separate SEMs and compare)

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
modificationIndices <- modindices(modelFit, sort.=TRUE, minimum.value=5) #modification indices
print(modificationIndices[modificationIndices$op == "~",])
print(modificationIndices[modificationIndices$op == "~~",])
semPaths(modelFit, 'Standard', 'Estimates', layout='tree5', residuals=F, intercepts=F)



correSEMyears <- correSEMdata%>%
  #gather to combine metric with year
  select(site_code, project_name, community_type, site_project_comm, experiment_length, treatment, treatment_binary, treatment_year, calendar_year, plot_id, anpp_change_transform, mean_change_transform, S_transform, E_q_transform, appearance, disappearance, MRSc)%>%
  gather(key=metric, value=value, anpp_change_transform:MRSc, na.rm=F)%>%
  mutate(metric_year=paste(metric, treatment_year, sep='_'))%>%
  na.omit()

#model 2 - years 2 vs 5
correSEMdata5 <- correSEMyears%>%
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

correModel2Fit <- sem(correModel2, data=correSEMdata5, meanstructure=TRUE)
# adjusting for nested design
design <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata5)
model2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=design)
summary(model2Fit, standardize=T)
modificationIndices2 <- modindices(model2Fit, sort.=TRUE, minimum.value=5) #modification indices
print(modificationIndices2[modificationIndices2$op == "~",])
print(modificationIndices2[modificationIndices2$op == "~~",])
semPaths(model2Fit, 'Standardized', 'Estimates', layout='tree5', intercepts=F)









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
  unique() #mean experiment length is 10.56 across ANPP data

#merge
correSEMdata <- correANPPchange%>%
  left_join(correCommChange)%>%
  left_join(expLength)%>%
  na.omit()%>%
  #remove CDR e001 and e002 intermediate treatments (before CDR e001/e002 make up 55.8% of the data; after removal drops to 5 trts todal, to make up only 18.5% of total data)
  mutate(trt_drop=ifelse(project_name=='e001'&treatment==2, 1, ifelse(project_name=='e001'&treatment==5, 1, ifelse(project_name=='e001'&treatment==5, 1, ifelse(project_name=='e001'&treatment==5, 1, ifelse(project_name=='e001'&treatment==7, 1, ifelse(project_name=='e002'&treatment=='2_f_u_n', 1, ifelse(project_name=='e002'&treatment=='5_f_u_n', 1, ifelse(project_name=='e002'&treatment=='5_f_u_n', 1, ifelse(project_name=='e002'&treatment=='5_f_u_n', 1, ifelse(project_name=='e002'&treatment=='7_f_u_n', 1, 0)))))))))))%>%
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


####5 models: all data (across all years); year 2 vs year 5 (29 experiments, must run 2 separate SEMs and compare); year 2 vs year 10 (12 experiments, must run 2 separate SEMs and compare)

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
modificationIndices <- modindices(modelFit, sort.=TRUE, minimum.value=5) #modification indices
print(modificationIndices[modificationIndices$op == "~",])
print(modificationIndices[modificationIndices$op == "~~",])
semPaths(modelFit, 'Standard', 'Estimates', layout='tree5', residuals=F, intercepts=F)




#model 2 - years 2 vs 5 - nothing matters
correSEMdata5 <- correSEMyears%>%
  filter(treatment_year!=1&treatment_year<5)%>%
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

correModel2Fit <- sem(correModel2, data=correSEMdata5, meanstructure=TRUE)
# adjusting for nested design
design <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata5)
model2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=design)
summary(model2Fit, standardize=T)
modificationIndices2 <- modindices(model2Fit, sort.=TRUE, minimum.value=5) #modification indices
print(modificationIndices2[modificationIndices2$op == "~",])
print(modificationIndices2[modificationIndices2$op == "~~",])
semPaths(model2Fit, 'Standardized', 'Estimates', layout='tree5')

correModel5 <- '
#year 5
anpp_change_transform_5 ~ mean_change_transform_5 + physiology_latent_5
mean_change_transform_5 ~ gain_5 + loss_5 + even_transform_5 + MRSc_5 + S_transform_5
physiology_latent_5 =~ treatment_binary
gain_5 ~ treatment_binary
loss_5 ~ treatment_binary
even_transform_5 ~ treatment_binary
MRSc_5 ~ treatment_binary
S_transform_5 ~ treatment_binary

#covariances
gain_5~~loss_5
gain_5~~even_transform_5
gain_5~~MRSc_5
gain_5~~S_transform_5
loss_5~~even_transform_5
loss_5~~MRSc_5
loss_5~~S_transform_5
even_transform_5~~MRSc_5
even_transform_5~~S_transform_5
MRSc_5~~S_transform_5'

correModel5Fit <- sem(correModel5, data=correSEMdata5, meanstructure=TRUE)
# adjusting for nested design
design <- svydesign(ids=~site_code, nest=TRUE, data=correSEMdata5)
model5Fit <- lavaan.survey(lavaan.fit=correModel5Fit, survey.design=design)
summary(model5Fit, standardize=T)
modificationIndices5 <- modindices(model5Fit, sort.=TRUE, minimum.value=5) #modification indices
print(modificationIndices5[modificationIndices5$op == "~",])
print(modificationIndices5[modificationIndices5$op == "~~",])
semPaths(model5Fit, 'Standardized', 'Estimates', layout='tree5')


#model 5 - years 2 vs 10
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
modificationIndices10 <- modindices(model10Fit, sort.=TRUE, minimum.value=5) #modification indices
print(modificationIndices10[modificationIndices10$op == "~",])
print(modificationIndices10[modificationIndices10$op == "~~",])
semPaths(model10Fit, 'Standardized', 'Estimates', layout='tree5')


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


##add treatment info to get binary treatments of various manipulation types
trt <- read.csv('ExperimentInformation_May2017.csv')%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(-site_code, -project_name, -community_type, -treatment_year, -public, -plot_mani)

#subset out anything without 10 years
numYears <- correSEMdata%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_years=length(treatment_year))

#subset out anything without years 1-10 of data
correSEMdataTrt10 <- correSEMdata%>%
  left_join(trt)%>%
  mutate(n_trt=ifelse(n>0, 1, 0), p_trt=ifelse(p>0, 1, 0), k_trt=ifelse(k>0, 1, 0), CO2_trt=ifelse(CO2>0, 1, 0), irr_trt=ifelse(precip>0, 1, 0), drought_trt=ifelse(precip<0, 1, 0), temp_trt=n_trt+p_trt+k_trt+CO2_trt+irr_trt, other_trt=ifelse((plot_mani-temp_trt)>0, 1, 0))%>%
  left_join(numYears)%>%
  filter(num_years>9&project_name!='TMECE'&project_name!='RaMPs')%>%
  mutate(test=1)

#write.csv(correSEMdataTrt10, 'SEM_10yr.csv')

#check the subset
ggplot(data=correSEMdataTrt, aes(x=treatment_year, y=test)) +
  geom_line() +
  facet_wrap(~site_project_comm)

#keep all data
correSEMdataTrt <- correSEMdata%>%
  left_join(trt)%>%
  mutate(n_trt=ifelse(n>0, 1, 0), p_trt=ifelse(p>0, 1, 0), k_trt=ifelse(k>0, 1, 0), CO2_trt=ifelse(CO2>0, 1, 0), irr_trt=ifelse(precip>0, 1, 0), drought_trt=ifelse(precip<0, 1, 0), temp_trt=n_trt+p_trt+k_trt+CO2_trt+irr_trt, other_trt=ifelse((plot_mani-temp_trt)>0, 1, 0))

#write.csv(correSEMdataTrt, 'SEM_allyr.csv')

###all data-----------
#put drought into "other_trt" because because long term experiments don't have drought
correModel <-  '
anpp_PC_transform ~ mean_change_transform + n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
mean_change_transform ~ PCSdiff + spdiffc + PCEvendiff_transform + MRSc_diff + n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
PCSdiff ~ n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
spdiffc ~ n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
PCEvendiff_transform ~ n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
MRSc_diff ~ n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt

#covariances
PCSdiff~~spdiffc
PCSdiff~~PCEvendiff_transform
PCSdiff~~MRSc_diff
spdiffc~~PCEvendiff_transform
spdiffc~~MRSc_diff
PCEvendiff_transform~~MRSc_diff
'
#year 1
correModelFit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==1), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==1))
survey1Fit <- lavaan.survey(lavaan.fit=correModelFit, survey.design=surveyDesign)
survey1Fit  #gives chi-square
# summary(survey1Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst1 <- parameterEstimates(survey1Fit, standardized=TRUE)%>%
  mutate(treatment_year=1)

#year 2
correModel2Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==2), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==2))
survey2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=surveyDesign)
survey2Fit  #gives chi-square
# summary(survey2Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst2 <- parameterEstimates(survey2Fit, standardized=TRUE)%>%
  mutate(treatment_year=2)

#year 3
correModel3Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==3), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==3))
survey3Fit <- lavaan.survey(lavaan.fit=correModel3Fit, survey.design=surveyDesign)
survey3Fit  #gives chi-square
# summary(survey3Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst3 <- parameterEstimates(survey3Fit, standardized=TRUE)%>%
  mutate(treatment_year=3)

#year 4
correModel4Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==4), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==4))
survey4Fit <- lavaan.survey(lavaan.fit=correModel4Fit, survey.design=surveyDesign)
survey4Fit  #gives chi-square
# summary(survey4Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst4 <- parameterEstimates(survey4Fit, standardized=TRUE)%>%
  mutate(treatment_year=4)

#year 5
correModel5Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==5), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==5))
survey5Fit <- lavaan.survey(lavaan.fit=correModel5Fit, survey.design=surveyDesign)
survey5Fit  #gives chi-square
# summary(survey5Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst5 <- parameterEstimates(survey5Fit, standardized=TRUE)%>%
  mutate(treatment_year=5)

#year 6
correModel6Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==6), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==6))
survey6Fit <- lavaan.survey(lavaan.fit=correModel6Fit, survey.design=surveyDesign)
survey6Fit  #gives chi-square
# summary(survey6Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst6 <- parameterEstimates(survey6Fit, standardized=TRUE)%>%
  mutate(treatment_year=6)

#year 7
correModel7Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==7), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==7))
survey7Fit <- lavaan.survey(lavaan.fit=correModel7Fit, survey.design=surveyDesign)
survey7Fit  #gives chi-square
# summary(survey7Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst7 <- parameterEstimates(survey7Fit, standardized=TRUE)%>%
  mutate(treatment_year=7)

#year 8
correModel8Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==8), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==8))
survey8Fit <- lavaan.survey(lavaan.fit=correModel8Fit, survey.design=surveyDesign)
survey8Fit  #gives chi-square
# summary(survey8Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst8 <- parameterEstimates(survey8Fit, standardized=TRUE)%>%
  mutate(treatment_year=8)

#year 9
correModel9Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==9), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==9))
survey9Fit <- lavaan.survey(lavaan.fit=correModel9Fit, survey.design=surveyDesign)
survey9Fit  #gives chi-square
# summary(survey9Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst9 <- parameterEstimates(survey9Fit, standardized=TRUE)%>%
  mutate(treatment_year=9)

#year 10
correModel10Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==10), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==10))
survey10Fit <- lavaan.survey(lavaan.fit=correModel10Fit, survey.design=surveyDesign)
survey10Fit  #gives chi-square
# summary(survey10Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst10 <- parameterEstimates(survey10Fit, standardized=TRUE)%>%
  mutate(treatment_year=10)

#bind together output from years 1-10
stdEst <- rbind(stdEst1, stdEst2, stdEst3, stdEst4, stdEst5, stdEst6, stdEst7, stdEst8, stdEst9, stdEst10)%>%
  na.omit()%>%
  mutate(std_alt=ifelse(pvalue>0.1, 0, std.all), std_weighted=std.all/se, est_alt=ifelse(pvalue>0.1, 0, est), est_weighted=est/se, trt_alt=ifelse(rhs=='CO2_trt'|rhs=='irr_trt'|rhs=='k_trt'|rhs=='n_trt'|rhs=='other_trt'|rhs=='p_trt', 'physiology', 'mean_change_transform'))

#all years figs-----------
#std.all includes non-sig effects
#std_alt sets non-sig effects to 0
#est is non-standardized effects
#est_alt sets non-sig non-standardized effects to 0
ggplot(data=subset(stdEst, lhs=='anpp_PC_transform'&rhs!=''&rhs!='anpp_PC_transform'), aes(x=treatment_year, y=est, color=rhs)) +
  geom_point(size=5) +
  geom_smooth(method='lm', size=3, se=F) +
  ylab('Effect Size') +
  xlab('Treatment Year')

ggplot(data=subset(stdEst, lhs=='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'), aes(x=treatment_year, y=std.all, color=rhs)) +
  geom_point(size=5) +
  geom_line(size=3) +
  ylab('Effect Size') +
  xlab('Treatment Year')

ggplot(data=subset(stdEst, lhs!='anpp_PC_transform'&lhs!='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'&rhs!='MRSc_diff'&rhs!='PCEvendiff_transform'&rhs!='spdiffc'&rhs!='PCSdiff'), aes(x=treatment_year, y=est_alt, color=rhs)) +
  geom_point(size=5) +
  geom_line(size=3) +
  facet_wrap(~lhs)



###only 10 year datasets-----------
#put drought into "other_trt" because because long term experiments don't have drought
correModel <-  '
anpp_PC_transform ~ mean_change_transform + n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
mean_change_transform ~ PCSdiff + spdiffc + PCEvendiff_transform + MRSc_diff + n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
PCSdiff ~ n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
spdiffc ~ n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
PCEvendiff_transform ~ n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt
MRSc_diff ~ n_trt + p_trt + k_trt + CO2_trt + irr_trt + other_trt

#covariances
PCSdiff~~spdiffc
PCSdiff~~PCEvendiff_transform
PCSdiff~~MRSc_diff
spdiffc~~PCEvendiff_transform
spdiffc~~MRSc_diff
PCEvendiff_transform~~MRSc_diff
'

#year 1
correModelFit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==1), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==1))
survey1Fit <- lavaan.survey(lavaan.fit=correModelFit, survey.design=surveyDesign)
survey1Fit  #gives chi-square
# summary(survey1Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst1 <- parameterEstimates(survey1Fit, standardized=TRUE)%>%
  mutate(treatment_year=1)

#year 2
correModel2Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==2), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==2))
survey2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=surveyDesign)
survey2Fit  #gives chi-square
# summary(survey2Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst2 <- parameterEstimates(survey2Fit, standardized=TRUE)%>%
  mutate(treatment_year=2)

#year 3
correModel3Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==3), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==3))
survey3Fit <- lavaan.survey(lavaan.fit=correModel3Fit, survey.design=surveyDesign)
survey3Fit  #gives chi-square
# summary(survey3Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst3 <- parameterEstimates(survey3Fit, standardized=TRUE)%>%
  mutate(treatment_year=3)

#year 4
correModel4Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==4), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==4))
survey4Fit <- lavaan.survey(lavaan.fit=correModel4Fit, survey.design=surveyDesign)
survey4Fit  #gives chi-square
# summary(survey4Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst4 <- parameterEstimates(survey4Fit, standardized=TRUE)%>%
  mutate(treatment_year=4)

#year 5
correModel5Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==5), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==5))
survey5Fit <- lavaan.survey(lavaan.fit=correModel5Fit, survey.design=surveyDesign)
survey5Fit  #gives chi-square
# summary(survey5Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst5 <- parameterEstimates(survey5Fit, standardized=TRUE)%>%
  mutate(treatment_year=5)

#year 6
correModel6Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==6), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==6))
survey6Fit <- lavaan.survey(lavaan.fit=correModel6Fit, survey.design=surveyDesign)
survey6Fit  #gives chi-square
# summary(survey6Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst6 <- parameterEstimates(survey6Fit, standardized=TRUE)%>%
  mutate(treatment_year=6)

#year 7
correModel7Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==7), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==7))
survey7Fit <- lavaan.survey(lavaan.fit=correModel7Fit, survey.design=surveyDesign)
survey7Fit  #gives chi-square
# summary(survey7Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst7 <- parameterEstimates(survey7Fit, standardized=TRUE)%>%
  mutate(treatment_year=7)

#year 8
correModel8Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==8), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==8))
survey8Fit <- lavaan.survey(lavaan.fit=correModel8Fit, survey.design=surveyDesign)
survey8Fit  #gives chi-square
# summary(survey8Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst8 <- parameterEstimates(survey8Fit, standardized=TRUE)%>%
  mutate(treatment_year=8)

#year 9
correModel9Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==9), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==9))
survey9Fit <- lavaan.survey(lavaan.fit=correModel9Fit, survey.design=surveyDesign)
survey9Fit  #gives chi-square
# summary(survey9Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst9 <- parameterEstimates(survey9Fit, standardized=TRUE)%>%
  mutate(treatment_year=9)

#year 10
correModel10Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==10), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==10))
survey10Fit <- lavaan.survey(lavaan.fit=correModel10Fit, survey.design=surveyDesign)
survey10Fit  #gives chi-square
# summary(survey10Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst10 <- parameterEstimates(survey10Fit, standardized=TRUE)%>%
  mutate(treatment_year=10)

#bind together output from years 1-10
stdEst <- rbind(stdEst1, stdEst2, stdEst3, stdEst4, stdEst5, stdEst6, stdEst7, stdEst8, stdEst9, stdEst10)%>%
  na.omit()%>%
  mutate(std_alt=ifelse(pvalue>0.05, 0, std.all))

#10 year figs-----------
ggplot(data=subset(stdEst, lhs=='anpp_PC_transform'&rhs!=''&rhs!='anpp_PC_transform'), aes(x=treatment_year, y=std_alt, color=rhs)) +
  geom_point(size=5) +
  geom_line(size=3)

ggplot(data=subset(stdEst, lhs=='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'), aes(x=treatment_year, y=std_alt, color=rhs)) +
  geom_point(size=5) +
  geom_line(size=3)

ggplot(data=subset(stdEst, lhs!='anpp_PC_transform'&lhs!='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'&rhs!='MRSc_diff'&rhs!='PCEvendiff_transform'&rhs!='spdiffc'&rhs!='PCSdiff'), aes(x=treatment_year, y=std_alt, color=rhs)) +
  geom_point(size=5) +
  geom_line(size=3) +
  facet_wrap(~lhs)
