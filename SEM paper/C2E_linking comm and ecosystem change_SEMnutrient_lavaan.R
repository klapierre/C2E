library(grid)
library(PerformanceAnalytics)
library(lavaan)
library(lavaan.survey)
library(semPlot)
library(tidyverse)



#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###merge CoRRE and NutNet datasets
correSEMdataTrt <- read.csv('CoRRE_comm and anpp diff_07160219.csv')%>%
  filter(trt_type %in% c('N','P','N*P','mult_nutrient'))%>%
  select(site_project_comm, treatment_year, treatment, n, p, k, anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference)%>%
  mutate(dataset='CoRRE')

allSEMdata <- read.csv('NutNet_comm and anpp diff_07160219.csv')%>%
  filter(! treatment %in% c('Fence', 'NPK+Fence'))%>%
  rename(site_project_comm=site_code)%>%
  select(site_project_comm, treatment_year, treatment, n, p, k, anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference)%>%
  mutate(n=ifelse(treatment %in% c('N', 'NP', 'NK', 'NPK'), 10, 0), p=ifelse(treatment %in% c('P', 'NP', 'PK', 'NPK'), 10, 0), k=ifelse(treatment %in% c('K', 'NK', 'PK', 'NPK'), 10, 0))%>%
  mutate(dataset='NutNet')%>%
  rbind(correSEMdataTrt)%>%
  separate(col=site_project_comm, into=c('site_code', 'project_name', 'community_type'), sep='::', remove=F)%>%
  mutate(anpp_pdiff_transform=log(anpp_pdiff+(1-min(anpp_pdiff))), composition_diff_transform=log(composition_diff))

dataVis <- (allSEMdata)%>%
  select(anpp_pdiff_transform, composition_diff_transform, richness_difference, evenness_diff, rank_difference, species_difference) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)
  


#######MODEL STRUCTURE: use differences between treatment and control plots in each year (horizontal arrows)-------------------

###all data, all years-----------

###difference metrics not through composition model
correModel <-  '
level: 1
anpp_pdiff_transform ~ composition_diff_transform + n + p + k
composition_diff_transform ~ n + p + k + richness_difference + evenness_diff + rank_difference + species_difference
richness_difference ~ n + p + k
evenness_diff ~ n + p + k
rank_difference ~ n + p + k
species_difference ~ n + p + k

level: 2
richness_difference ~~ evenness_diff + rank_difference + species_difference
evenness_diff ~~ rank_difference + species_difference
rank_difference ~~ species_difference
'


#year 1
summary(correModelFit <- sem(correModel, data=subset(allSEMdata, treatment_year==1), meanstructure=TRUE, cluster='site_code', optim.method = "em"))
stdEst1 <- parameterEstimates(correModelFit, standardized=TRUE)%>%
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
  mutate(std_alt=ifelse(pvalue>0.1, 0, std.all), std_weighted=std.all/se, est_alt=ifelse(pvalue>0.1, 0, est), est_weighted=est/se, trt_alt=ifelse(rhs=='CO2'|rhs=='precip'|rhs=='k'|rhs=='n'|rhs=='other_trt'|rhs=='p', 'physiology', 'composition_diff'))

#all years figs-----------
#std.all includes non-sig effects
#std_alt sets non-sig effects to 0
#est is non-standardized effects
#est_alt sets non-sig non-standardized effects to 0
ggplot(data=subset(stdEst, lhs=='anpp_pdiff'&rhs!=''&rhs!='anpp_pdiff'), aes(x=treatment_year, y=std.all)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se), width=0.5) +
  geom_line(size=2) +
  # geom_smooth(method='lm', size=3, se=F) +
  ylab('ANPP (%) Difference Effect Size') +
  xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~rhs) +
  theme(strip.text=element_text(size=24))

ggplot(data=subset(stdEst, lhs=='composition_diff'&rhs!=''&rhs!='anpp_pdiff'&rhs!='composition_diff'), aes(x=treatment_year, y=std.all)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se), width=0.2) +
  geom_line(size=2) +
  ylab('Composition Difference Effect Size') +
  xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~rhs) +
  theme(strip.text=element_text(size=24))

ggplot(data=subset(stdEst, lhs!='anpp_pdiff'&lhs!='composition_diff'&rhs!=''&rhs!='anpp_pdiff'&rhs!='composition_diff'&rhs!='richness_difference'&rhs!='evenness_diff'&rhs!='rank_difference'&rhs!='species_difference'), aes(x=treatment_year, y=std.all, color=rhs)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se)) +
  geom_line(size=2) +
  ylab('Effect Size') +
  xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~lhs) +
  theme(strip.text=element_text(size=24))
