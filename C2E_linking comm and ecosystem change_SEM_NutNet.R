library(ggplot2)
library(grid)
library(PerformanceAnalytics)
library(lavaan)
library(lavaan.survey)
library(semPlot)
library(tidyverse)

#kim's desktop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


#######MODEL STRUCTURE: use differences between treatment and control plots in each year (horizontal arrows); include treatments as factors at the bottom, with multi-factor treatments having more than one greater than 0 (and make them categorical response variables)-------------------
####data input

###community data
nutnetCommChange <- read.csv('NutNet_community differences_12052017.csv')%>%
  mutate(treatment_year=time)%>%
  select(-X, -time)%>%
  mutate(n=ifelse(treatment=='N'|treatment=='NP'|treatment=='NK'|treatment=='NPK'|treatment=='NPK+Fence', 10, 0), p=ifelse(treatment=='P'|treatment=='NP'|treatment=='PK'|treatment=='NPK'|treatment=='NPK+Fence', 10, 0), k=ifelse(treatment=='K'|treatment=='NK'|treatment=='PK'|treatment=='NPK'|treatment=='NPK+Fence', 10, 0), fence=ifelse(treatment=='Fence'|treatment=='NPK+Fence', 1, 0))

#checking community data for outliers
ggplot(nutnetCommChange, aes(comp_diff)) + geom_histogram() + facet_wrap(~site_code, scales='free')

###anpp outliers
nutnetANPPoutliers <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\NutNet data\\La Pierre_NutNet_anpp_potential outliers_12122017.csv')%>%
  filter(checked.with.PI=='incorrect')%>%
  select(-live, -notes)

###anpp data
#remove outliers
nutnetANPP <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\NutNet data\\full-biomass-04-December-2017.csv')%>%
  filter(year_trt!=0, live==1)%>%
  merge(nutnetANPPoutliers, by=c("year", "year_trt", "trt", "site_name", "site_code", "block", "plot", "subplot", "mass", "category"), all.x=T)%>%
  filter(is.na(checked.with.PI))%>%
  group_by(site_code, plot, year_trt, trt)%>%
  summarise(anpp=sum(mass))


# #checking site level data for outliers in anpp
# ggplot(subset(nutnetANPP, category=='GRAMINOID'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='FORB'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='LEGUME'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='BRYOPHYTE'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='CACTUS'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='LIVE'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='PERENNIAL'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='WOODY'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='VASCULAR'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')


#calculating anpp difference
#anpp ctl data
nutnetANPPctl <- nutnetANPP%>%
  filter(trt=='Control')%>%
  select(site_code, year_trt, anpp)%>%
  group_by(site_code, year_trt)%>%
  summarise(anpp_ctl=mean(anpp))%>%
  ungroup()
#anpp change
nutnetANPPchange <- nutnetANPP%>%
  filter(trt!='Control')%>%
  group_by(site_code, year_trt, trt)%>%
  summarise(anpp=mean(anpp))%>%
  ungroup()%>%
  left_join(nutnetANPPctl)%>%
  #calculate anpp change as percent change from ctl in each year
  mutate(anpp_PC=(anpp-anpp_ctl)/anpp_ctl, treatment=trt, treatment_year=year_trt)%>%
  select(-anpp, -anpp_ctl, -trt, -year_trt)

#calculate dataset lengths
numYears <- nutnetANPPchange%>%
  select(treatment_year, site_code)%>%
  unique()%>%
  group_by(site_code)%>%
  summarise(num_years=length(treatment_year))

#merge
nutnetSEMdata <- nutnetANPPchange%>%
  left_join(nutnetCommChange)%>%
  na.omit()%>%
  #filter out experiments less than 3 years old
  left_join(numYears)%>%
  filter(num_years>2)%>%
  #drop fencing treatment
  filter(fence==0)%>%
  #transform to improve normality
  mutate(mean_change_transform=sqrt(comp_diff), Ed_transform=log10(as.numeric(Ed)), Sd_transform=log10(as.numeric(Sd)+1), anpp_PC_transform=log10(anpp_PC+1-min(anpp_PC)))
#all of these compare treatment to control!


###exploratory correlations and histograms (all variables compare treatment to control plots)
dataVis <- nutnetSEMdata%>%
  select(anpp_PC_transform, mean_change_transform, Sd_transform, Ed_transform, Rd, spd) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)


##anpp change through time by trt type
#all yrs
ggplot(data=nutnetSEMdata, aes(x=treatment_year, y=anpp_PC_transform)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=treatment)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('log ANPP (%) Difference')
#export at 1000 x 800



###anpp change correlations
#with mean_change
anpp_mean_y1 <- ggplot(data=subset(nutnetSEMdata, treatment_year==1), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,1.2) + xlim(0,1)
anpp_mean_y2 <- ggplot(data=subset(nutnetSEMdata, treatment_year==2), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,1.2) + xlim(0,1)
anpp_mean_y3 <- ggplot(data=subset(nutnetSEMdata, treatment_year==3), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,1.2) + xlim(0,1)
anpp_mean_y4 <- ggplot(data=subset(nutnetSEMdata, treatment_year==4), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,1.2) + xlim(0,1)
anpp_mean_y5 <- ggplot(data=subset(nutnetSEMdata, treatment_year==5), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,1.2) + xlim(0,1)
anpp_mean_y6 <- ggplot(data=subset(nutnetSEMdata, treatment_year==6), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,1.2) + xlim(0,1)
anpp_mean_y7 <- ggplot(data=subset(nutnetSEMdata, treatment_year==7), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,1.2) + xlim(0,1)
anpp_mean_y8 <- ggplot(data=subset(nutnetSEMdata, treatment_year==8), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,1.2) + xlim(0,1)
anpp_mean_y9 <- ggplot(data=subset(nutnetSEMdata, treatment_year==9), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,1.2) + xlim(0,1)
#with N
anpp_n_y1 <- ggplot(data=subset(nutnetSEMdata, treatment_year==1), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,1.2) + xlim(0,10)
anpp_n_y2 <- ggplot(data=subset(nutnetSEMdata, treatment_year==2), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,1.2) + xlim(0,10)
anpp_n_y3 <- ggplot(data=subset(nutnetSEMdata, treatment_year==3), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,1.2) + xlim(0,10)
anpp_n_y4 <- ggplot(data=subset(nutnetSEMdata, treatment_year==4), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,1.2) + xlim(0,10)
anpp_n_y5 <- ggplot(data=subset(nutnetSEMdata, treatment_year==5), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,1.2) + xlim(0,10)
anpp_n_y6 <- ggplot(data=subset(nutnetSEMdata, treatment_year==6), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,1.2) + xlim(0,10)
anpp_n_y7 <- ggplot(data=subset(nutnetSEMdata, treatment_year==7), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,1.2) + xlim(0,10)
anpp_n_y8 <- ggplot(data=subset(nutnetSEMdata, treatment_year==8), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,1.2) + xlim(0,10)
anpp_n_y9 <- ggplot(data=subset(nutnetSEMdata, treatment_year==9), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,1.2) + xlim(0,10)

pushViewport(viewport(layout=grid.layout(9,2)))
print(anpp_mean_y1, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(anpp_mean_y2, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(anpp_mean_y3, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(anpp_mean_y4, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(anpp_mean_y5, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))
print(anpp_mean_y6, vp=viewport(layout.pos.row = 6, layout.pos.col = 1))
print(anpp_mean_y7, vp=viewport(layout.pos.row = 7, layout.pos.col = 1))
print(anpp_mean_y8, vp=viewport(layout.pos.row = 8, layout.pos.col = 1))
print(anpp_mean_y9, vp=viewport(layout.pos.row = 9, layout.pos.col = 1))
print(anpp_n_y1, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(anpp_n_y2, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(anpp_n_y3, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(anpp_n_y4, vp=viewport(layout.pos.row = 4, layout.pos.col = 2))
print(anpp_n_y5, vp=viewport(layout.pos.row = 5, layout.pos.col = 2))
print(anpp_n_y6, vp=viewport(layout.pos.row = 6, layout.pos.col = 2))
print(anpp_n_y7, vp=viewport(layout.pos.row = 7, layout.pos.col = 2))
print(anpp_n_y8, vp=viewport(layout.pos.row = 8, layout.pos.col = 2))
print(anpp_n_y9, vp=viewport(layout.pos.row = 9, layout.pos.col = 2))
#export at 2000 x 4000


###Community Change correlations
#with mean_change
mean_Rd_y1 <- ggplot(data=subset(nutnetSEMdata, treatment_year==1), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y2 <- ggplot(data=subset(nutnetSEMdata, treatment_year==2), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y3 <- ggplot(data=subset(nutnetSEMdata, treatment_year==3), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y4 <- ggplot(data=subset(nutnetSEMdata, treatment_year==4), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y5 <- ggplot(data=subset(nutnetSEMdata, treatment_year==5), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y6 <- ggplot(data=subset(nutnetSEMdata, treatment_year==6), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y7 <- ggplot(data=subset(nutnetSEMdata, treatment_year==7), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y8 <- ggplot(data=subset(nutnetSEMdata, treatment_year==8), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y9 <- ggplot(data=subset(nutnetSEMdata, treatment_year==9), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
#with N
mean_n_y1 <- ggplot(data=subset(nutnetSEMdata, treatment_year==1), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,10)
mean_n_y2 <- ggplot(data=subset(nutnetSEMdata, treatment_year==2), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,10)
mean_n_y3 <- ggplot(data=subset(nutnetSEMdata, treatment_year==3), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,10)
mean_n_y4 <- ggplot(data=subset(nutnetSEMdata, treatment_year==4), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,10)
mean_n_y5 <- ggplot(data=subset(nutnetSEMdata, treatment_year==5), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,10)
mean_n_y6 <- ggplot(data=subset(nutnetSEMdata, treatment_year==6), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,10)
mean_n_y7 <- ggplot(data=subset(nutnetSEMdata, treatment_year==7), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,10)
mean_n_y8 <- ggplot(data=subset(nutnetSEMdata, treatment_year==8), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,10)
mean_n_y9 <- ggplot(data=subset(nutnetSEMdata, treatment_year==9), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,10)

pushViewport(viewport(layout=grid.layout(9,2)))
print(mean_Rd_y1, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(mean_Rd_y2, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(mean_Rd_y3, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(mean_Rd_y4, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(mean_Rd_y5, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))
print(mean_Rd_y6, vp=viewport(layout.pos.row = 6, layout.pos.col = 1))
print(mean_Rd_y7, vp=viewport(layout.pos.row = 7, layout.pos.col = 1))
print(mean_Rd_y8, vp=viewport(layout.pos.row = 8, layout.pos.col = 1))
print(mean_Rd_y9, vp=viewport(layout.pos.row = 9, layout.pos.col = 1))
print(mean_n_y1, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(mean_n_y2, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(mean_n_y3, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(mean_n_y4, vp=viewport(layout.pos.row = 4, layout.pos.col = 2))
print(mean_n_y5, vp=viewport(layout.pos.row = 5, layout.pos.col = 2))
print(mean_n_y6, vp=viewport(layout.pos.row = 6, layout.pos.col = 2))
print(mean_n_y7, vp=viewport(layout.pos.row = 7, layout.pos.col = 2))
print(mean_n_y8, vp=viewport(layout.pos.row = 8, layout.pos.col = 2))
print(mean_n_y9, vp=viewport(layout.pos.row = 9, layout.pos.col = 2))
#export at 2000 x 4000

###all data-----------
nutnetModel <-  '
anpp_PC_transform ~ mean_change_transform + n + p + k
mean_change_transform ~ Sd_transform + Ed_transform + Rd + spd + n + p + k
Sd_transform ~ n + p + k
Ed_transform ~ n + p + k
Rd ~ n + p + k
spd ~ n + p + k

#covariances
Sd_transform~~Ed_transform
Sd_transform~~Rd
Sd_transform~~spd
Ed_transform~~Rd
Ed_transform~~spd
Rd~~spd
'
#year 1
nutnetModelFit <- sem(nutnetModel, data=subset(nutnetSEMdata, treatment_year==1), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=subset(nutnetSEMdata, treatment_year==1))
survey1Fit <- lavaan.survey(lavaan.fit=nutnetModelFit, survey.design=surveyDesign)
survey1Fit  #gives chi-square
# summary(survey1Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst1 <- parameterEstimates(survey1Fit, standardized=TRUE)%>%
  mutate(treatment_year=1)

#year 2
nutnetModel2Fit <- sem(nutnetModel, data=subset(nutnetSEMdata, treatment_year==2), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=subset(nutnetSEMdata, treatment_year==2))
survey2Fit <- lavaan.survey(lavaan.fit=nutnetModel2Fit, survey.design=surveyDesign)
survey2Fit  #gives chi-square
# summary(survey2Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst2 <- parameterEstimates(survey2Fit, standardized=TRUE)%>%
  mutate(treatment_year=2)

#year 3
nutnetModel3Fit <- sem(nutnetModel, data=subset(nutnetSEMdata, treatment_year==3), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=subset(nutnetSEMdata, treatment_year==3))
survey3Fit <- lavaan.survey(lavaan.fit=nutnetModel3Fit, survey.design=surveyDesign)
survey3Fit  #gives chi-square
# summary(survey3Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst3 <- parameterEstimates(survey3Fit, standardized=TRUE)%>%
  mutate(treatment_year=3)

#year 4
nutnetModel4Fit <- sem(nutnetModel, data=subset(nutnetSEMdata, treatment_year==4), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=subset(nutnetSEMdata, treatment_year==4))
survey4Fit <- lavaan.survey(lavaan.fit=nutnetModel4Fit, survey.design=surveyDesign)
survey4Fit  #gives chi-square
# summary(survey4Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst4 <- parameterEstimates(survey4Fit, standardized=TRUE)%>%
  mutate(treatment_year=4)

#year 5
nutnetModel5Fit <- sem(nutnetModel, data=subset(nutnetSEMdata, treatment_year==5), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=subset(nutnetSEMdata, treatment_year==5))
survey5Fit <- lavaan.survey(lavaan.fit=nutnetModel5Fit, survey.design=surveyDesign)
survey5Fit  #gives chi-square
# summary(survey5Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst5 <- parameterEstimates(survey5Fit, standardized=TRUE)%>%
  mutate(treatment_year=5)

#year 6
nutnetModel6Fit <- sem(nutnetModel, data=subset(nutnetSEMdata, treatment_year==6), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=subset(nutnetSEMdata, treatment_year==6))
survey6Fit <- lavaan.survey(lavaan.fit=nutnetModel6Fit, survey.design=surveyDesign)
survey6Fit  #gives chi-square
# summary(survey6Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst6 <- parameterEstimates(survey6Fit, standardized=TRUE)%>%
  mutate(treatment_year=6)

#year 7
nutnetModel7Fit <- sem(nutnetModel, data=subset(nutnetSEMdata, treatment_year==7), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=subset(nutnetSEMdata, treatment_year==7))
survey7Fit <- lavaan.survey(lavaan.fit=nutnetModel7Fit, survey.design=surveyDesign)
survey7Fit  #gives chi-square
# summary(survey7Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst7 <- parameterEstimates(survey7Fit, standardized=TRUE)%>%
  mutate(treatment_year=7)

#year 8
nutnetModel8Fit <- sem(nutnetModel, data=subset(nutnetSEMdata, treatment_year==8), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=subset(nutnetSEMdata, treatment_year==8))
survey8Fit <- lavaan.survey(lavaan.fit=nutnetModel8Fit, survey.design=surveyDesign)
survey8Fit  #gives chi-square
# summary(survey8Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst8 <- parameterEstimates(survey8Fit, standardized=TRUE)%>%
  mutate(treatment_year=8)

#year 9
nutnetModel9Fit <- sem(nutnetModel, data=subset(nutnetSEMdata, treatment_year==9), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_code, nest=TRUE, data=subset(nutnetSEMdata, treatment_year==9))
survey9Fit <- lavaan.survey(lavaan.fit=nutnetModel9Fit, survey.design=surveyDesign)
survey9Fit  #gives chi-square
# summary(survey9Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst9 <- parameterEstimates(survey9Fit, standardized=TRUE)%>%
  mutate(treatment_year=9)

#bind together output from years 1-10
stdEst <- rbind(stdEst1, stdEst2, stdEst3, stdEst4, stdEst5, stdEst6, stdEst7, stdEst8, stdEst9)%>%
  na.omit()%>%
  mutate(std_alt=ifelse(pvalue>0.1, 0, std.all), std_weighted=std.all/se, est_alt=ifelse(pvalue>0.1, 0, est), est_weighted=est/se, trt_alt=ifelse(rhs=='n'|rhs=='p'|rhs=='k', 'physiology', 'mean_change_transform'))

#all years figs-----------
#std.all includes non-sig effects
#std_alt sets non-sig effects to 0
#est is non-standardized effects
#est_alt sets non-sig non-standardized effects to 0
ggplot(data=subset(stdEst, lhs=='anpp_PC_transform'&rhs!=''&rhs!='anpp_PC_transform'), aes(x=treatment_year, y=std.all)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se), width=0.5) +
  geom_line(size=2) +
  # geom_smooth(method='lm', size=3, se=F) +
  ylab('ANPP Difference Effect Size') +
  xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~rhs) +
  theme(strip.text=element_text(size=24))

ggplot(data=subset(stdEst, lhs=='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'), aes(x=treatment_year, y=std.all)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se), width=0.2) +
  geom_line(size=2) +
  ylab('Community Difference Effect Size') +
  xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~rhs) +
  theme(strip.text=element_text(size=24))

ggplot(data=subset(stdEst, lhs!='anpp_PC_transform'&lhs!='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'&rhs!='Sd_transform'&rhs!='Ed_transform'&rhs!='Rd'&rhs!='spd'), aes(x=treatment_year, y=std.all, color=rhs)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se)) +
  geom_line(size=2) +
  ylab('Effect Size') +
  xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~lhs) +
  theme(strip.text=element_text(size=24))

#bin into first and last 5 years for bar graphs
stdEstBin <- stdEst%>%
  mutate(temporal_bin=ifelse(treatment_year>5, '1-5', '6-10'))%>%
  group_by(lhs, rhs, temporal_bin)%>%
  summarise(est_mean=mean(est), est_sd=sd(est), est_N=n(), std.all_mean=mean(std.all), std.all_sd=sd(std.all), std.all_N=n())%>%
  ungroup()

ggplot(data=subset(stdEstBin, lhs=='anpp_PC_transform'&rhs!=''&rhs!='anpp_PC_transform'), aes(x=temporal_bin, y=std.all_mean, color=rhs)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=std.all_mean-(std.all_sd/sqrt(std.all_N)), ymax=std.all_mean+(std.all_sd/sqrt(std.all_N))), position=position_dodge(0.9), width=0.2) +
  ylab('Effect Size') +
  xlab('Treatment Year')

ggplot(data=subset(stdEstBin, lhs=='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'), aes(x=temporal_bin, y=std.all_mean, color=rhs)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=std.all_mean-(std.all_sd/sqrt(std.all_N)), ymax=std.all_mean+(std.all_sd/sqrt(std.all_N))), position=position_dodge(0.9), width=0.2) +
  geom_line(size=3) +
  ylab('Effect Size') +
  xlab('Treatment Year')

ggplot(data=subset(stdEstBin, lhs!='anpp_PC_transform'&lhs!='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'&rhs!='Sd_transform'&rhs!='Ed_transform'&rhs!='Rd'&rhs!='spd'), aes(x=temporal_bin, y=std.all_mean, color=rhs)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=std.all_mean-(std.all_sd/sqrt(std.all_N)), ymax=std.all_mean+(std.all_sd/sqrt(std.all_N))), position=position_dodge(0.9), width=0.2) +
  geom_line(size=3) +
  facet_wrap(~lhs)
