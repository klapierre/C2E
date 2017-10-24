library(ggplot2)
library(grid)
library(mgcv)
library(lsmeans)
library(tidyverse)
library(nlme)

#kim's laptop
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's desktop
setwd("C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

###functions
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}

colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)

##################################################################################
##################################################################################
#experiment information --------------------------------------------------------
expRaw <- read.csv('ExperimentInformation_May2017.csv')

expInfo <- expRaw%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

rawData <- read.csv('ANPP_Spatail_ForAnalysis.csv')

rawData2<- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  summarise(anpp_mean=mean(anpp_sp_cv), std_mean=sd(anpp_sp_cv)) #to backtransform

#select just data in this analysis
expInfo2 <- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani))

#for table of experiment summarizing various factors
rawDataAll <- read.csv('ANPP_Spatail_ForAnalysis.csv')
expInfoSummary <- rawDataAll%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich), anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  summarise(length_mean=mean(experiment_length), length_min=min(experiment_length), length_max=max(experiment_length),
            plot_mani_median=mean(plot_mani), plot_mani_min=min(plot_mani), plot_mani_max=max(plot_mani),
            rrich_mean=mean(rrich), rrich_min=min(rrich), rrich_max=max(rrich),
            anpp_mean=mean(anpp), anpp_min=min(anpp), anpp_max=max(anpp),
            MAP_mean=mean(MAP), MAP_min=min(MAP), MAP_max=max(MAP),
            MAT_mean=mean(MAT), MAT_min=min(MAT), MAT_max=max(MAT))%>%
  gather(variable, estimate)

#treatment info
trtInfo <- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich),
            anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  left_join(expInfo)%>%
  mutate(resource_mani=(nutrients+carbon+irrigation+drought), id=1:length(treatment))%>%
  filter(id!=211&id!=212&id!=213&id!=214)

################################################################################
################################################################################

chainsSpatial <- read.csv('bayesian_anpp_spatial_output_mean sd_102117.csv')

#merge together with experiment list
chainsSpatialTrt <- chainsSpatial%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsSpatialTrt%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  #mutate(alt_length=ifelse(alt_length>=10, 9, alt_length))%>%
  mutate(yr10=(intercept+linear*9+quadratic*9^2)*(16.73752)+(29.78862))%>%
  mutate(yr_final=(intercept+linear*alt_length+quadratic*alt_length^2)*(16.73752)+(29.78862))
 


###by magnitude of resource manipulated---------------------------------
trtDetail <- expRaw%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))%>%
  ungroup()%>%
  mutate(site_code=as.character(site_code), project_name=as.character(project_name), community_type=as.character(community_type), treatment=as.character(treatment))

trtInteractions <- read.csv('treatment interactions_09072017.csv')%>%
  select(site_code, project_name, community_type, treatment, trt_type)

chainsEquationsTrts <- chainsEquations%>%
  select(-irrigation, -drought)%>%
  left_join(trtDetail)%>%
  select(-resource_mani)%>%
  left_join(trtInteractions)%>%
  mutate(trt_type=as.character(trt_type), trt_type_2=ifelse(plot_mani==0, 'control', trt_type))%>%
  mutate(trt_type_2=as.factor(trt_type_2))

nPlotFinal <- ggplot(data=subset(chainsEquationsTrts, n>0&plot_mani==1), aes(x=n, y=yr10)) +
  geom_point(size=5) +
  # scale_x_log10() +
  scale_y_continuous(name='ANPP Spatial Stability') +
  xlab('Nitrogen Amount (g/m2)') +
  geom_smooth(method='lm')

h2oPlotFinal <- ggplot(data=subset(chainsEquationsTrts, precip!=0&plot_mani==1), aes(x=precip, y=yr10)) +
  geom_point(size=5) +
  # scale_x_log10() +
  scale_y_continuous(name='ANPP Spatial Stability') +
  xlab('Precip Deviation (%))') +
  geom_smooth(method='lm')

pushViewport(viewport(layout=grid.layout(1,2)))
print(nPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(h2oPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 1800 x 800



###comparing different resource manipulation types (Figure 3)
controls <- chainsEquationsTrts%>%
  filter(trt_type_2=='control')%>%
  mutate(yr10ctl=yr10)%>%
  select(site_code, project_name, community_type, yr10ctl)

lnRR <- chainsEquationsTrts%>%
  filter(trt_type!='control')%>%
  left_join(controls)%>%
  mutate(lnRR=log(yr10/yr10ctl))

ggplot(barGraphStats(data=subset(lnRR, trt_type_2!='NA'), variable="lnRR", byFactorNames=c("trt_type_2")), aes(x=trt_type_2, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.2) +
  scale_y_continuous(name='ANPP Spatial Stability') +
  xlab('') +
  scale_x_discrete(limits=c('CO2', 'drought', 'irr', 'N', 'P', 'N+CO2', 'N+irr', 'N+P', 'N+P+K', 'N+P+K+irr', 'other'), labels=c(expression(paste(CO[2], '(2)')), 'Drought (4)', 'Irrigation (16)', 'N (38)', 'P (12)', expression(paste(CO[2],'*N (2)')), 'N*Irr (5)', 'N*P (14)', 'N*P*K (7)', 'N*P*K*Irr (1)', 'Other (101)')) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
