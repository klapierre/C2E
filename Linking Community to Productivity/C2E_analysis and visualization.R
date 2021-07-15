################################################################################
##  C2E_analysis and visualization.R: Data anaylsis and visualization.
##
##  Author: Kimberly Komatsu
##  Date created: July 12, 2021
################################################################################

library(lme4)
library(emmeans)
library(performance)
library(tidyverse)
library(ggpubr)

#options, themes, and functions

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\C2E\\Products\\testing HRF\\data')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=20),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=20),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=16), legend.text=element_text(size=16))

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


options(contrasts=c('contr.sum','contr.poly'))

#not in function
`%!in%` = Negate(`%in%`)


##### import data #####
correAllDataTrt <- read.csv('CoRRE_comm and anpp diff_07122021.csv')%>%
  select(site_project_comm, site_project_comm_trt, site_code, project_name, community_type, treatment_year, calendar_year, treatment, richness_diff, evenness_diff, rank_diff, species_diff, composition_diff, abs_dispersion_diff, richness_change, evenness_change, rank_change, gains, losses, composition_change, dispersion_change, anpp_pdiff, plot_mani)%>%
  drop_na()%>%
  filter(anpp_pdiff<4)

# ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)
# AIC greater than 6 to prefer polynomical model, due to sample size of 193 (trt within site_project_comm)

##### ANPP differences through time #####
summary(anppTemporal <- lmer(abs(anpp_pdiff)~poly(treatment_year,2) + (1+treatment_year|site_project_comm/treatment), data=correAllDataTrt))  #significant quadratic effect
check_model(anppTemporal)
r.squaredGLMM(anppTemporal) #marginal R2=0.008211761, conditional R2=0.8174605 
AIC(anppTemporal) #polynomial model (AIC=799.3321) better than liner model (AIC=815.9347)

ggplot(data=correAllDataTrt, aes(x=treatment_year, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Treatment Year') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

##### ANPP difference related to community difference #####
###overall community change
summary(anppCommChange <- lmer(abs(anpp_pdiff)~composition_change + (1+treatment_year|site_project_comm/treatment), data=correAllDataTrt)) #significant effect
check_model(anppCommChange)
r.squaredGLMM(anppCommChange) #marginal R2=0.02663398, conditional R2=0.7998936 
AIC(anppCommChange) #polynomial model (AIC=785.5084) NOT better than liner model (AIC=788.8264)

ggplot(data=correAllDataTrt, aes(x=composition_change, y=abs(anpp_pdiff))) +
    geom_point(color='grey') +
    geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
    geom_smooth(method='lm', se=F, size=3, color='black') +
    xlab('Composition Difference') + ylab('|ANPP % Difference|') +
    theme(legend.position='none')

###richness only
summary(anppRichChange <- lmer(abs(anpp_pdiff)~poly(richness_change,2) + (1+treatment_year|site_project_comm/treatment), data=correAllDataTrt)) #significant effect
check_model(anppRichChange)
r.squaredGLMM(anppRichChange) #marginal R2=0.06465337, conditional R2=0.8239859 
AIC(anppRichChange) #polynomial model (AIC=734.9237) better than liner model (AIC=751.6983)

ggplot(data=correAllDataTrt, aes(x=richness_change, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Richness Change') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

###multiple change metrics
#individually to determine quadratic vs linear
AIC(anppEvenChange <- lmer(abs(anpp_pdiff)~evenness_change + (1+treatment_year|site_project_comm/treatment), data=correAllDataTrt)) #polynomial model (AIC=809.2267) NOT better than liner model (AIC=809.5304)
AIC(anppRankChange <- lmer(abs(anpp_pdiff)~rank_change + (1+treatment_year|site_project_comm/treatment), data=correAllDataTrt)) #polynomial model (AIC=787.1062) NOT better than liner model (AIC=794.8822)
AIC(anppGainChange <- lmer(abs(anpp_pdiff)~gains + (1+treatment_year|site_project_comm/treatment), data=correAllDataTrt)) #polynomial model (AIC=793.1192) better than liner model (AIC=809.5304)
AIC(anppEvenChange <- lmer(abs(anpp_pdiff)~poly(losses,2) + (1+treatment_year|site_project_comm/treatment), data=correAllDataTrt)) #polynomial model (AIC=714.3914) better than liner model (AIC=725.6287)

#all together
summary(anppCommMetricChange <- lmer(abs(anpp_pdiff)~evenness_change+rank_change+poly(gains,2)+poly(losses,2) + (1+treatment_year|site_project_comm/treatment), data=correAllDataTrt)) #significant effect
check_model(anppCommMetricChange)
r.squaredGLMM(anppCommMetricChange) #marginal R2=0.06465337, conditional R2=0.8239859 
AIC(anppCommMetricChange) #polynomial model (AIC=734.9237) better than liner model (AIC=751.6983)

ggplot(data=correAllDataTrt, aes(x=richness_change, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Richness Change') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))